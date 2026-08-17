# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with the
# License. You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Methods for merging existing features to create a new example.

Covers:
- Merging features across chains.
- Merging the paired and unpaired parts of the MSA.
"""

from typing import TypeAlias

from alphafold3.model import data_constants
import jax.numpy as jnp
import numpy as np

NUM_SEQ_NUM_RES_MSA_FEATURES = data_constants.NUM_SEQ_NUM_RES_MSA_FEATURES
NUM_SEQ_MSA_FEATURES = data_constants.NUM_SEQ_MSA_FEATURES
MSA_PAD_VALUES = data_constants.MSA_PAD_VALUES


xnp_ndarray: TypeAlias = np.ndarray | jnp.ndarray  # pylint: disable=invalid-name
BatchDict: TypeAlias = dict[str, xnp_ndarray]


def _pad_features_to_max(feat_name: str, chains: list[BatchDict], axis: int):
  """Pad a set of features to the maximum size amongst all chains.

  Args:
    feat_name: The feature name to pad.
    chains: A list of chains with associated features.
    axis: Which axis to pad to the max.

  Returns:
    A list of features, all with the same size on the given axis.
  """
  max_num_seq = np.max([chain[feat_name].shape[axis] for chain in chains])

  padded_feats = []
  for chain in chains:
    feat = chain[feat_name]

    padding = np.zeros_like(feat.shape)  # pytype: disable=attribute-error
    padding[axis] = max_num_seq - feat.shape[axis]  # pytype: disable=attribute-error
    padding = [(0, p) for p in padding]
    padded_feats.append(
        np.pad(
            feat,
            padding,
            mode='constant',
            constant_values=MSA_PAD_VALUES[feat_name],
        )
    )
  return padded_feats


def merge_msa_features(feat_name: str, chains: list[BatchDict]) -> np.ndarray:
  """Merges MSA features with shape (NUM_SEQ, NUM_RES) across chains."""
  expected_dtype = chains[0][feat_name].dtype
  if '_all_seq' in feat_name:
    return np.concatenate(
        [c.get(feat_name, np.array([], expected_dtype)) for c in chains], axis=1
    )
  else:
    # Since each MSA can be of different lengths, we first need to pad them
    # all to the size of the largest MSA before concatenating.
    padded_feats = _pad_features_to_max(feat_name, chains, axis=0)
    return np.concatenate(padded_feats, axis=1)


def merge_paired_and_unpaired_msa(example: BatchDict) -> BatchDict:
  """Concatenates the paired (all_seq) MSA features with the unpaired ones."""
  new_example = dict(example)

  for feature_name in NUM_SEQ_NUM_RES_MSA_FEATURES + NUM_SEQ_MSA_FEATURES:
    if feature_name in example and feature_name + '_all_seq' in example:
      feat = example[feature_name]
      feat_all_seq = example[feature_name + '_all_seq']
      merged_feat = np.concatenate([feat_all_seq, feat], axis=0)
      new_example[feature_name] = merged_feat

  new_example['num_alignments'] = np.array(
      new_example['msa'].shape[0], dtype=np.int32
  )
  return new_example
