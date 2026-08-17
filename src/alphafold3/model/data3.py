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

"""Protein features that are computed from parsed mmCIF objects."""

from collections.abc import Mapping
import datetime
from typing import TypeAlias

from alphafold3.constants import residue_names
from alphafold3.cpp import msa_profile
from alphafold3.model import protein_data_processing
import numpy as np


FeatureDict: TypeAlias = Mapping[str, np.ndarray]


def get_profile_features(
    msa: np.ndarray, deletion_matrix: np.ndarray
) -> FeatureDict:
  """Returns the MSA profile and deletion_mean features."""
  num_restypes = residue_names.POLYMER_TYPES_NUM_WITH_UNKNOWN_AND_GAP
  profile = msa_profile.compute_msa_profile(
      msa=msa, num_residue_types=num_restypes
  )

  return {
      'profile': profile.astype(np.float32),
      'deletion_mean': np.mean(deletion_matrix, axis=0),
  }


def fix_template_features(
    template_features: FeatureDict, num_res: int
) -> FeatureDict:
  """Convert template features to AlphaFold 3 format.

  Args:
    template_features: Template features for the protein.
    num_res: The length of the amino acid sequence of the protein.

  Returns:
    Updated template_features for the chain.
  """
  if not template_features['template_aatype'].shape[0]:
    template_features = empty_template_features(num_res)
  else:
    template_release_timestamp = [
        _get_timestamp(x.decode('utf-8'))
        for x in template_features['template_release_date']
    ]

    # Convert from atom37 to dense atom
    dense_atom_indices = np.take(
        protein_data_processing.PROTEIN_AATYPE_DENSE_ATOM_TO_ATOM37,
        template_features['template_aatype'],
        axis=0,
    )

    atom_mask = np.take_along_axis(
        template_features['template_all_atom_masks'], dense_atom_indices, axis=2
    )
    atom_positions = np.take_along_axis(
        template_features['template_all_atom_positions'],
        dense_atom_indices[..., None],
        axis=2,
    )
    atom_positions *= atom_mask[..., None]

    template_features = {
        'template_aatype': template_features['template_aatype'],
        'template_atom_mask': atom_mask.astype(np.int32),
        'template_atom_positions': atom_positions.astype(np.float32),
        'template_domain_names': np.array(
            template_features['template_domain_names'], dtype=object
        ),
        'template_release_timestamp': np.array(
            template_release_timestamp, dtype=np.float32
        ),
    }
  return template_features


def empty_template_features(num_res: int) -> FeatureDict:
  """Creates a fully masked out template features to allow padding to work.

  Args:
    num_res: The length of the target chain.

  Returns:
    Empty template features for the chain.
  """
  template_features = {
      'template_aatype': np.zeros(num_res, dtype=np.int32)[None, ...],
      'template_atom_mask': np.zeros(
          (num_res, protein_data_processing.NUM_DENSE), dtype=np.int32
      )[None, ...],
      'template_atom_positions': np.zeros(
          (num_res, protein_data_processing.NUM_DENSE, 3), dtype=np.float32
      )[None, ...],
      'template_domain_names': np.array([b''], dtype=object),
      'template_release_timestamp': np.array([0.0], dtype=np.float32),
  }
  return template_features


def _get_timestamp(date_str: str):
  dt = datetime.datetime.fromisoformat(date_str)
  dt = dt.replace(tzinfo=datetime.timezone.utc)
  return dt.timestamp()
