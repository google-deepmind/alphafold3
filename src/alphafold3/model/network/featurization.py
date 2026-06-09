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

"""Model-side of the input features processing."""

import functools

from alphafold3.constants import residue_names
from alphafold3.model import feat_batch
from alphafold3.model import features
from alphafold3.model.components import utils
import jax
import jax.numpy as jnp


def _grid_keys(key, shape):
  """Generate a grid of rng keys that is consistent with different padding.

  Generate random keys such that the keys will be identical, regardless of
  how much padding is added to any dimension.

  Args:
    key: A PRNG key.
    shape: The shape of the output array of keys that will be generated.

  Returns:
    An array of shape `shape` consisting of random keys.
  """
  if not shape:
    return key
  new_keys = jax.vmap(functools.partial(jax.random.fold_in, key))(
      jnp.arange(shape[0])
  )
  return jax.vmap(functools.partial(_grid_keys, shape=shape[1:]))(new_keys)


def _padding_consistent_rng(f):
  """Modify any element-wise random function to be consistent with padding.

  Normally if you take a function like jax.random.normal and generate an array,
  say of size (10,10), you will get a different set of random numbers to if you
  add padding and take the first (10,10) sub-array.

  This function makes a random function that is consistent regardless of the
  amount of padding added.

  Note: The padding-consistent function is likely to be slower to compile and
  run than the function it is wrapping, but these slowdowns are likely to be
  negligible in a large network.

  Args:
    f: Any element-wise function that takes (PRNG key, shape) as the first 2
      arguments.

  Returns:
    An equivalent function to f, that is now consistent for different amounts of
    padding.
  """

  def inner(key, shape, **kwargs):
    keys = _grid_keys(key, shape)
    signature = (
        '()->()'
        if jax.dtypes.issubdtype(keys.dtype, jax.dtypes.prng_key)
        else '(2)->()'
    )
    return jnp.vectorize(
        functools.partial(f, shape=(), **kwargs), signature=signature
    )(keys)

  return inner


def gumbel_argsort_sample_idx(
    key: jnp.ndarray, logits: jnp.ndarray
) -> jnp.ndarray:
  """Samples with replacement from a distribution given by 'logits'.

  This uses Gumbel trick to implement the sampling an efficient manner. For a
  distribution over k items this samples k times without replacement, so this
  is effectively sampling a random permutation with probabilities over the
  permutations derived from the logprobs.

  Args:
    key: prng key
    logits: logarithm of probabilities to sample from, probabilities can be
      unnormalized.

  Returns:
    Sample from logprobs in one-hot form.
  """
  gumbel = _padding_consistent_rng(jax.random.gumbel)
  z = gumbel(key, logits.shape)
  # This construction is equivalent to jnp.argsort, but using a non stable sort,
  # since stable sort's aren't supported by jax2tf
  axis = len(logits.shape) - 1
  iota = jax.lax.broadcasted_iota(jnp.int64, logits.shape, axis)
  _, perm = jax.lax.sort_key_val(
      logits + z, iota, dimension=-1, is_stable=False
  )
  return perm[::-1]


def create_msa_feat(msa: features.MSA) -> jax.Array:
  """Create and concatenate MSA features."""
  msa_1hot = jax.nn.one_hot(
      msa.rows, residue_names.POLYMER_TYPES_NUM_WITH_UNKNOWN_AND_GAP + 1
  )
  deletion_matrix = msa.deletion_matrix
  has_deletion = jnp.clip(deletion_matrix, 0.0, 1.0)[..., None]
  deletion_value = (jnp.arctan(deletion_matrix / 3.0) * (2.0 / jnp.pi))[
      ..., None
  ]

  msa_feat = [
      msa_1hot,
      has_deletion,
      deletion_value,
  ]

  return jnp.concatenate(msa_feat, axis=-1)


def truncate_msa_batch(msa: features.MSA, num_msa: int) -> features.MSA:
  indices = jnp.arange(num_msa)
  return msa.index_msa_rows(indices)


def create_target_feat(
    batch: feat_batch.Batch,
    append_per_atom_features: bool,
) -> jax.Array:
  """Make target feat."""
  token_features = batch.token_features
  target_features = []
  target_features.append(
      jax.nn.one_hot(
          token_features.aatype,
          residue_names.POLYMER_TYPES_NUM_WITH_UNKNOWN_AND_GAP,
      )
  )
  target_features.append(batch.msa.profile)
  target_features.append(batch.msa.deletion_mean[..., None])

  # Reference structure features
  if append_per_atom_features:
    ref_mask = batch.ref_structure.mask
    element_feat = jax.nn.one_hot(batch.ref_structure.element, 128)
    element_feat = utils.mask_mean(
        mask=ref_mask[..., None], value=element_feat, axis=-2, eps=1e-6
    )
    target_features.append(element_feat)
    pos_feat = batch.ref_structure.positions
    pos_feat = pos_feat.reshape([pos_feat.shape[0], -1])
    target_features.append(pos_feat)
    target_features.append(ref_mask)

  return jnp.concatenate(target_features, axis=-1)


def create_relative_encoding(
    seq_features: features.TokenFeatures,
    max_relative_idx: int,
    max_relative_chain: int,
) -> jax.Array:
  """Add relative position encodings."""
  rel_feats = []
  token_index = seq_features.token_index
  residue_index = seq_features.residue_index
  asym_id = seq_features.asym_id
  entity_id = seq_features.entity_id
  sym_id = seq_features.sym_id

  left_asym_id = asym_id[:, None]
  right_asym_id = asym_id[None, :]

  left_residue_index = residue_index[:, None]
  right_residue_index = residue_index[None, :]

  left_token_index = token_index[:, None]
  right_token_index = token_index[None, :]

  left_entity_id = entity_id[:, None]
  right_entity_id = entity_id[None, :]

  left_sym_id = sym_id[:, None]
  right_sym_id = sym_id[None, :]

  # Embed relative positions using a one-hot embedding of distance along chain
  offset = left_residue_index - right_residue_index
  clipped_offset = jnp.clip(
      offset + max_relative_idx, min=0, max=2 * max_relative_idx
  )
  asym_id_same = left_asym_id == right_asym_id
  final_offset = jnp.where(
      asym_id_same,
      clipped_offset,
      (2 * max_relative_idx + 1) * jnp.ones_like(clipped_offset),
  )
  rel_pos = jax.nn.one_hot(final_offset, 2 * max_relative_idx + 2)
  rel_feats.append(rel_pos)

  # Embed relative token index as a one-hot embedding of distance along residue
  token_offset = left_token_index - right_token_index
  clipped_token_offset = jnp.clip(
      token_offset + max_relative_idx, min=0, max=2 * max_relative_idx
  )
  residue_same = (left_asym_id == right_asym_id) & (
      left_residue_index == right_residue_index
  )
  final_token_offset = jnp.where(
      residue_same,
      clipped_token_offset,
      (2 * max_relative_idx + 1) * jnp.ones_like(clipped_token_offset),
  )
  rel_token = jax.nn.one_hot(final_token_offset, 2 * max_relative_idx + 2)
  rel_feats.append(rel_token)

  # Embed same entity ID
  entity_id_same = left_entity_id == right_entity_id
  rel_feats.append(entity_id_same.astype(rel_pos.dtype)[..., None])

  # Embed relative chain ID inside each symmetry class
  rel_sym_id = left_sym_id - right_sym_id

  max_rel_chain = max_relative_chain

  clipped_rel_chain = jnp.clip(
      rel_sym_id + max_rel_chain, min=0, max=2 * max_rel_chain
  )

  final_rel_chain = jnp.where(
      entity_id_same,
      clipped_rel_chain,
      (2 * max_rel_chain + 1) * jnp.ones_like(clipped_rel_chain),
  )
  rel_chain = jax.nn.one_hot(final_rel_chain, 2 * max_relative_chain + 2)

  rel_feats.append(rel_chain)

  return jnp.concatenate(rel_feats, axis=-1)


def shuffle_msa(
    key: jax.Array, msa: features.MSA
) -> tuple[features.MSA, jax.Array]:
  """Shuffle MSA randomly, return batch with shuffled MSA.

  Args:
    key: rng key for random number generation.
    msa: MSA object to sample msa from.

  Returns:
    Protein with sampled msa.
  """
  key, sample_key = jax.random.split(key)
  # Sample uniformly among sequences with at least one non-masked position.
  logits = (jnp.clip(jnp.sum(msa.mask, axis=-1), 0.0, 1.0) - 1.0) * 1e6
  index_order = gumbel_argsort_sample_idx(sample_key, logits)

  return msa.index_msa_rows(index_order), key
