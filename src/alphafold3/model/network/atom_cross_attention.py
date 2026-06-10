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

"""Per-atom cross attention."""
import dataclasses

from alphafold3.common import base_config
from alphafold3.model import feat_batch
from alphafold3.model import model_config
from alphafold3.model.atom_layout import atom_layout
from alphafold3.model.components import haiku_modules as hm
from alphafold3.model.components import utils
from alphafold3.model.network import diffusion_transformer
import jax
import jax.numpy as jnp


class AtomCrossAttEncoderConfig(base_config.BaseConfig):
  per_token_channels: int = 768
  per_atom_channels: int = 128
  atom_transformer: diffusion_transformer.CrossAttTransformer.Config = (
      base_config.autocreate(num_intermediate_factor=2, num_blocks=3)
  )
  per_atom_pair_channels: int = 16


def _per_atom_conditioning(
    config: AtomCrossAttEncoderConfig, batch: feat_batch.Batch, name: str
) -> tuple[jnp.ndarray, jnp.ndarray]:
  """computes single and pair conditioning for all atoms in each token."""

  c = config
  # Compute per-atom single conditioning
  # Shape (num_tokens, num_dense, channels)
  act = hm.Linear(
      c.per_atom_channels, precision='highest', name=f'{name}_embed_ref_pos'
  )(batch.ref_structure.positions)
  act += hm.Linear(c.per_atom_channels, name=f'{name}_embed_ref_mask')(
      batch.ref_structure.mask.astype(jnp.float32)[:, :, None]
  )
  # Element is encoded as atomic number if the periodic table, so
  # 128 should be fine.
  act += hm.Linear(c.per_atom_channels, name=f'{name}_embed_ref_element')(
      jax.nn.one_hot(batch.ref_structure.element, 128)
  )
  act += hm.Linear(c.per_atom_channels, name=f'{name}_embed_ref_charge')(
      jnp.arcsinh(batch.ref_structure.charge)[:, :, None]
  )
  # Characters are encoded as ASCII code minus 32, so we need 64 classes,
  # to encode all standard ASCII characters between 32 and 96.
  atom_name_chars_1hot = jax.nn.one_hot(batch.ref_structure.atom_name_chars, 64)
  num_token, num_dense, _ = act.shape
  act += hm.Linear(c.per_atom_channels, name=f'{name}_embed_ref_atom_name')(
      atom_name_chars_1hot.reshape(num_token, num_dense, -1)
  )
  act *= batch.ref_structure.mask.astype(jnp.float32)[:, :, None]

  # Compute pair conditioning
  # shape (num_tokens, num_dense, num_dense, channels)
  # Embed single features
  row_act = hm.Linear(
      c.per_atom_pair_channels, name=f'{name}_single_to_pair_cond_row'
  )(jax.nn.relu(act))
  col_act = hm.Linear(
      c.per_atom_pair_channels, name=f'{name}_single_to_pair_cond_col'
  )(jax.nn.relu(act))
  pair_act = row_act[:, :, None, :] + col_act[:, None, :, :]
  # Embed pairwise offsets
  pair_act += hm.Linear(
      c.per_atom_pair_channels,
      precision='highest',
      name=f'{name}_embed_pair_offsets',
  )(
      batch.ref_structure.positions[:, :, None, :]
      - batch.ref_structure.positions[:, None, :, :]
  )
  # Embed pairwise inverse squared distances
  sq_dists = jnp.sum(
      jnp.square(
          batch.ref_structure.positions[:, :, None, :]
          - batch.ref_structure.positions[:, None, :, :]
      ),
      axis=-1,
  )
  pair_act += hm.Linear(
      c.per_atom_pair_channels, name=f'{name}_embed_pair_distances'
  )(1.0 / (1 + sq_dists[:, :, :, None]))

  return act, pair_act


@dataclasses.dataclass(frozen=True)
class AtomCrossAttEncoderOutput:
  token_act: jnp.ndarray  # (num_tokens, ch)
  skip_connection: jnp.ndarray  # (num_subsets, num_queries, ch)
  queries_mask: jnp.ndarray  # (num_subsets, num_queries)
  queries_single_cond: jnp.ndarray  # (num_subsets, num_queries, ch)
  keys_mask: jnp.ndarray  # (num_subsets, num_keys)
  keys_single_cond: jnp.ndarray  # (num_subsets, num_keys, ch)
  pair_cond: jnp.ndarray  # (num_subsets, num_queries, num_keys, ch)


jax.tree_util.register_dataclass(
    AtomCrossAttEncoderOutput,
    data_fields=[f.name for f in dataclasses.fields(AtomCrossAttEncoderOutput)],
    meta_fields=[],
)


def atom_cross_att_encoder(
    token_atoms_act: jnp.ndarray | None,  # (num_tokens, max_atoms_per_token, 3)
    trunk_single_cond: jnp.ndarray | None,  # (num_tokens, ch)
    trunk_pair_cond: jnp.ndarray | None,  # (num_tokens, num_tokens, ch)
    config: AtomCrossAttEncoderConfig,
    global_config: model_config.GlobalConfig,
    batch: feat_batch.Batch,
    name: str,
) -> AtomCrossAttEncoderOutput:
  """Cross-attention on flat atom subsets and mapping to per-token features."""
  c = config

  # Compute single conditioning from atom meta data and convert to queries
  # layout.
  # (num_subsets, num_queries, channels)
  token_atoms_single_cond, _ = _per_atom_conditioning(config, batch, name)
  token_atoms_mask = batch.predicted_structure_info.atom_mask
  queries_single_cond = atom_layout.convert(
      batch.atom_cross_att.token_atoms_to_queries,
      token_atoms_single_cond,
      layout_axes=(-3, -2),
  )
  queries_mask = atom_layout.convert(
      batch.atom_cross_att.token_atoms_to_queries,
      token_atoms_mask,
      layout_axes=(-2, -1),
  )

  # If provided, broadcast single conditioning from trunk to all queries
  if trunk_single_cond is not None:
    trunk_single_cond = hm.Linear(
        c.per_atom_channels,
        precision='highest',
        initializer=global_config.final_init,
        name=f'{name}_embed_trunk_single_cond',
    )(
        hm.LayerNorm(
            use_fast_variance=False,
            create_offset=False,
            name=f'{name}_lnorm_trunk_single_cond',
        )(trunk_single_cond)
    )
    queries_single_cond += atom_layout.convert(
        batch.atom_cross_att.tokens_to_queries,
        trunk_single_cond,
        layout_axes=(-2,),
    )

  if token_atoms_act is None:
    # if no token_atoms_act is given (e.g. begin of evoformer), we use the
    # static conditioning only
    queries_act = queries_single_cond
  else:
    # Convert token_atoms_act to queries layout and map to per_atom_channels
    # (num_subsets, num_queries, channels)
    queries_act = atom_layout.convert(
        batch.atom_cross_att.token_atoms_to_queries,
        token_atoms_act,
        layout_axes=(-3, -2),
    )
    queries_act = hm.Linear(
        c.per_atom_channels,
        precision='highest',
        name=f'{name}_atom_positions_to_features',
    )(queries_act)
    queries_act *= queries_mask[..., None]
    queries_act += queries_single_cond

  # Gather the keys from the queries.
  keys_single_cond = atom_layout.convert(
      batch.atom_cross_att.queries_to_keys,
      queries_single_cond,
      layout_axes=(-3, -2),
  )
  keys_mask = atom_layout.convert(
      batch.atom_cross_att.queries_to_keys, queries_mask, layout_axes=(-2, -1)
  )

  # Embed single features into the pair conditioning.
  # shape (num_subsets, num_queries, num_keys, ch)
  row_act = hm.Linear(
      c.per_atom_pair_channels, name=f'{name}_single_to_pair_cond_row'
  )(jax.nn.relu(queries_single_cond))
  pair_cond_keys_input = atom_layout.convert(
      batch.atom_cross_att.queries_to_keys,
      queries_single_cond,
      layout_axes=(-3, -2),
  )
  col_act = hm.Linear(
      c.per_atom_pair_channels, name=f'{name}_single_to_pair_cond_col'
  )(jax.nn.relu(pair_cond_keys_input))
  pair_act = row_act[:, :, None, :] + col_act[:, None, :, :]

  if trunk_pair_cond is not None:
    # If provided, broadcast the pair conditioning for the trunk (evoformer
    # pairs) to the atom pair activations. This should boost ligands, but also
    # help for cross attention within proteins, because we always have atoms
    # from multiple residues in a subset.
    # Map trunk pair conditioning to per_atom_pair_channels
    # (num_tokens, num_tokens, per_atom_pair_channels)
    trunk_pair_cond = hm.Linear(
        c.per_atom_pair_channels,
        precision='highest',
        initializer=global_config.final_init,
        name=f'{name}_embed_trunk_pair_cond',
    )(
        hm.LayerNorm(
            use_fast_variance=False,
            create_offset=False,
            name=f'{name}_lnorm_trunk_pair_cond',
        )(trunk_pair_cond)
    )

    # Create the GatherInfo into a flattened trunk_pair_cond from the
    # queries and keys gather infos.
    num_tokens = trunk_pair_cond.shape[0]
    # (num_subsets, num_queries)
    tokens_to_queries = batch.atom_cross_att.tokens_to_queries
    # (num_subsets, num_keys)
    tokens_to_keys = batch.atom_cross_att.tokens_to_keys
    # (num_subsets, num_queries, num_keys)
    trunk_pair_to_atom_pair = atom_layout.GatherInfo(
        gather_idxs=(
            num_tokens * tokens_to_queries.gather_idxs[:, :, None]
            + tokens_to_keys.gather_idxs[:, None, :]
        ),
        gather_mask=(
            tokens_to_queries.gather_mask[:, :, None]
            & tokens_to_keys.gather_mask[:, None, :]
        ),
        input_shape=jnp.array((num_tokens, num_tokens)),
    )
    # Gather the conditioning and add it to the atom-pair activations.
    pair_act += atom_layout.convert(
        trunk_pair_to_atom_pair, trunk_pair_cond, layout_axes=(-3, -2)
    )

  # Embed pairwise offsets
  queries_ref_pos = atom_layout.convert(
      batch.atom_cross_att.token_atoms_to_queries,
      batch.ref_structure.positions,
      layout_axes=(-3, -2),
  )
  queries_ref_space_uid = atom_layout.convert(
      batch.atom_cross_att.token_atoms_to_queries,
      batch.ref_structure.ref_space_uid,
      layout_axes=(-2, -1),
  )
  keys_ref_pos = atom_layout.convert(
      batch.atom_cross_att.queries_to_keys,
      queries_ref_pos,
      layout_axes=(-3, -2),
  )
  keys_ref_space_uid = atom_layout.convert(
      batch.atom_cross_att.queries_to_keys,
      batch.ref_structure.ref_space_uid,
      layout_axes=(-2, -1),
  )

  offsets_valid = (
      queries_ref_space_uid[:, :, None] == keys_ref_space_uid[:, None, :]
  )
  offsets = queries_ref_pos[:, :, None, :] - keys_ref_pos[:, None, :, :]
  pair_act += (
      hm.Linear(
          c.per_atom_pair_channels,
          precision='highest',
          name=f'{name}_embed_pair_offsets',
      )(offsets)
      * offsets_valid[:, :, :, None]
  )

  # Embed pairwise inverse squared distances
  sq_dists = jnp.sum(jnp.square(offsets), axis=-1)
  pair_act += (
      hm.Linear(c.per_atom_pair_channels, name=f'{name}_embed_pair_distances')(
          1.0 / (1 + sq_dists[:, :, :, None])
      )
      * offsets_valid[:, :, :, None]
  )
  # Embed offsets valid mask
  pair_act += hm.Linear(
      c.per_atom_pair_channels, name=f'{name}_embed_pair_offsets_valid'
  )(offsets_valid[:, :, :, None].astype(jnp.float32))

  # Run a small MLP on the pair acitvations
  pair_act2 = hm.Linear(
      c.per_atom_pair_channels, initializer='relu', name=f'{name}_pair_mlp_1'
  )(jax.nn.relu(pair_act))
  pair_act2 = hm.Linear(
      c.per_atom_pair_channels, initializer='relu', name=f'{name}_pair_mlp_2'
  )(jax.nn.relu(pair_act2))
  pair_act += hm.Linear(
      c.per_atom_pair_channels,
      initializer=global_config.final_init,
      name=f'{name}_pair_mlp_3',
  )(jax.nn.relu(pair_act2))

  # Run the atom cross attention transformer.
  queries_act = diffusion_transformer.CrossAttTransformer(
      c.atom_transformer, global_config, name=f'{name}_atom_transformer_encoder'
  )(
      queries_act=queries_act,
      queries_mask=queries_mask,
      queries_to_keys=batch.atom_cross_att.queries_to_keys,
      keys_mask=keys_mask,
      queries_single_cond=queries_single_cond,
      keys_single_cond=keys_single_cond,
      pair_cond=pair_act,
  )
  queries_act *= queries_mask[..., None]
  skip_connection = queries_act

  # Convert back to token-atom layout and aggregate to tokens
  queries_act = hm.Linear(
      c.per_token_channels, name=f'{name}_project_atom_features_for_aggr'
  )(queries_act)
  token_atoms_act = atom_layout.convert(
      batch.atom_cross_att.queries_to_token_atoms,
      queries_act,
      layout_axes=(-3, -2),
  )
  token_act = utils.mask_mean(
      token_atoms_mask[..., None], jax.nn.relu(token_atoms_act), axis=-2
  )

  return AtomCrossAttEncoderOutput(
      token_act=token_act,
      skip_connection=skip_connection,
      queries_mask=queries_mask,
      queries_single_cond=queries_single_cond,
      keys_mask=keys_mask,
      keys_single_cond=keys_single_cond,
      pair_cond=pair_act,
  )


class AtomCrossAttDecoderConfig(base_config.BaseConfig):
  per_atom_channels: int = 128
  atom_transformer: diffusion_transformer.CrossAttTransformer.Config = (
      base_config.autocreate(num_intermediate_factor=2, num_blocks=3)
  )


def atom_cross_att_decoder(
    token_act: jnp.ndarray,  # (num_tokens, ch)
    enc: AtomCrossAttEncoderOutput,
    config: AtomCrossAttDecoderConfig,
    global_config: model_config.GlobalConfig,
    batch: feat_batch.Batch,
    name: str,
):  # (num_tokens, max_atoms_per_token, 3)
  """Mapping to per-atom features and self-attention on subsets."""
  c = config
  # map per-token act down to per_atom channels
  token_act = hm.Linear(
      c.per_atom_channels, name=f'{name}_project_token_features_for_broadcast'
  )(token_act)
  # Broadcast to token-atoms layout and convert to queries layout.
  num_token, max_atoms_per_token = (
      batch.atom_cross_att.queries_to_token_atoms.shape
  )
  token_atom_act = jnp.broadcast_to(
      token_act[:, None, :],
      (num_token, max_atoms_per_token, c.per_atom_channels),
  )
  queries_act = atom_layout.convert(
      batch.atom_cross_att.token_atoms_to_queries,
      token_atom_act,
      layout_axes=(-3, -2),
  )
  queries_act += enc.skip_connection
  queries_act *= enc.queries_mask[..., None]

  # Run the atom cross attention transformer.
  queries_act = diffusion_transformer.CrossAttTransformer(
      c.atom_transformer, global_config, name=f'{name}_atom_transformer_decoder'
  )(
      queries_act=queries_act,
      queries_mask=enc.queries_mask,
      queries_to_keys=batch.atom_cross_att.queries_to_keys,
      keys_mask=enc.keys_mask,
      queries_single_cond=enc.queries_single_cond,
      keys_single_cond=enc.keys_single_cond,
      pair_cond=enc.pair_cond,
  )
  queries_act *= enc.queries_mask[..., None]
  queries_act = hm.LayerNorm(
      use_fast_variance=False,
      create_offset=False,
      name=f'{name}_atom_features_layer_norm',
  )(queries_act)
  queries_position_update = hm.Linear(
      3,
      initializer=global_config.final_init,
      precision='highest',
      name=f'{name}_atom_features_to_position_update',
  )(queries_act)
  position_update = atom_layout.convert(
      batch.atom_cross_att.queries_to_token_atoms,
      queries_position_update,
      layout_axes=(-3, -2),
  )
  return position_update
