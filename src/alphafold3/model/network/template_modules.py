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

"""Modules for embedding templates."""

from alphafold3.common import base_config
from alphafold3.constants import residue_names
from alphafold3.jax import geometry
from alphafold3.model import features
from alphafold3.model import model_config
from alphafold3.model import protein_data_processing
from alphafold3.model.components import haiku_modules as hm
from alphafold3.model.network import modules
from alphafold3.model.scoring import scoring
import haiku as hk
import jax
import jax.numpy as jnp


class DistogramFeaturesConfig(base_config.BaseConfig):
  # The left edge of the first bin.
  min_bin: float = 3.25
  # The left edge of the final bin. The final bin catches everything larger than
  # `max_bin`.
  max_bin: float = 50.75
  # The number of bins in the distogram.
  num_bins: int = 39


def dgram_from_positions(positions, config: DistogramFeaturesConfig):
  """Compute distogram from amino acid positions.

  Args:
    positions: (num_res, 3) Position coordinates.
    config: Distogram bin configuration.

  Returns:
    Distogram with the specified number of bins.
  """
  lower_breaks = jnp.linspace(config.min_bin, config.max_bin, config.num_bins)
  lower_breaks = jnp.square(lower_breaks)
  upper_breaks = jnp.concatenate(
      [lower_breaks[1:], jnp.array([1e8], dtype=jnp.float32)], axis=-1
  )
  dist2 = jnp.sum(
      jnp.square(
          jnp.expand_dims(positions, axis=-2)
          - jnp.expand_dims(positions, axis=-3)
      ),
      axis=-1,
      keepdims=True,
  )

  dgram = (dist2 > lower_breaks).astype(jnp.float32) * (
      dist2 < upper_breaks
  ).astype(jnp.float32)
  return dgram


def make_backbone_rigid(
    positions: geometry.Vec3Array,
    mask: jnp.ndarray,
    group_indices: jnp.ndarray,
) -> tuple[geometry.Rigid3Array, jnp.ndarray]:
  """Make backbone Rigid3Array and mask.

  Args:
    positions: (num_res, num_atoms) of atom positions as Vec3Array.
    mask: (num_res, num_atoms) for atom mask.
    group_indices: (num_res, num_group, 3) for atom indices forming groups.

  Returns:
    tuple of backbone Rigid3Array and mask (num_res,).
  """
  backbone_indices = group_indices[:, 0]

  # main backbone frames differ in sidechain frame convention.
  # for sidechain it's (C, CA, N), for backbone it's (N, CA, C)
  # Hence using c, b, a, each of shape (num_res,).
  c, b, a = [backbone_indices[..., i] for i in range(3)]

  slice_index = jax.vmap(lambda x, i: x[i])
  rigid_mask = (
      slice_index(mask, a) * slice_index(mask, b) * slice_index(mask, c)
  ).astype(jnp.float32)

  frame_positions = []
  for indices in [a, b, c]:
    frame_positions.append(
        jax.tree.map(lambda x, idx=indices: slice_index(x, idx), positions)
    )

  rotation = geometry.Rot3Array.from_two_vectors(
      frame_positions[2] - frame_positions[1],
      frame_positions[0] - frame_positions[1],
  )
  rigid = geometry.Rigid3Array(rotation, frame_positions[1])

  return rigid, rigid_mask


class TemplateEmbedding(hk.Module):
  """Embed a set of templates."""

  class Config(base_config.BaseConfig):
    num_channels: int = 64
    template_stack: modules.PairFormerIteration.Config = base_config.autocreate(
        num_layer=2,
        pair_transition=base_config.autocreate(num_intermediate_factor=2),
    )
    dgram_features: DistogramFeaturesConfig = base_config.autocreate()

  def __init__(
      self,
      config: Config,
      global_config: model_config.GlobalConfig,
      name='template_embedding',
  ):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(
      self,
      query_embedding: jnp.ndarray,
      templates: features.Templates,
      padding_mask_2d: jnp.ndarray,
      multichain_mask_2d: jnp.ndarray,
      key: jnp.ndarray,
  ) -> jnp.ndarray:
    """Generate an embedding for a set of templates.

    Args:
      query_embedding: [num_res, num_res, num_channel] a query tensor that will
        be used to attend over the templates to remove the num_templates
        dimension.
      templates: A 'Templates' object.
      padding_mask_2d: [num_res, num_res] Pair mask for attention operations.
      multichain_mask_2d: [num_res, num_res] Pair mask for multichain.
      key: random key generator.

    Returns:
      An embedding of size [num_res, num_res, num_channels]
    """
    c = self.config
    num_residues = query_embedding.shape[0]
    num_templates = templates.aatype.shape[0]
    query_num_channels = query_embedding.shape[2]
    num_atoms = 24
    assert query_embedding.shape == (
        num_residues,
        num_residues,
        query_num_channels,
    )
    assert templates.aatype.shape == (num_templates, num_residues)
    assert templates.atom_positions.shape == (
        num_templates,
        num_residues,
        num_atoms,
        3,
    )
    assert templates.atom_mask.shape == (num_templates, num_residues, num_atoms)
    assert padding_mask_2d.shape == (num_residues, num_residues)

    num_templates = templates.aatype.shape[0]
    num_res, _, query_num_channels = query_embedding.shape

    # Embed each template separately.
    template_embedder = SingleTemplateEmbedding(self.config, self.global_config)

    subkeys = jnp.array(jax.random.split(key, num_templates))

    def scan_fn(carry, x):
      templates, key = x
      embedding = template_embedder(
          query_embedding,
          templates,
          padding_mask_2d,
          multichain_mask_2d,
          key,
      )
      return carry + embedding, None

    scan_init = jnp.zeros(
        (num_res, num_res, c.num_channels), dtype=query_embedding.dtype
    )
    summed_template_embeddings, _ = hk.scan(
        scan_fn, scan_init, (templates, subkeys)
    )

    embedding = summed_template_embeddings / (1e-7 + num_templates)
    embedding = jax.nn.relu(embedding)
    embedding = hm.Linear(
        query_num_channels, initializer='relu', name='output_linear'
    )(embedding)

    assert embedding.shape == (num_residues, num_residues, query_num_channels)
    return embedding


class SingleTemplateEmbedding(hk.Module):
  """Embed a single template."""

  def __init__(
      self,
      config: TemplateEmbedding.Config,
      global_config: model_config.GlobalConfig,
      name='single_template_embedding',
  ):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(
      self,
      query_embedding: jnp.ndarray,
      templates: features.Templates,
      padding_mask_2d: jnp.ndarray,
      multichain_mask_2d: jnp.ndarray,
      key: jnp.ndarray,
  ) -> jnp.ndarray:
    """Build the single template embedding graph.

    Args:
      query_embedding: (num_res, num_res, num_channels) - embedding of the query
        sequence/msa.
      templates: 'Templates' object containing single Template.
      padding_mask_2d: Padding mask (Note: this doesn't care if a template
        exists, unlike the template_pseudo_beta_mask).
      multichain_mask_2d: A mask indicating intra-chain residue pairs, used to
        mask out between chain distances/features when templates are for single
        chains.
      key: Random key generator.

    Returns:
      A template embedding (num_res, num_res, num_channels).
    """
    gc = self.global_config
    c = self.config
    assert padding_mask_2d.dtype == query_embedding.dtype
    dtype = query_embedding.dtype
    num_channels = self.config.num_channels

    def construct_input(
        query_embedding, templates: features.Templates, multichain_mask_2d
    ):

      # Compute distogram feature for the template.
      aatype = templates.aatype
      dense_atom_mask = templates.atom_mask

      dense_atom_positions = templates.atom_positions
      dense_atom_positions *= dense_atom_mask[..., None]

      pseudo_beta_positions, pseudo_beta_mask = scoring.pseudo_beta_fn(
          templates.aatype, dense_atom_positions, dense_atom_mask
      )
      pseudo_beta_mask_2d = (
          pseudo_beta_mask[:, None] * pseudo_beta_mask[None, :]
      )
      pseudo_beta_mask_2d *= multichain_mask_2d
      dgram = dgram_from_positions(
          pseudo_beta_positions, self.config.dgram_features
      )
      dgram *= pseudo_beta_mask_2d[..., None]
      dgram = dgram.astype(dtype)
      pseudo_beta_mask_2d = pseudo_beta_mask_2d.astype(dtype)
      to_concat = [(dgram, 1), (pseudo_beta_mask_2d, 0)]

      aatype = jax.nn.one_hot(
          aatype,
          residue_names.POLYMER_TYPES_NUM_WITH_UNKNOWN_AND_GAP,
          axis=-1,
          dtype=dtype,
      )
      to_concat.append((aatype[None, :, :], 1))
      to_concat.append((aatype[:, None, :], 1))

      # Compute a feature representing the normalized vector between each
      # backbone affine - i.e. in each residues local frame, what direction are
      # each of the other residues.

      template_group_indices = jnp.take(
          protein_data_processing.RESTYPE_RIGIDGROUP_DENSE_ATOM_IDX,
          templates.aatype,
          axis=0,
      )
      rigid, backbone_mask = make_backbone_rigid(
          geometry.Vec3Array.from_array(dense_atom_positions),
          dense_atom_mask,
          template_group_indices.astype(jnp.int32),
      )
      points = rigid.translation
      rigid_vec = rigid[:, None].inverse().apply_to_point(points)
      unit_vector = rigid_vec.normalized()
      unit_vector = [unit_vector.x, unit_vector.y, unit_vector.z]

      unit_vector = [x.astype(dtype) for x in unit_vector]
      backbone_mask = backbone_mask.astype(dtype)

      backbone_mask_2d = backbone_mask[:, None] * backbone_mask[None, :]
      backbone_mask_2d *= multichain_mask_2d
      unit_vector = [x * backbone_mask_2d for x in unit_vector]

      # Note that the backbone_mask takes into account C, CA and N (unlike
      # pseudo beta mask which just needs CB) so we add both masks as features.
      to_concat.extend([(x, 0) for x in unit_vector])
      to_concat.append((backbone_mask_2d, 0))

      query_embedding = hm.LayerNorm(name='query_embedding_norm')(
          query_embedding
      )
      # Allow the template embedder to see the query embedding.  Note this
      # contains the position relative feature, so this is how the network knows
      # which residues are next to each other.
      to_concat.append((query_embedding, 1))

      act = 0

      for i, (x, n_input_dims) in enumerate(to_concat):
        act += hm.Linear(
            num_channels,
            num_input_dims=n_input_dims,
            initializer='relu',
            name=f'template_pair_embedding_{i}',
        )(x)
      return act

    act = construct_input(query_embedding, templates, multichain_mask_2d)

    if c.template_stack.num_layer:

      def template_iteration_fn(x):
        return modules.PairFormerIteration(
            c.template_stack, gc, name='template_embedding_iteration'
        )(act=x, pair_mask=padding_mask_2d)

      template_stack = hk.experimental.layer_stack(c.template_stack.num_layer)(
          template_iteration_fn
      )
      act = template_stack(act)

    act = hm.LayerNorm(name='output_layer_norm')(act)
    return act
