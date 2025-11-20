# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Evoformer network."""

import functools

from alphafold3.common import base_config
from alphafold3.model import feat_batch
from alphafold3.model import features
from alphafold3.model import model_config
from alphafold3.model.components import haiku_modules as hm
from alphafold3.model.components import utils
from alphafold3.model.network import atom_cross_attention
from alphafold3.model.network import featurization
from alphafold3.model.network import modules
from alphafold3.model.network import template_modules
import haiku as hk
import jax
import jax.numpy as jnp


class Evoformer(hk.Module):
  """Creates 'single' and 'pair' embeddings."""

  class PairformerConfig(modules.PairFormerIteration.Config):  # pytype: disable=invalid-function-definition
    block_remat: bool = False
    remat_block_size: int = 8

  class Config(base_config.BaseConfig):
    """Configuration for Evoformer."""

    max_relative_chain: int = 2
    msa_channel: int = 64
    seq_channel: int = 384
    max_relative_idx: int = 32
    num_msa: int = 1024
    pair_channel: int = 128
    pairformer: 'Evoformer.PairformerConfig' = base_config.autocreate(
        single_transition=base_config.autocreate(),
        single_attention=base_config.autocreate(),
        num_layer=48,
    )
    per_atom_conditioning: atom_cross_attention.AtomCrossAttEncoderConfig = (
        base_config.autocreate(
            per_token_channels=384,
            per_atom_channels=128,
            atom_transformer=base_config.autocreate(
                num_intermediate_factor=2,
                num_blocks=3,
            ),
            per_atom_pair_channels=16,
        )
    )
    template: template_modules.TemplateEmbedding.Config = (
        base_config.autocreate()
    )
    msa_stack: modules.EvoformerIteration.Config = base_config.autocreate()

  def __init__(
      self,
      config: Config,
      global_config: model_config.GlobalConfig,
      name='evoformer',
  ):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def _relative_encoding(
      self, batch: feat_batch.Batch, pair_activations: jnp.ndarray
  ) -> jnp.ndarray:
    """Add relative position encodings."""
    rel_feat = featurization.create_relative_encoding(
        seq_features=batch.token_features,
        max_relative_idx=self.config.max_relative_idx,
        max_relative_chain=self.config.max_relative_chain,
    )
    rel_feat = rel_feat.astype(pair_activations.dtype)

    pair_activations += hm.Linear(
        self.config.pair_channel, name='position_activations'
    )(rel_feat)
    return pair_activations

  @hk.transparent
  def _seq_pair_embedding(
      self,
      token_features: features.TokenFeatures,
      target_feat: jnp.ndarray,
  ) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Generated Pair embedding from sequence."""
    left_single = hm.Linear(self.config.pair_channel, name='left_single')(
        target_feat
    )[:, None]
    right_single = hm.Linear(self.config.pair_channel, name='right_single')(
        target_feat
    )[None]
    dtype = left_single.dtype
    pair_activations = left_single + right_single
    num_residues = pair_activations.shape[0]
    assert pair_activations.shape == (
        num_residues,
        num_residues,
        self.config.pair_channel,
    )
    mask = token_features.mask
    pair_mask = (mask[:, None] * mask[None, :]).astype(dtype)
    assert pair_mask.shape == (num_residues, num_residues)
    return pair_activations, pair_mask  # pytype: disable=bad-return-type  # jax-ndarray
  @hk.transparent
  def _embed_bonds(
      self,
      batch: feat_batch.Batch,
      pair_activations: jnp.ndarray,
  ) -> jnp.ndarray:
    """Embeds bond features and merges into pair activations."""
    # Construct contact matrix.
    num_tokens = batch.token_features.token_index.shape[0]
    contact_matrix = jnp.zeros((num_tokens, num_tokens))

    tokens_to_polymer_ligand_bonds = (
        batch.polymer_ligand_bond_info.tokens_to_polymer_ligand_bonds
    )
    gather_idxs_polymer_ligand = tokens_to_polymer_ligand_bonds.gather_idxs
    # gather_mask entries are per-row masks that indicate valid rows.
    gather_mask_polymer_ligand = (
        tokens_to_polymer_ligand_bonds.gather_mask.prod(axis=1).astype(
            gather_idxs_polymer_ligand.dtype
        )
    )  # shape (N,), dtype same as gather_idxs
    # Keep boolean mask for selecting valid rows later.
    valid_mask_polymer_ligand = gather_mask_polymer_ligand.astype(bool)

    tokens_to_ligand_ligand_bonds = (
        batch.ligand_ligand_bond_info.tokens_to_ligand_ligand_bonds
    )
    gather_idxs_ligand_ligand = tokens_to_ligand_ligand_bonds.gather_idxs
    gather_mask_ligand_ligand = tokens_to_ligand_ligand_bonds.gather_mask.prod(
        axis=1
    ).astype(gather_idxs_ligand_ligand.dtype)
    valid_mask_ligand_ligand = gather_mask_ligand_ligand.astype(bool)

    # Concatenate index arrays and validity masks.
    gather_idxs = jnp.concatenate(
        [gather_idxs_polymer_ligand, gather_idxs_ligand_ligand], axis=0
    )  # shape (M, 2)
    valid_rows = jnp.concatenate(
        [valid_mask_polymer_ligand, valid_mask_ligand_ligand], axis=0
    )  # shape (M,)

    # Select only valid rows. This avoids treating padded (zero) rows as real
    # indices, and handles the empty case safely.
    valid_gather_idxs = gather_idxs[valid_rows]

    # Only perform the scatter if we have at least one valid index.
    def _set_contacts(idx_array):
      return contact_matrix.at[idx_array[:, 0], idx_array[:, 1]].set(1.0)

    # jnp.shape(valid_gather_idxs)[0] is an int-like array; use lax.cond
    # so the function remains jittable.
    contact_matrix = jax.lax.cond(
        valid_gather_idxs.shape[0] > 0,
        _set_contacts,
        lambda _idx: contact_matrix,
        operand=valid_gather_idxs,
    )

    # Ensure padded-index zero entry doesn't accidentally indicate a contact when
    # there were no valid gathers for that entry. (This mirrors previous
    # behavior but now only runs after we populated valid ones.)
    contact_matrix = contact_matrix.at[0, 0].set(0.0)

    bonds_act = hm.Linear(self.config.pair_channel, name='bond_embedding')(
        contact_matrix[:, :, None].astype(pair_activations.dtype)
    )
    return pair_activations + bonds_act


  @hk.transparent
  def _embed_template_pair(
      self,
      batch: feat_batch.Batch,
      pair_activations: jnp.ndarray,
      pair_mask: jnp.ndarray,
      key: jnp.ndarray,
  ) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Embeds Templates and merges into pair activations."""
    dtype = pair_activations.dtype
    key, subkey = jax.random.split(key)
    template_module = template_modules.TemplateEmbedding(
        self.config.template, self.global_config
    )
    templates = batch.templates
    asym_id = batch.token_features.asym_id
    # Construct a mask such that only intra-chain template features are
    # computed, since all templates are for each chain individually.
    multichain_mask = (asym_id[:, None] == asym_id[None, :]).astype(dtype)

    template_fn = functools.partial(template_module, key=subkey)
    template_act = template_fn(
        query_embedding=pair_activations,
        templates=templates,
        multichain_mask_2d=multichain_mask,
        padding_mask_2d=pair_mask,
    )
    return pair_activations + template_act, key

  @hk.transparent
  def _embed_process_msa(
      self,
      msa_batch: features.MSA,
      pair_activations: jnp.ndarray,
      pair_mask: jnp.ndarray,
      key: jnp.ndarray,
      target_feat: jnp.ndarray,
  ) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Processes MSA and returns updated pair activations."""
    dtype = pair_activations.dtype
    msa_batch, key = featurization.shuffle_msa(key, msa_batch)
    msa_batch = featurization.truncate_msa_batch(msa_batch, self.config.num_msa)
    msa_feat = featurization.create_msa_feat(msa_batch).astype(dtype)

    msa_activations = hm.Linear(
        self.config.msa_channel, name='msa_activations'
    )(msa_feat)

    msa_activations += hm.Linear(
        self.config.msa_channel, name='extra_msa_target_feat'
    )(target_feat)[None]
    msa_mask = msa_batch.mask.astype(dtype)

    # Evoformer MSA stack.
    evoformer_input = {'msa': msa_activations, 'pair': pair_activations}
    masks = {'msa': msa_mask, 'pair': pair_mask}

    def evoformer_fn(x):
      return modules.EvoformerIteration(
          self.config.msa_stack, self.global_config, name='msa_stack'
      )(
          activations=x,
          masks=masks,
      )

    evoformer_stack = hk.experimental.layer_stack(
        self.config.msa_stack.num_layer
    )(evoformer_fn)

    evoformer_output = evoformer_stack(evoformer_input)

    return evoformer_output['pair'], key

  def __call__(
      self,
      batch: feat_batch.Batch,
      prev: dict[str, jnp.ndarray],
      target_feat: jnp.ndarray,
      key: jnp.ndarray,
  ) -> dict[str, jnp.ndarray]:

    assert self.global_config.bfloat16 in {'all', 'none'}

    num_residues = target_feat.shape[0]
    assert batch.token_features.aatype.shape == (num_residues,)

    dtype = (
        jnp.bfloat16 if self.global_config.bfloat16 == 'all' else jnp.float32
    )

    with utils.bfloat16_context():
      pair_activations, pair_mask = self._seq_pair_embedding(
          batch.token_features, target_feat
      )

      pair_activations += hm.Linear(
          pair_activations.shape[-1],
          name='prev_embedding',
          initializer=self.global_config.final_init,
      )(
          hm.LayerNorm(name='prev_embedding_layer_norm')(
              prev['pair'].astype(pair_activations.dtype)
          )
      )

      pair_activations = self._relative_encoding(batch, pair_activations)

      pair_activations = self._embed_bonds(
          batch=batch, pair_activations=pair_activations
      )

      pair_activations, key = self._embed_template_pair(
          batch=batch,
          pair_activations=pair_activations,
          pair_mask=pair_mask,
          key=key,
      )
      pair_activations, key = self._embed_process_msa(
          msa_batch=batch.msa,
          pair_activations=pair_activations,
          pair_mask=pair_mask,
          key=key,
          target_feat=target_feat,
      )
      del key  # Unused after this point.

      single_activations = hm.Linear(
          self.config.seq_channel, name='single_activations'
      )(target_feat)

      single_activations += hm.Linear(
          single_activations.shape[-1],
          name='prev_single_embedding',
          initializer=self.global_config.final_init,
      )(
          hm.LayerNorm(name='prev_single_embedding_layer_norm')(
              prev['single'].astype(single_activations.dtype)
          )
      )

      def pairformer_fn(x):
        pairformer_iteration = modules.PairFormerIteration(
            self.config.pairformer,
            self.global_config,
            with_single=True,
            name='trunk_pairformer',
        )
        pair_act, single_act = x
        return pairformer_iteration(
            act=pair_act,
            single_act=single_act,
            pair_mask=pair_mask,
            seq_mask=batch.token_features.mask.astype(dtype),
        )

      pairformer_stack = hk.experimental.layer_stack(
          self.config.pairformer.num_layer
      )(pairformer_fn)

      pair_activations, single_activations = pairformer_stack(
          (pair_activations, single_activations)
      )

      assert pair_activations.shape == (
          num_residues,
          num_residues,
          self.config.pair_channel,
      )
      assert single_activations.shape == (num_residues, self.config.seq_channel)
      assert len(target_feat.shape) == 2
      assert target_feat.shape[0] == num_residues
      output = {
          'single': single_activations,
          'pair': pair_activations,
          'target_feat': target_feat,
      }

    return output
