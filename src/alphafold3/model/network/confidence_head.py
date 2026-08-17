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

"""Confidence Head."""

from alphafold3.common import base_config
from alphafold3.model import model_config
from alphafold3.model.atom_layout import atom_layout
from alphafold3.model.components import haiku_modules as hm
from alphafold3.model.components import utils
from alphafold3.model.network import modules
from alphafold3.model.network import template_modules
import haiku as hk
import jax
import jax.numpy as jnp


def _safe_norm(x, keepdims, axis, eps=1e-8):
  return jnp.sqrt(eps + jnp.sum(jnp.square(x), axis=axis, keepdims=keepdims))


class ConfidenceHead(hk.Module):
  """Head to predict the distance errors in a prediction."""

  class PAEConfig(base_config.BaseConfig):
    max_error_bin: float = 31.0
    num_bins: int = 64

  class Config(base_config.BaseConfig):
    """Configuration for ConfidenceHead."""

    pairformer: modules.PairFormerIteration.Config = base_config.autocreate(
        single_attention=base_config.autocreate(),
        single_transition=base_config.autocreate(),
        num_layer=4,
    )
    max_error_bin: float = 31.0
    num_plddt_bins: int = 50
    num_bins: int = 64
    no_embedding_prob: float = 0.2
    pae: 'ConfidenceHead.PAEConfig' = base_config.autocreate()
    dgram_features: template_modules.DistogramFeaturesConfig = (
        base_config.autocreate()
    )

  def __init__(
      self,
      config: Config,
      global_config: model_config.GlobalConfig,
      name='confidence_head',
  ):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def _embed_features(
      self,
      dense_atom_positions,
      token_atoms_to_pseudo_beta,
      pair_mask,
      pair_act,
      target_feat,
  ):
    out = hm.Linear(pair_act.shape[-1], name='left_target_feat_project')(
        target_feat
    ).astype(pair_act.dtype)
    out += hm.Linear(pair_act.shape[-1], name='right_target_feat_project')(
        target_feat
    ).astype(pair_act.dtype)[:, None]
    positions = atom_layout.convert(
        token_atoms_to_pseudo_beta,
        dense_atom_positions,
        layout_axes=(-3, -2),
    )
    dgram = template_modules.dgram_from_positions(
        positions, self.config.dgram_features
    )
    dgram *= pair_mask[..., None]

    out += hm.Linear(pair_act.shape[-1], name='distogram_feat_project')(
        dgram.astype(pair_act.dtype)
    )
    return out

  def __call__(
      self,
      dense_atom_positions: jnp.ndarray,
      embeddings: dict[str, jnp.ndarray],
      seq_mask: jnp.ndarray,
      token_atoms_to_pseudo_beta: atom_layout.GatherInfo,
      asym_id: jnp.ndarray,
  ) -> dict[str, jnp.ndarray]:
    """Builds ConfidenceHead module.

    Arguments:
      dense_atom_positions: [N_res, N_atom, 3] array of positions.
      embeddings: Dictionary of representations.
      seq_mask: Sequence mask.
      token_atoms_to_pseudo_beta: Pseudo beta info for atom tokens.
      asym_id: Asym ID token features.

    Returns:
      Dictionary of results.
    """
    dtype = (
        jnp.bfloat16 if self.global_config.bfloat16 == 'all' else jnp.float32
    )
    with utils.bfloat16_context():
      seq_mask_cast = seq_mask.astype(dtype)
      pair_mask = seq_mask_cast[:, None] * seq_mask_cast[None, :]
      pair_mask = pair_mask.astype(dtype)

      pair_act = embeddings['pair'].astype(dtype)
      single_act = embeddings['single'].astype(dtype)
      target_feat = embeddings['target_feat'].astype(dtype)

      num_residues = seq_mask.shape[0]
      num_pair_channels = pair_act.shape[2]

      pair_act += self._embed_features(
          dense_atom_positions,
          token_atoms_to_pseudo_beta,
          pair_mask,
          pair_act,
          target_feat,
      )

      def pairformer_fn(act):
        pair_act, single_act = act
        return modules.PairFormerIteration(
            self.config.pairformer,
            self.global_config,
            with_single=True,
            name='confidence_pairformer',
        )(
            act=pair_act,
            single_act=single_act,
            pair_mask=pair_mask,
            seq_mask=seq_mask,
        )

      pairformer_stack = hk.experimental.layer_stack(
          self.config.pairformer.num_layer
      )(pairformer_fn)

      pair_act, single_act = pairformer_stack((pair_act, single_act))
      pair_act = pair_act.astype(jnp.float32)
      assert pair_act.shape == (num_residues, num_residues, num_pair_channels)

      # Produce logits to predict a distogram of pairwise distance errors
      # between the input prediction and the ground truth.

      # Shape (num_res, num_res, num_bins)
      left_distance_logits = hm.Linear(
          self.config.num_bins,
          initializer=self.global_config.final_init,
          name='left_half_distance_logits',
      )(hm.LayerNorm(name='logits_ln')(pair_act))
      right_distance_logits = left_distance_logits
      distance_logits = left_distance_logits + jnp.swapaxes(  # Symmetrize.
          right_distance_logits, -2, -3
      )
      # Shape (num_bins,)
      distance_breaks = jnp.linspace(
          0.0, self.config.max_error_bin, self.config.num_bins - 1
      )

      step = distance_breaks[1] - distance_breaks[0]

      # Add half-step to get the center
      bin_centers = distance_breaks + step / 2
      # Add a catch-all bin at the end.
      bin_centers = jnp.concatenate(
          [bin_centers, bin_centers[-1:] + step], axis=0
      )

      distance_probs = jax.nn.softmax(distance_logits, axis=-1)

      pred_distance_error = (
          jnp.sum(distance_probs * bin_centers, axis=-1) * pair_mask
      )
      average_pred_distance_error = jnp.sum(
          pred_distance_error, axis=[-2, -1]
      ) / jnp.sum(pair_mask, axis=[-2, -1])

      # Predicted aligned error
      pae_outputs = {}
      # Shape (num_res, num_res, num_bins)
      pae_logits = hm.Linear(
          self.config.pae.num_bins,
          initializer=self.global_config.final_init,
          name='pae_logits',
      )(hm.LayerNorm(name='pae_logits_ln')(pair_act))
      # Shape (num_bins,)
      pae_breaks = jnp.linspace(
          0.0, self.config.pae.max_error_bin, self.config.pae.num_bins - 1
      )
      step = pae_breaks[1] - pae_breaks[0]
      # Add half-step to get the center
      bin_centers = pae_breaks + step / 2
      # Add a catch-all bin at the end.
      bin_centers = jnp.concatenate(
          [bin_centers, bin_centers[-1:] + step], axis=0
      )
      pae_probs = jax.nn.softmax(pae_logits, axis=-1)

      seq_mask_bool = seq_mask.astype(bool)
      pair_mask_bool = seq_mask_bool[:, None] * seq_mask_bool[None, :]
      pae = jnp.sum(pae_probs * bin_centers, axis=-1) * pair_mask_bool
      pae_outputs.update({
          'full_pae': pae,
      })

    # The pTM is computed outside of bfloat16 context.
    tmscore_adjusted_pae_global, tmscore_adjusted_pae_interface = (
        self._get_tmscore_adjusted_pae(
            asym_id=asym_id,
            seq_mask=seq_mask,
            pair_mask=pair_mask_bool,
            bin_centers=bin_centers,
            pae_probs=pae_probs,
        )
    )
    pae_outputs.update({
        'tmscore_adjusted_pae_global': tmscore_adjusted_pae_global,
        'tmscore_adjusted_pae_interface': tmscore_adjusted_pae_interface,
    })
    single_act = single_act.astype('float32')

    # pLDDT
    # Shape (num_res, num_atom, num_bins)
    plddt_logits = hm.Linear(
        (dense_atom_positions.shape[-2], self.config.num_plddt_bins),
        initializer=self.global_config.final_init,
        name='plddt_logits',
    )(hm.LayerNorm(name='plddt_logits_ln')(single_act))

    bin_width = 1.0 / self.config.num_plddt_bins
    bin_centers = jnp.arange(0.5 * bin_width, 1.0, bin_width)
    predicted_lddt = jnp.sum(
        jax.nn.softmax(plddt_logits, axis=-1) * bin_centers, axis=-1
    )
    predicted_lddt = predicted_lddt * 100.0

    # Experimentally resolved
    # Shape (num_res, num_atom, 2)
    experimentally_resolved_logits = hm.Linear(
        (dense_atom_positions.shape[-2], 2),
        initializer=self.global_config.final_init,
        name='experimentally_resolved_logits',
    )(hm.LayerNorm(name='experimentally_resolved_ln')(single_act))

    predicted_experimentally_resolved = jax.nn.softmax(
        experimentally_resolved_logits, axis=-1
    )[..., 1]

    return {
        'predicted_lddt': predicted_lddt,
        'predicted_experimentally_resolved': predicted_experimentally_resolved,
        'full_pde': pred_distance_error,
        'average_pde': average_pred_distance_error,
        **pae_outputs,
    }

  def _get_tmscore_adjusted_pae(
      self,
      asym_id: jnp.ndarray,
      seq_mask: jnp.ndarray,
      pair_mask: jnp.ndarray,
      bin_centers: jnp.ndarray,
      pae_probs: jnp.ndarray,
  ):
    def get_tmscore_adjusted_pae(num_interface_tokens, bin_centers, pae_probs):
      # Clip to avoid negative/undefined d0.
      clipped_num_res = jnp.maximum(num_interface_tokens, 19)

      # Compute d_0(num_res) as defined by TM-score, eqn. (5) in
      # http://zhanglab.ccmb.med.umich.edu/papers/2004_3.pdf
      # Yang & Skolnick "Scoring function for automated
      # assessment of protein structure template quality" 2004.
      d0 = 1.24 * (clipped_num_res - 15) ** (1.0 / 3) - 1.8

      # Make compatible with [num_tokens, num_tokens, num_bins]
      d0 = d0[:, :, None]
      bin_centers = bin_centers[None, None, :]

      # TM-Score term for every bin.
      tm_per_bin = 1.0 / (1 + jnp.square(bin_centers) / jnp.square(d0))
      # E_distances tm(distance).
      predicted_tm_term = jnp.sum(pae_probs * tm_per_bin, axis=-1)
      return predicted_tm_term

    # Interface version
    x = asym_id[None, :] == asym_id[:, None]
    num_chain_tokens = jnp.sum(x * pair_mask, axis=-1)
    num_interface_tokens = num_chain_tokens[None, :] + num_chain_tokens[:, None]
    # Don't double-count within a single chain
    num_interface_tokens -= x * (num_interface_tokens // 2)
    num_interface_tokens = num_interface_tokens * pair_mask

    num_global_tokens = jnp.full(
        shape=pair_mask.shape, fill_value=seq_mask.sum()
    )

    assert num_global_tokens.dtype == 'int32'
    assert num_interface_tokens.dtype == 'int32'
    global_apae = get_tmscore_adjusted_pae(
        num_global_tokens, bin_centers, pae_probs
    )
    interface_apae = get_tmscore_adjusted_pae(
        num_interface_tokens, bin_centers, pae_probs
    )
    return global_apae, interface_apae
