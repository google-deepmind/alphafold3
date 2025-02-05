# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Distogram head."""

from typing import Final

from alphafold3.common import base_config
from alphafold3.model import feat_batch
from alphafold3.model import model_config
from alphafold3.model.components import haiku_modules as hm
import haiku as hk
import jax
import jax.numpy as jnp


_CONTACT_THRESHOLD: Final[float] = 8.0
_CONTACT_EPSILON: Final[float] = 1e-3


class DistogramHead(hk.Module):
  """Distogram head."""

  class Config(base_config.BaseConfig):
    first_break: float = 2.3125
    last_break: float = 21.6875
    num_bins: int = 64

  def __init__(
      self,
      config: Config,
      global_config: model_config.GlobalConfig,
      name='distogram_head',
  ):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(
      self,
      batch: feat_batch.Batch,
      embeddings: dict[str, jnp.ndarray],
      return_distogram: bool = False
  ) -> dict[str, jnp.ndarray]:
    pair_act = embeddings['pair']
    seq_mask = batch.token_features.mask.astype(bool)
    pair_mask = seq_mask[:, None] * seq_mask[None, :]

    left_half_logits = hm.Linear(
        self.config.num_bins,
        initializer=self.global_config.final_init,
        name='half_logits',
    )(pair_act)

    right_half_logits = left_half_logits
    logits = left_half_logits + jnp.swapaxes(right_half_logits, -2, -3)
    probs = jax.nn.softmax(logits, axis=-1)
    breaks = jnp.linspace(
        self.config.first_break,
        self.config.last_break,
        self.config.num_bins - 1,
    )

    bin_tops = jnp.append(breaks, breaks[-1] + (breaks[-1] - breaks[-2]))
    threshold = _CONTACT_THRESHOLD + _CONTACT_EPSILON
    is_contact_bin = 1.0 * (bin_tops <= threshold)
    contact_probs = jnp.einsum(
        'ijk,k->ij', probs, is_contact_bin, precision=jax.lax.Precision.HIGHEST
    )
    contact_probs = pair_mask * contact_probs
    
    if return_distogram:
      return {
          'bin_edges': breaks,
          'contact_probs': contact_probs,
          'distogram': logits,
      }
    return {
        'bin_edges': breaks,
        'contact_probs': contact_probs,
    }
