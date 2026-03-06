# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Diffusion transformer model."""

from alphafold3.common import base_config
from alphafold3.model import model_config
from alphafold3.model.atom_layout import atom_layout
from alphafold3.model.components import haiku_modules as hm
import haiku as hk
import jax
from jax import numpy as jnp
import tokamax


def adaptive_layernorm(x, single_cond, name):
  """Adaptive LayerNorm."""
  # Adopted from Scalable Diffusion Models with Transformers
  # https://arxiv.org/abs/2212.09748
  if single_cond is None:
    x = hm.LayerNorm(name=f'{name}layer_norm', use_fast_variance=False)(x)
  else:
    x = hm.LayerNorm(
        name=f'{name}layer_norm',
        use_fast_variance=False,
        create_scale=False,
        create_offset=False,
    )(x)
    single_cond = hm.LayerNorm(
        name=f'{name}single_cond_layer_norm',
        use_fast_variance=False,
        create_offset=False,
    )(single_cond)
    single_scale = hm.Linear(
        x.shape[-1],
        initializer='zeros',
        use_bias=True,
        name=f'{name}single_cond_scale',
    )(single_cond)
    single_bias = hm.Linear(
        x.shape[-1], initializer='zeros', name=f'{name}single_cond_bias'
    )(single_cond)
    x = jax.nn.sigmoid(single_scale) * x + single_bias
  return x


def adaptive_zero_init(
    x, num_channels, single_cond, global_config: model_config.GlobalConfig, name
):
  """Adaptive zero init, from AdaLN-zero."""
  if single_cond is None:
    output = hm.Linear(
        num_channels,
        initializer=global_config.final_init,
        name=f'{name}transition2',
    )(x)
  else:
    output = hm.Linear(num_channels, name=f'{name}transition2')(x)
    # Init to a small gain, sigmoid(-2) ~ 0.1
    cond = hm.Linear(
        output.shape[-1],
        initializer='zeros',
        use_bias=True,
        bias_init=-2.0,
        name=f'{name}adaptive_zero_cond',
    )(single_cond)
    output = jax.nn.sigmoid(cond) * output
  return output


def transition_block(
    x: jnp.ndarray,
    num_intermediate_factor: int,
    global_config: model_config.GlobalConfig,
    single_cond: jnp.ndarray | None = None,
    use_glu_kernel: bool | None = None,
    name: str = '',
) -> jnp.ndarray:
  """Transition Block."""
  num_channels = x.shape[-1]
  num_intermediates = num_intermediate_factor * num_channels

  x = adaptive_layernorm(x, single_cond, name=f'{name}ffw_')

  if use_glu_kernel is None:
    use_glu_kernel = global_config.use_glu_kernel

  if use_glu_kernel:
    weights, _ = hm.haiku_linear_get_params(
        x,
        num_output=num_intermediates * 2,
        initializer='relu',
        name=f'{name}ffw_transition1',
    )
    weights = jnp.reshape(weights, (len(weights), 2, num_intermediates))
    c = tokamax.gated_linear_unit(x=x, weights=weights, activation=jax.nn.swish)
  else:
    x = hm.Linear(
        num_intermediates * 2, initializer='relu', name=f'{name}ffw_transition1'
    )(x)
    a, b = jnp.split(x, 2, axis=-1)
    c = jax.nn.swish(a) * b

  output = adaptive_zero_init(
      c, num_channels, single_cond, global_config, f'{name}ffw_'
  )
  return output


class SelfAttentionConfig(base_config.BaseConfig):
  num_head: int = 16
  key_dim: int | None = None
  value_dim: int | None = None


def self_attention(
    x: jnp.ndarray,  # (num_tokens, ch)
    mask: jnp.ndarray,  # (num_tokens,)
    pair_logits: jnp.ndarray | None,  # (num_heads, num_tokens, num_tokens)
    config: SelfAttentionConfig,
    global_config: model_config.GlobalConfig,
    single_cond: jnp.ndarray | None = None,  # (num_tokens, ch)
    name: str = '',
) -> jnp.ndarray:
  """Multihead self-attention."""
  assert len(mask.shape) == len(x.shape) - 1, f'{mask.shape}, {x.shape}'
  # bias: ... x heads (1) x query (1) x key
  bias = (1e9 * (mask - 1.0))[..., None, None, :]

  x = adaptive_layernorm(x, single_cond, name=name)

  num_channels = x.shape[-1]
  # Sensible default for when the config keys are missing
  key_dim = config.key_dim if config.key_dim is not None else num_channels
  value_dim = config.value_dim if config.value_dim is not None else num_channels
  num_head = config.num_head
  assert key_dim % num_head == 0, f'{key_dim=} % {num_head=} != 0'
  assert value_dim % num_head == 0, f'{value_dim=} % {num_head=} != 0'
  key_dim = key_dim // num_head
  value_dim = value_dim // num_head

  qk_shape = (num_head, key_dim)
  q = hm.Linear(qk_shape, use_bias=True, name=f'{name}q_projection')(x)
  k = hm.Linear(qk_shape, use_bias=False, name=f'{name}k_projection')(x)

  # In some situations the gradient norms can blow up without running this
  # einsum in float32.
  q = q.astype(jnp.float32)
  k = k.astype(jnp.float32)
  bias = bias.astype(jnp.float32)
  logits = jnp.einsum('...qhc,...khc->...hqk', q * key_dim ** (-0.5), k) + bias
  if pair_logits is not None:
    logits += pair_logits  # (num_heads, seq_len, seq_len)
  weights = jax.nn.softmax(logits, axis=-1)
  weights = jnp.asarray(weights, dtype=x.dtype)

  v_shape = (num_head, value_dim)
  v = hm.Linear(v_shape, use_bias=False, name=f'{name}v_projection')(x)
  weighted_avg = jnp.einsum('...hqk,...khc->...qhc', weights, v)
  weighted_avg = jnp.reshape(weighted_avg, weighted_avg.shape[:-2] + (-1,))

  gate_logits = hm.Linear(
      num_head * value_dim,
      bias_init=1.0,
      initializer='zeros',
      name=f'{name}gating_query',
  )(x)
  weighted_avg *= jax.nn.sigmoid(gate_logits)

  output = adaptive_zero_init(
      weighted_avg, num_channels, single_cond, global_config, name
  )
  return output


class Transformer(hk.Module):
  """Simple transformer stack."""

  class Config(base_config.BaseConfig):
    attention: SelfAttentionConfig = base_config.autocreate()
    num_blocks: int = 24
    block_remat: bool = False
    super_block_size: int = 4
    num_intermediate_factor: int = 2

  def __init__(
      self,
      config: Config,
      global_config: model_config.GlobalConfig,
      name: str = 'transformer',
  ):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(
      self,
      act: jnp.ndarray,
      mask: jnp.ndarray,
      single_cond: jnp.ndarray,
      pair_cond: jnp.ndarray | None,
  ) -> jnp.ndarray:
    def block(act, pair_logits):
      act += self_attention(
          act,
          mask,
          pair_logits,
          self.config.attention,
          self.global_config,
          single_cond,
          name=self.name,
      )
      act += transition_block(
          act,
          self.config.num_intermediate_factor,
          self.global_config,
          single_cond,
          name=self.name,
      )
      return act, None

    # Precompute pair logits for performance
    if pair_cond is None:
      pair_act = None
    else:
      pair_act = hm.LayerNorm(
          name='pair_input_layer_norm',
          use_fast_variance=False,
          create_offset=False,
      )(pair_cond)

    assert self.config.num_blocks % self.config.super_block_size == 0
    num_super_blocks = self.config.num_blocks // self.config.super_block_size

    def super_block(act):
      if pair_act is None:
        pair_logits = None
      else:
        pair_logits = hm.Linear(
            (self.config.super_block_size, self.config.attention.num_head),
            name='pair_logits_projection',
        )(pair_act)
        pair_logits = jnp.transpose(pair_logits, [2, 3, 0, 1])
      return hk.experimental.layer_stack(
          self.config.super_block_size, with_per_layer_inputs=True
      )(block)(act, pair_logits)

    return hk.experimental.layer_stack(
        num_super_blocks, with_per_layer_inputs=True
    )(super_block)(act)[0]


class CrossAttentionConfig(base_config.BaseConfig):
  num_head: int = 4
  key_dim: int = 128
  value_dim: int = 128


def cross_attention(
    x_q: jnp.ndarray,  # (..., Q, C)
    x_k: jnp.ndarray,  # (..., K, C)
    mask_q: jnp.ndarray,  # (..., Q)
    mask_k: jnp.ndarray,  # (..., K)
    config: CrossAttentionConfig,
    global_config: model_config.GlobalConfig,
    pair_logits: jnp.ndarray | None = None,  # (..., Q, K)
    single_cond_q: jnp.ndarray | None = None,  # (..., Q, C)
    single_cond_k: jnp.ndarray | None = None,  # (..., K, C)
    name: str = '',
) -> jnp.ndarray:
  """Multihead self-attention."""
  assert len(mask_q.shape) == len(x_q.shape) - 1, f'{mask_q.shape}, {x_q.shape}'
  assert len(mask_k.shape) == len(x_k.shape) - 1, f'{mask_k.shape}, {x_k.shape}'
  # bias: ... x heads (1) x query x key
  bias = (
      1e9
      * (mask_q - 1.0)[..., None, :, None]
      * (mask_k - 1.0)[..., None, None, :]
  )

  x_q = adaptive_layernorm(x_q, single_cond_q, name=f'{name}q')
  x_k = adaptive_layernorm(x_k, single_cond_k, name=f'{name}k')

  assert config.key_dim % config.num_head == 0
  assert config.value_dim % config.num_head == 0
  key_dim = config.key_dim // config.num_head
  value_dim = config.value_dim // config.num_head

  q = hm.Linear(
      (config.num_head, key_dim), use_bias=True, name=f'{name}q_projection'
  )(x_q)
  k = hm.Linear(
      (config.num_head, key_dim), use_bias=False, name=f'{name}k_projection'
  )(x_k)

  # In some situations the gradient norms can blow up without running this
  # einsum in float32.
  q = q.astype(jnp.float32)
  k = k.astype(jnp.float32)
  bias = bias.astype(jnp.float32)
  logits = jnp.einsum('...qhc,...khc->...hqk', q * key_dim ** (-0.5), k) + bias
  if pair_logits is not None:
    logits += pair_logits
  weights = jax.nn.softmax(logits, axis=-1)
  weights = jnp.asarray(weights, dtype=x_q.dtype)

  v = hm.Linear(
      (config.num_head, value_dim), use_bias=False, name=f'{name}v_projection'
  )(x_k)
  weighted_avg = jnp.einsum('...hqk,...khc->...qhc', weights, v)
  weighted_avg = jnp.reshape(weighted_avg, weighted_avg.shape[:-2] + (-1,))

  gate_logits = hm.Linear(
      config.num_head * value_dim,
      bias_init=1.0,
      initializer='zeros',
      name=f'{name}gating_query',
  )(x_q)
  weighted_avg *= jax.nn.sigmoid(gate_logits)

  output = adaptive_zero_init(
      weighted_avg, x_q.shape[-1], single_cond_q, global_config, name
  )
  return output


class CrossAttTransformer(hk.Module):
  """Transformer that applies cross attention between two sets of subsets."""

  class Config(base_config.BaseConfig):
    num_intermediate_factor: int
    num_blocks: int
    attention: CrossAttentionConfig = base_config.autocreate()

  def __init__(
      self,
      config: Config,
      global_config: model_config.GlobalConfig,
      name: str = 'transformer',
  ):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(
      self,
      queries_act: jnp.ndarray,  # (num_subsets, num_queries, ch)
      queries_mask: jnp.ndarray,  # (num_subsets, num_queries)
      queries_to_keys: atom_layout.GatherInfo,  # (num_subsets, num_keys)
      keys_mask: jnp.ndarray,  # (num_subsets, num_keys)
      queries_single_cond: jnp.ndarray,  # (num_subsets, num_queries, ch)
      keys_single_cond: jnp.ndarray,  # (num_subsets, num_keys, ch)
      pair_cond: jnp.ndarray,  # (num_subsets, num_queries, num_keys, ch)
  ) -> jnp.ndarray:
    def block(queries_act, pair_logits):
      # copy the queries activations to the keys layout
      keys_act = atom_layout.convert(
          queries_to_keys, queries_act, layout_axes=(-3, -2)
      )
      # cross attention
      queries_act += cross_attention(
          x_q=queries_act,
          x_k=keys_act,
          mask_q=queries_mask,
          mask_k=keys_mask,
          config=self.config.attention,
          global_config=self.global_config,
          pair_logits=pair_logits,
          single_cond_q=queries_single_cond,
          single_cond_k=keys_single_cond,
          name=self.name,
      )
      queries_act += transition_block(
          queries_act,
          self.config.num_intermediate_factor,
          self.global_config,
          queries_single_cond,
          name=self.name,
      )
      return queries_act, None

    # Precompute pair logits for performance
    pair_act = hm.LayerNorm(
        name='pair_input_layer_norm',
        use_fast_variance=False,
        create_offset=False,
    )(pair_cond)
    # (num_subsets, num_queries, num_keys, num_blocks, num_heads)
    pair_logits = hm.Linear(
        (self.config.num_blocks, self.config.attention.num_head),
        name='pair_logits_projection',
    )(pair_act)
    # (num_block, num_subsets, num_heads, num_queries, num_keys)
    pair_logits = jnp.transpose(pair_logits, [3, 0, 4, 1, 2])

    return hk.experimental.layer_stack(
        self.config.num_blocks, with_per_layer_inputs=True
    )(block)(queries_act, pair_logits)[0]
