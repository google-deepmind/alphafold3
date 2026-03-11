# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Haiku modules for the Diffuser model."""

from collections.abc import Sequence
from typing import Literal

from alphafold3.common import base_config
from alphafold3.model import model_config
from alphafold3.model.components import haiku_modules as hm
from alphafold3.model.components import mapping
from alphafold3.model.network import diffusion_transformer
import haiku as hk
import jax
import jax.numpy as jnp
import tokamax


def get_shard_size(
    num_residues: int, shard_spec: Sequence[tuple[int | None, int | None]]
) -> int | None:
  shard_size = shard_spec[0][-1]
  for num_residues_upper_bound, num_residues_shard_size in shard_spec:
    shard_size = num_residues_shard_size
    if (
        num_residues_upper_bound is None
        or num_residues <= num_residues_upper_bound
    ):
      break
  return shard_size


class TransitionBlock(hk.Module):
  """Transition block for transformer."""

  class Config(base_config.BaseConfig):
    num_intermediate_factor: int = 4
    use_glu_kernel: bool = True

  def __init__(
      self, config: Config, global_config: model_config.GlobalConfig, *, name
  ):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self, act, broadcast_dim=0):
    num_channels = act.shape[-1]

    num_intermediate = int(num_channels * self.config.num_intermediate_factor)

    act = hm.LayerNorm(name='input_layer_norm')(act)

    if self.global_config.use_glu_kernel:
      weights, _ = hm.haiku_linear_get_params(
          act,
          num_output=num_intermediate * 2,
          initializer='relu',
          name='transition1',
      )
      weights = jnp.reshape(weights, (len(weights), 2, num_intermediate))
      c = tokamax.gated_linear_unit(
          x=act, weights=weights, activation=jax.nn.swish
      )
    else:
      act = hm.Linear(
          num_intermediate * 2, initializer='relu', name='transition1'
      )(act)
      a, b = jnp.split(act, 2, axis=-1)
      c = jax.nn.swish(a) * b

    return hm.Linear(
        num_channels,
        initializer=self.global_config.final_init,
        name='transition2',
    )(c)


class MSAAttention(hk.Module):
  """MSA Attention."""

  class Config(base_config.BaseConfig):
    num_head: int = 8

  def __init__(
      self, config: Config, global_config: model_config.GlobalConfig, *, name
  ):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self, act, mask, pair_act):
    act = hm.LayerNorm(name='act_norm')(act)
    pair_act = hm.LayerNorm(name='pair_norm')(pair_act)
    logits = hm.Linear(
        self.config.num_head, use_bias=False, name='pair_logits'
    )(pair_act)
    logits = jnp.transpose(logits, [2, 0, 1])
    logits += 1e9 * (jnp.max(mask, axis=0) - 1.0)
    weights = jax.nn.softmax(logits, axis=-1)
    num_channels = act.shape[-1]
    value_dim = num_channels // self.config.num_head
    v = hm.Linear(
        [self.config.num_head, value_dim], use_bias=False, name='v_projection'
    )(act)
    v_avg = jnp.einsum('hqk, bkhc -> bqhc', weights, v)
    v_avg = jnp.reshape(v_avg, v_avg.shape[:-2] + (-1,))
    gate_values = hm.Linear(
        self.config.num_head * value_dim,
        bias_init=1.0,
        initializer='zeros',
        name='gating_query',
    )(act)
    v_avg *= jax.nn.sigmoid(gate_values)

    return hm.Linear(
        num_channels,
        initializer=self.global_config.final_init,
        name='output_projection',
    )(v_avg)


class GridSelfAttention(hk.Module):
  """Self attention that is either per-sequence or per-residue."""

  class Config(base_config.BaseConfig):
    num_head: int = 4

  def __init__(
      self,
      config: Config,
      global_config: model_config.GlobalConfig,
      transpose: bool,
      *,
      name: str,
  ):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config
    self.transpose = transpose

  @hk.transparent
  def _attention(
      self,
      act,
      mask,
      bias,
  ):
    num_channels = act.shape[-1]
    assert num_channels % self.config.num_head == 0
    # Triton requires a minimum dimension of 16 for doing matmul.
    qkv_dim = max(num_channels // self.config.num_head, 16)

    qkv_shape = (self.config.num_head, qkv_dim)
    q = hm.Linear(
        qkv_shape, use_bias=False, name='q_projection', transpose_weights=True
    )(act)
    k = hm.Linear(
        qkv_shape, use_bias=False, name='k_projection', transpose_weights=True
    )(act)
    v = hm.Linear(qkv_shape, use_bias=False, name='v_projection')(act)

    # Dot product attention requires the bias term to have a batch dimension.
    bias = jnp.expand_dims(bias, 0)

    weighted_avg = tokamax.dot_product_attention(
        q,
        k,
        v,
        mask=mask,
        bias=bias,
        implementation=self.global_config.flash_attention_implementation,
    )
    weighted_avg = jnp.reshape(weighted_avg, weighted_avg.shape[:-2] + (-1,))

    gate_values = hm.Linear(
        self.config.num_head * qkv_dim,
        bias_init=1.0,
        initializer='zeros',
        transpose_weights=True,
        name='gating_query',
    )(act)
    weighted_avg *= jax.nn.sigmoid(gate_values)

    return hm.Linear(
        num_channels,
        initializer=self.global_config.final_init,
        name='output_projection',
    )(weighted_avg)

  def __call__(self, act, pair_mask):
    """Builds a module.

    Arguments:
      act: [num_seq, num_res, channels] activations tensor
      pair_mask: [num_seq, num_res] mask of non-padded regions in the tensor.
        Only used in inducing points attention currently.

    Returns:
      Result of the self-attention operation.
    """
    assert len(act.shape) == 3
    assert len(pair_mask.shape) == 2

    pair_mask = jnp.swapaxes(pair_mask, -1, -2)
    act = hm.LayerNorm(name='act_norm')(act)

    nonbatched_bias = hm.Linear(
        self.config.num_head, use_bias=False, name='pair_bias_projection'
    )(act)
    nonbatched_bias = jnp.transpose(nonbatched_bias, [2, 0, 1])

    num_residues = act.shape[0]

    chunk_size = get_shard_size(
        num_residues, self.global_config.pair_attention_chunk_size
    )

    if self.transpose:
      act = jnp.swapaxes(act, -2, -3)

    pair_mask = pair_mask[:, None, None, :].astype(jnp.bool_)

    act = mapping.inference_subbatch(
        self._attention,
        chunk_size,
        batched_args=[act, pair_mask],
        nonbatched_args=[nonbatched_bias],
    )

    if self.transpose:
      act = jnp.swapaxes(act, -2, -3)

    return act


class TriangleMultiplication(hk.Module):
  """Triangle Multiplication."""

  class Config(base_config.BaseConfig):
    equation: Literal['ikc,jkc->ijc', 'kjc,kic->ijc']
    use_glu_kernel: bool = True

  def __init__(
      self, config: Config, global_config: model_config.GlobalConfig, *, name
  ):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self, act, mask):
    """Applies Module.

    Args:
      act: The activation.
      mask: The mask.

    Returns:
      Outputs, should have same shape/type as output_act
    """
    num_residues = act.shape[0]

    chunk_size = get_shard_size(
        num_residues, self.global_config.pair_attention_chunk_size
    )

    return mapping.inference_subbatch(
        self._apply,
        chunk_size,
        batched_args=[act, mask],
        nonbatched_args=[],
    )

  def _apply(self, act, mask):
    """Applies Module.

    Args:
      act: The activation.
      mask: The mask.

    Returns:
      Outputs, should have same shape/type as output_act
    """
    mask = mask[None, ...]
    num_channels = act.shape[-1]
    equation = {
        'ikc,jkc->ijc': 'cik,cjk->cij',
        'kjc,kic->ijc': 'ckj,cki->cij',
    }[self.config.equation]

    act = hm.LayerNorm(name='left_norm_input')(act)
    input_act = act

    if self.global_config.use_glu_kernel:
      weights_projection, _ = hm.haiku_linear_get_params(
          act, num_output=num_channels * 2, name='projection'
      )
      weights_gate, _ = hm.haiku_linear_get_params(
          act,
          num_output=num_channels * 2,
          initializer=self.global_config.final_init,
          name='gate',
      )
      weights_glu = jnp.stack([weights_gate, weights_projection], axis=1)

      projection = tokamax.gated_linear_unit(
          act, weights_glu, activation=jax.nn.sigmoid
      )
      projection = jnp.transpose(projection, (2, 0, 1))
      projection *= mask
    else:
      projection = hm.Linear(num_channels * 2, name='projection')(act)
      projection = jnp.transpose(projection, (2, 0, 1))
      projection *= mask

      gate = hm.Linear(
          num_channels * 2,
          name='gate',
          bias_init=1.0,
          initializer=self.global_config.final_init,
      )(act)
      gate = jnp.transpose(gate, (2, 0, 1))
      projection *= jax.nn.sigmoid(gate)

    projection = projection.reshape(num_channels, 2, *projection.shape[1:])
    a, b = jnp.split(projection, 2, axis=1)
    a, b = jnp.squeeze(a, axis=1), jnp.squeeze(b, axis=1)
    act = jnp.einsum(equation, a, b)
    act = hm.LayerNorm(name='center_norm', axis=0, param_axis=0)(act)

    act = jnp.transpose(act, (1, 2, 0))
    act = hm.Linear(
        num_channels,
        initializer=self.global_config.final_init,
        name='output_projection',
    )(act)

    gate_out = hm.Linear(
        num_channels,
        name='gating_linear',
        bias_init=1.0,
        initializer=self.global_config.final_init,
    )(input_act)
    act *= jax.nn.sigmoid(gate_out)

    return act


class OuterProductMean(hk.Module):
  """Computed mean outer product."""

  class Config(base_config.BaseConfig):
    chunk_size: int = 128
    num_outer_channel: int = 32

  def __init__(
      self,
      config: Config,
      global_config: model_config.GlobalConfig,
      num_output_channel,
      *,
      name,
  ):
    super().__init__(name=name)
    self.global_config = global_config
    self.config = config
    self.num_output_channel = num_output_channel

  def __call__(self, act, mask):
    mask = mask[..., None]
    act = hm.LayerNorm(name='layer_norm_input')(act)

    left_act = mask * hm.Linear(
        self.config.num_outer_channel,
        initializer='linear',
        name='left_projection',
    )(act)

    right_act = mask * hm.Linear(
        self.config.num_outer_channel,
        initializer='linear',
        name='right_projection',
    )(act)

    if self.global_config.final_init == 'zeros':
      w_init = hk.initializers.Constant(0.0)
    else:
      w_init = hk.initializers.VarianceScaling(scale=2.0, mode='fan_in')

    output_w = hk.get_parameter(
        'output_w',
        shape=(
            self.config.num_outer_channel,
            self.config.num_outer_channel,
            self.num_output_channel,
        ),
        dtype=act.dtype,
        init=w_init,
    )
    output_b = hk.get_parameter(
        'output_b',
        shape=(self.num_output_channel,),
        dtype=act.dtype,
        init=hk.initializers.Constant(0.0),
    )

    def compute_chunk(left_act):
      # Make sure that the 'b' dimension is the most minor batch like dimension
      # so it will be treated as the real batch by XLA (both during the forward
      # and the backward pass)
      left_act = jnp.transpose(left_act, [0, 2, 1])
      act = jnp.einsum('acb,ade->dceb', left_act, right_act)
      act = jnp.einsum('dceb,cef->dbf', act, output_w) + output_b
      return jnp.transpose(act, [1, 0, 2])

    act = mapping.inference_subbatch(
        compute_chunk,
        self.config.chunk_size,
        batched_args=[left_act],
        nonbatched_args=[],
        input_subbatch_dim=1,
        output_subbatch_dim=0,
    )

    epsilon = 1e-3
    norm = jnp.einsum('abc,adc->bdc', mask, mask)
    return act / (epsilon + norm)


class PairFormerIteration(hk.Module):
  """Single Iteration of Pair Former."""

  class Config(base_config.BaseConfig):
    """Config for PairFormerIteration."""

    num_layer: int
    pair_attention: GridSelfAttention.Config = base_config.autocreate()
    pair_transition: TransitionBlock.Config = base_config.autocreate()
    single_attention: diffusion_transformer.SelfAttentionConfig | None = None
    single_transition: TransitionBlock.Config | None = None
    triangle_multiplication_incoming: TriangleMultiplication.Config = (
        base_config.autocreate(equation='kjc,kic->ijc')
    )
    triangle_multiplication_outgoing: TriangleMultiplication.Config = (
        base_config.autocreate(equation='ikc,jkc->ijc')
    )
    shard_transition_blocks: bool = True

  def __init__(
      self,
      config: Config,
      global_config: model_config.GlobalConfig,
      with_single=False,
      *,
      name,
  ):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config
    self.with_single = with_single

  def __call__(
      self,
      act,
      pair_mask,
      single_act=None,
      seq_mask=None,
  ):
    """Build a single iteration of the pair former.

    Args:
      act: [num_res, num_res, num_channel] Input pairwise activations.
      pair_mask: [num_res, num_res] padding mask.
      single_act: [num_res, single_channel] Single Input activations, optional
      seq_mask: [num_res] Sequence Mask, optional.

    Returns:
      [num_res, num_res, num_channel] tensor of activations.
    """

    num_residues = act.shape[0]

    act += TriangleMultiplication(
        self.config.triangle_multiplication_outgoing,
        self.global_config,
        name='triangle_multiplication_outgoing',
    )(act, pair_mask)

    act += TriangleMultiplication(
        self.config.triangle_multiplication_incoming,
        self.global_config,
        name='triangle_multiplication_incoming',
    )(act, pair_mask)

    act += GridSelfAttention(
        self.config.pair_attention,
        self.global_config,
        name='pair_attention1',
        transpose=False,
    )(act, pair_mask)

    act += GridSelfAttention(
        self.config.pair_attention,
        self.global_config,
        name='pair_attention2',
        transpose=True,
    )(act, pair_mask)

    transition_block = TransitionBlock(
        self.config.pair_transition, self.global_config, name='pair_transition'
    )
    if self.config.shard_transition_blocks:
      transition_block = mapping.sharded_apply(
          transition_block,
          get_shard_size(
              num_residues, self.global_config.pair_transition_shard_spec
          ),
      )
    act += transition_block(act)

    if self.with_single:
      assert self.config.single_attention is not None
      pair_logits = hm.Linear(
          self.config.single_attention.num_head,
          name='single_pair_logits_projection',
      )(hm.LayerNorm(name='single_pair_logits_norm')(act))

      pair_logits = jnp.transpose(pair_logits, [2, 0, 1])

      single_act += diffusion_transformer.self_attention(
          single_act,
          seq_mask,
          pair_logits=pair_logits,
          config=self.config.single_attention,
          global_config=self.global_config,
          name='single_attention_',
      )

      single_act += TransitionBlock(
          self.config.single_transition,
          self.global_config,
          name='single_transition',
      )(single_act, broadcast_dim=None)

      return act, single_act
    else:
      return act


class EvoformerIteration(hk.Module):
  """Single Iteration of Evoformer Main Stack."""

  class Config(base_config.BaseConfig):
    """Configuration for EvoformerIteration."""

    num_layer: int = 4
    msa_attention: MSAAttention.Config = base_config.autocreate()
    outer_product_mean: OuterProductMean.Config = base_config.autocreate()
    msa_transition: TransitionBlock.Config = base_config.autocreate()
    pair_attention: GridSelfAttention.Config = base_config.autocreate()
    pair_transition: TransitionBlock.Config = base_config.autocreate()
    triangle_multiplication_incoming: TriangleMultiplication.Config = (
        base_config.autocreate(equation='kjc,kic->ijc')
    )
    triangle_multiplication_outgoing: TriangleMultiplication.Config = (
        base_config.autocreate(equation='ikc,jkc->ijc')
    )
    shard_transition_blocks: bool = True

  def __init__(
      self,
      config: Config,
      global_config: model_config.GlobalConfig,
      name='evoformer_iteration',
  ):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self, activations, masks):

    msa_act, pair_act = activations['msa'], activations['pair']

    num_residues = pair_act.shape[0]

    msa_mask, pair_mask = masks['msa'], masks['pair']

    pair_act += OuterProductMean(
        config=self.config.outer_product_mean,
        global_config=self.global_config,
        num_output_channel=int(pair_act.shape[-1]),
        name='outer_product_mean',
    )(msa_act, msa_mask)

    msa_act += MSAAttention(
        self.config.msa_attention, self.global_config, name='msa_attention1'
    )(msa_act, msa_mask, pair_act=pair_act)

    msa_act += TransitionBlock(
        self.config.msa_transition, self.global_config, name='msa_transition'
    )(msa_act)

    pair_act += TriangleMultiplication(
        self.config.triangle_multiplication_outgoing,
        self.global_config,
        name='triangle_multiplication_outgoing',
    )(pair_act, pair_mask)

    pair_act += TriangleMultiplication(
        self.config.triangle_multiplication_incoming,
        self.global_config,
        name='triangle_multiplication_incoming',
    )(pair_act, pair_mask)

    pair_act += GridSelfAttention(
        self.config.pair_attention,
        self.global_config,
        name='pair_attention1',
        transpose=False,
    )(pair_act, pair_mask)

    pair_act += GridSelfAttention(
        self.config.pair_attention,
        self.global_config,
        name='pair_attention2',
        transpose=True,
    )(pair_act, pair_mask)

    transition_block = TransitionBlock(
        self.config.pair_transition, self.global_config, name='pair_transition'
    )
    if self.config.shard_transition_blocks:
      transition_block = mapping.sharded_apply(
          transition_block,
          get_shard_size(
              num_residues, self.global_config.pair_transition_shard_spec
          ),
      )
    pair_act += transition_block(pair_act)

    return {'msa': msa_act, 'pair': pair_act}
