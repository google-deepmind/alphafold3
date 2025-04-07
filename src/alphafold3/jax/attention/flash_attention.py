# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Triton FlashAttention implementation."""

import dataclasses
import functools

from alphafold3.jax.attention import attention_base as base
from alphafold3.jax.common import array_view
from alphafold3.jax.common import precision as precision_lib
from alphafold3.jax.common import triton_utils
import jax
import jax.numpy as jnp
import jax_triton as jt
import jaxtyping
from jaxtyping import Array, Bool, Float, Int  # pylint: disable=g-multiple-import,g-importing-member
import triton
import triton.language as tl
import typeguard


@triton.jit
def _fwd_kernel_inner(
    start_loop,
    end_loop,
    q,
    span_q,
    k_block_ptr,
    v_block_ptr,
    bias_block_ptr,
    mask_block_ptr,
    k_start,
    k_end,
    seq_len_k,
    acc,
    m_i,
    l_i,
    bias_advance: tl.constexpr,
    mask_advance: tl.constexpr,
    is_causal: tl.constexpr,
    use_attention_mask: tl.constexpr,
    use_k_start: tl.constexpr,
    use_k_end: tl.constexpr,
    use_bias: tl.constexpr,
    block_k: tl.constexpr,
    use_mask_k: tl.constexpr,
    k_boundary_check: tl.constexpr,
    v_boundary_check: tl.constexpr,
    dot_fn_qk: tl.constexpr,
    dot_fn_kv: tl.constexpr,
):
  """Triton MHA forward kernel's inner loop."""

  for start_k in range(start_loop, end_loop, block_k):
    start_k = tl.multiple_of(start_k, block_k)
    span_k = start_k + tl.arange(0, block_k)

    k = tl.load(
        k_block_ptr,
        boundary_check=k_boundary_check,
        padding_option="zero" if len(k_boundary_check.value) else "",
    )
    v = tl.load(
        v_block_ptr,
        boundary_check=v_boundary_check,
        padding_option="zero" if len(v_boundary_check.value) else "",
    )

    if use_bias:
      bias = tl.load(bias_block_ptr)

    qk = dot_fn_qk(q.to(k.dtype), k)  # [block_q, block_k]

    if use_bias:
      # Prevent dot accumulating into the bias tensor. It appears that Triton
      # doesn't pipeline the bias load as it does the `k` load, so the bias load
      # blocks the matmul if the add is merged.
      qk = qk.to(tl.uint32, bitcast=True) & 0xFFFFFFFF
      qk = qk.to(tl.float32, bitcast=True)
      qk += bias

    if use_attention_mask | use_k_start | use_k_end:
      mask_value = float(jnp.finfo(jnp.float32).min)

    if use_attention_mask:
      mask = tl.load(mask_block_ptr)
      qk = tl.where(mask, qk, mask_value)

    if use_k_start:
      # This check is there to work around a triton compiler bug, but it
      # shouldn't be strictly needed.
      if tl.sum(k_start) != 0:
        qk = tl.where(k_start[:, None] <= span_k[None, :], qk, mask_value)
    if is_causal:
      qk = tl.where(span_q[:, None] >= span_k[None, :], qk, float("-inf"))
    elif use_k_end:
      # When called with k_end and is_causal=True, the causal mask gets folded
      # into k_end and is_causal is set to False.
      qk = tl.where(k_end[:, None] > span_k[None, :], qk, mask_value)

    if use_mask_k:
      qk = tl.where((span_k < seq_len_k)[None, :], qk, float("-inf"))

    m_ij = tl.maximum(m_i, tl.max(qk, axis=1))  # Shape [block_q].
    p = tl.exp(qk - m_ij[:, None])  # Shape [block_q, block_k].
    alpha = tl.exp(m_i - m_ij)
    m_i = m_ij
    acc *= alpha[:, None]
    l_i *= alpha
    l_i += tl.sum(p, axis=1)

    # Add the new block of attention weights.
    acc += dot_fn_kv(p.to(v.dtype), v)

    k_block_ptr = tl.advance(k_block_ptr, (0, block_k))
    v_block_ptr = tl.advance(v_block_ptr, (block_k, 0))
    bias_block_ptr = tl.advance(bias_block_ptr, bias_advance.value)
    mask_block_ptr = tl.advance(mask_block_ptr, mask_advance.value)

  return (
      k_block_ptr,
      v_block_ptr,
      bias_block_ptr,
      mask_block_ptr,
      acc,
      m_i,
      l_i,
  )


# Based on Algorithm 1 of https://arxiv.org/abs/2205.14135.
# Inspired by the official Triton tutorial implementation
# https://triton-lang.org/main/getting-started/tutorials/06-fused-attention.html
@triton.jit
def _fwd_kernel(
    # Input arrays.
    q_ptr,
    k_ptr,
    v_ptr,
    bias_ptr,
    mask_ptr,
    k_start_ptr,
    k_end_ptr,
    # Scalar inputs.
    q_offset,
    k_offset,
    v_offset,
    q_stride_b,
    q_stride_s,
    q_stride_h,
    q_stride_d,
    k_stride_b,
    k_stride_s,
    k_stride_h,
    k_stride_d,
    v_stride_b,
    v_stride_s,
    v_stride_h,
    v_stride_d,
    bias_stride_b,
    bias_stride_h,
    bias_stride_sq,
    bias_stride_sk,
    mask_stride_b,
    mask_stride_h,
    mask_stride_sq,
    mask_stride_sk,
    k_start_stride_b,
    k_start_stride_h,
    k_start_stride_sq,
    k_end_stride_b,
    k_end_stride_h,
    k_end_stride_sq,
    o_stride_b,
    o_stride_s,
    o_stride_h,
    o_stride_d,
    num_heads_q,
    num_heads_k,
    seq_len_q,
    seq_len_k,
    # Output arrays.
    o_ptr,
    # Compile-time constants.
    is_causal: tl.constexpr,
    use_attention_mask: tl.constexpr,
    use_k_start: tl.constexpr,
    use_k_end: tl.constexpr,
    use_bias: tl.constexpr,
    sm_scale: tl.constexpr,
    block_q: tl.constexpr,
    block_k: tl.constexpr,
    head_dim: tl.constexpr,
    use_mask_q: tl.constexpr,
    use_mask_k: tl.constexpr,
    bias_bcast_sq: tl.constexpr,
    mask_bcast_sq: tl.constexpr,
    dot_fn_qk: tl.constexpr,
    dot_fn_kv: tl.constexpr,
):
  """Triton MHA forward kernel."""
  # pytype: disable=annotation-type-mismatch,unsupported-operands
  block_d: tl.constexpr = jt.utils.next_power_of_2(head_dim.value)

  # Each thread block processes one batch element (b) and one head (h).
  start_q = tl.program_id(1) * block_q
  off_h = tl.program_id(0)  # int in [0, num_heads_o).
  off_b = tl.program_id(2)  # int in [0, batch_size)

  off_h_k = off_h // (num_heads_q // num_heads_k)

  q_ptr += off_h * q_stride_h + off_b * q_stride_b + q_offset
  k_ptr += off_h_k * k_stride_h + off_b * k_stride_b + k_offset
  v_ptr += off_h_k * v_stride_h + off_b * v_stride_b + v_offset
  o_ptr += off_h * o_stride_h + off_b * o_stride_b

  if use_bias:
    bias_ptr += off_b * bias_stride_b + off_h * bias_stride_h
  if use_attention_mask:
    mask_ptr += off_b * mask_stride_b + off_h * mask_stride_h
  if use_k_start:
    k_start_ptr += off_b * k_start_stride_b + off_h * k_start_stride_h
  if use_k_end:
    k_end_ptr += off_b * k_end_stride_b + off_h * k_end_stride_h

  q_block_ptr = tl.make_block_ptr(
      q_ptr,
      shape=(seq_len_q, head_dim),
      strides=(q_stride_s, q_stride_d),
      offsets=(start_q, 0),
      block_shape=(block_q, block_d),
      order=(1, 0),
  )
  k_block_ptr = tl.make_block_ptr(
      k_ptr,
      shape=(head_dim, seq_len_k),
      strides=(k_stride_d, k_stride_s),
      offsets=(0, 0),
      block_shape=(block_d, block_k),
      order=(0, 1),
  )
  v_block_ptr = tl.make_block_ptr(
      v_ptr,
      shape=(seq_len_k, head_dim),
      strides=(v_stride_s, v_stride_d),
      offsets=(0, 0),
      block_shape=(block_k, block_d),
      order=(1, 0),
  )

  q_boundary_check0: tl.constexpr = (0,) if use_mask_q else ()
  q_boundary_check1: tl.constexpr = (1,) if head_dim != block_d else ()
  q_boundary_check: tl.constexpr = q_boundary_check0 + q_boundary_check1
  q_padding_option: tl.constexpr = "zero" if len(q_boundary_check.value) else ""
  k_boundary_check: tl.constexpr = (0,) if head_dim != block_d else ()
  v_boundary_check: tl.constexpr = (0,) if use_mask_k else ()

  # If broadcasting in a given dim, use a 1D block (observed to be faster).
  bias_start_dim: tl.constexpr = 1 if bias_bcast_sq else 0
  bias_block_ptr = tl.make_block_ptr(
      bias_ptr,
      shape=(seq_len_q, seq_len_k)[bias_start_dim:],
      strides=(bias_stride_sq, bias_stride_sk)[bias_start_dim:],
      offsets=(start_q, 0)[bias_start_dim:],
      block_shape=(block_q, block_k)[bias_start_dim:],
      order=(1, 0)[bias_start_dim:],
  )
  bias_advance: tl.constexpr = (0, block_k)[bias_start_dim:]

  mask_start_dim: tl.constexpr = 1 if mask_bcast_sq else 0
  mask_block_ptr = tl.make_block_ptr(
      mask_ptr,
      shape=(seq_len_q, seq_len_k)[mask_start_dim:],
      strides=(mask_stride_sq, mask_stride_sk)[mask_start_dim:],
      offsets=(start_q, 0)[mask_start_dim:],
      block_shape=(block_q, block_k)[mask_start_dim:],
      order=(1, 0)[mask_start_dim:],
  )
  mask_advance: tl.constexpr = (0, block_k)[mask_start_dim:]

  k_start_block_ptr = tl.make_block_ptr(
      k_start_ptr,
      shape=(seq_len_q,),
      strides=(k_start_stride_sq,),
      offsets=(start_q,),
      block_shape=(block_q,),
      order=(0,),
  )
  k_end_block_ptr = tl.make_block_ptr(
      k_end_ptr,
      shape=(seq_len_q,),
      strides=(k_end_stride_sq,),
      offsets=(start_q,),
      block_shape=(block_q,),
      order=(0,),
  )
  # pytype: enable=annotation-type-mismatch,unsupported-operands

  # Each thread block processes a block of block_q queries.
  span_q = start_q + tl.arange(0, block_q)

  # m_i and l_i (see FlashAttention paper) are updated during the k,v loop.
  m_i = tl.full([block_q], float("-inf"), dtype=tl.float32)
  l_i = tl.zeros([block_q], dtype=tl.float32)
  # acc is the buffer where we accumulate the output on sram.
  acc = tl.zeros([block_q, block_d], dtype=tl.float32)

  # Load q: it will stay in smem throughout. Indices form a matrix because we
  # read, compute, and write all in 2d chunks. 1 element ~= 1 CUDA thread index.
  q = tl.load(
      q_block_ptr,
      boundary_check=q_boundary_check,
      padding_option=q_padding_option,
  )
  q *= sm_scale

  # In FlashAttention algorithm 1 there are 2 loops: slow over tiles of kv (size
  # (Bc == block_k here), and fast over blocks of q (size Br == block_q here).
  # Here we only loop over blocks of kv to process entire seq_len, the loop over
  # blocks of q is carried out by the grid.
  k_start = None
  if use_k_start:
    k_start = tl.load(k_start_block_ptr)
    start_loop = tl.maximum(tl.min(k_start), 0)
    blocks_to_skip = start_loop // block_k
    start_loop = block_k * blocks_to_skip  # Floor to multiple of block_k.
    for _ in range(blocks_to_skip):
      # Advance all block pointers to the first valid block.
      k_block_ptr = tl.advance(k_block_ptr, (0, block_k))
      v_block_ptr = tl.advance(v_block_ptr, (block_k, 0))
      bias_block_ptr = tl.advance(bias_block_ptr, bias_advance.value)
      mask_block_ptr = tl.advance(mask_block_ptr, mask_advance.value)
  else:
    start_loop = 0

  k_end = None
  if is_causal:
    end_loop = tl.minimum((start_q // block_k) * block_k, seq_len_k)
  elif use_k_end:
    k_end = tl.load(k_end_block_ptr)
    end_loop = tl.minimum(tl.max(k_end), seq_len_k)
  else:
    end_loop = seq_len_k

  (
      k_block_ptr,
      v_block_ptr,
      bias_block_ptr,
      mask_block_ptr,
      acc,
      m_i,
      l_i,
  ) = _fwd_kernel_inner(
      start_loop,
      end_loop,
      q,
      span_q,
      k_block_ptr,
      v_block_ptr,
      bias_block_ptr,
      mask_block_ptr,
      k_start,
      k_end,
      seq_len_k,
      acc,
      m_i,
      l_i,
      bias_advance,
      mask_advance,
      False,  # is_causal
      use_attention_mask,
      use_k_start,
      use_k_end,
      use_bias,
      block_k,
      use_mask_k,
      k_boundary_check,
      v_boundary_check,
      dot_fn_qk,
      dot_fn_kv,
  )

  if is_causal:
    tl.debug_barrier()  # Help compiler schedule loops independently.
    start_loop, end_loop = end_loop, tl.minimum(end_loop + block_k, seq_len_k)

    _, _, _, _, acc, _, l_i = _fwd_kernel_inner(
        start_loop,
        end_loop,
        q,
        span_q,
        k_block_ptr,
        v_block_ptr,
        bias_block_ptr,
        mask_block_ptr,
        k_start,
        k_end,
        seq_len_k,
        acc,
        m_i,
        l_i,
        bias_advance,
        mask_advance,
        True,  # is_causal
        use_attention_mask,
        use_k_start,
        use_k_end,
        use_bias,
        block_k,
        use_mask_k,
        k_boundary_check,
        v_boundary_check,
        dot_fn_qk,
        dot_fn_kv,
    )

  # It is possible that every value in a row was masked to f32 min or that the
  # main loop has been completely optimised out, and that `l_i` is `0` for that
  # row. Add epsilon value to avoid NaNs from `0 / 0`.
  l_i += float(jnp.finfo(jnp.float32).tiny)

  acc /= l_i[:, None]

  # Write output to dram.
  o_block_ptr = tl.make_block_ptr(
      o_ptr,
      shape=(seq_len_q, head_dim),
      strides=(o_stride_s, o_stride_d),
      offsets=(start_q, 0),
      block_shape=(block_q, block_d),
      order=(1, 0),
  )
  acc = acc.to(o_ptr.dtype.element_ty)
  tl.store(o_block_ptr, acc, boundary_check=q_boundary_check)


@jaxtyping.jaxtyped(typechecker=typeguard.typechecked)
def _fwd(
    q: Float[array_view.ArrayView, "*B T H D"],
    k: Float[array_view.ArrayView, "*B t h D"],
    v: Float[array_view.ArrayView, "*B t h D"],
    bias: Float[Array, "*#B #H #T #t"] | None,
    mask: Bool[Array, "*#B #H #T #t"] | None,
    k_start: Int[Array, "*#B #H #T"] | None,
    k_end: Int[Array, "*#B #H #T"] | None,
    *,
    logits_scale: float,
    is_causal: bool,
    q_k_dot_precision: precision_lib.DotPrecision,
    weights_v_dot_precision: precision_lib.DotPrecision,
) -> Float[Array, "*B T H D"]:
  """Forward pass of Triton FlashAttention."""

  orig_q_shape = q.shape
  q = q.collapse(0, -3, allow_copy=True)
  batch_size, seq_len_q, num_heads_q, head_dim = q.shape
  *_, seq_len_k, num_heads_kv, _ = k.shape
  # Maybe broadcast `k`/`v` heads dimension.
  kv_shape = (batch_size, seq_len_k, num_heads_kv, head_dim)
  k = k.collapse(0, -3, allow_copy=True).broadcast_to(kv_shape)
  v = v.collapse(0, -3, allow_copy=True).broadcast_to(kv_shape)

  def get_bias_mask_view(x, dtype):
    if x is None:
      x = jnp.array([], dtype=dtype)
      return array_view.ArrayView(x, shape=(0, 0, 0, 0), strides=(0, 0, 0, 0))

    shape = orig_q_shape[:-3] + (num_heads_q, seq_len_q, seq_len_k)
    return (
        array_view.ArrayView(x)
        .broadcast_to(shape)
        .collapse(0, -3, allow_copy=True)
    )

  bias = get_bias_mask_view(bias, dtype=q.dtype)
  mask = get_bias_mask_view(mask, dtype=jnp.bool_)

  def get_range_view(x, seq_len):
    if x is None:
      x = jnp.array([], dtype=jnp.int32)
      return array_view.ArrayView(x, shape=(0, 0, 0), strides=(0, 0, 0))

    shape = orig_q_shape[:-3] + (num_heads_q, seq_len)
    return (
        array_view.ArrayView(x)
        .broadcast_to(shape)
        .collapse(0, -2, allow_copy=True)
    )

  k_start = get_range_view(k_start, seq_len_q)
  k_end = get_range_view(k_end, seq_len_q)

  block_q = 64
  block_k = 64

  return jt.triton_call(
      q.base,
      k.base,
      v.base,
      bias.base,
      mask.base,
      k_start.base,
      k_end.base,
      q.offset,
      k.offset,
      v.offset,
      *q.strides,
      *k.strides,
      *v.strides,
      *bias.strides,
      *mask.strides,
      k_start.strides,
      k_end.strides,
      *jt.utils.strides_from_shape(q.shape),  # out strides.
      num_heads_q,
      num_heads_kv,
      seq_len_q,
      seq_len_k,
      kernel=_fwd_kernel,
      name="triton_flash_attention",
      out_shape=jax.ShapeDtypeStruct(q.shape, q.dtype),
      grid=(num_heads_q, triton.cdiv(seq_len_q, block_q), batch_size),
      num_stages=2,
      num_warps=4,
      is_causal=is_causal,
      use_attention_mask=(mask.size != 0),
      use_k_start=(k_start.size != 0),
      use_k_end=(k_end.size != 0),
      use_bias=(bias.size != 0),
      sm_scale=logits_scale,
      block_q=block_q,
      block_k=block_k,
      head_dim=head_dim,
      use_mask_q=(seq_len_q % block_q != 0),
      use_mask_k=(seq_len_k % block_q != 0),
      bias_bcast_sq=(bias.strides[-2] == 0),
      mask_bcast_sq=(mask.strides[-2] == 0),
      dot_fn_qk=triton_utils.get_tl_dot_fn(q_k_dot_precision),
      dot_fn_kv=triton_utils.get_tl_dot_fn(weights_v_dot_precision),
  ).reshape(orig_q_shape)


def _as_batched_array_view(x, axis_size):
  batched_shape = (axis_size,) + x.shape
  batched_strides = (x.base.size // axis_size,) + x.strides
  return dataclasses.replace(x, shape=batched_shape, strides=batched_strides)


def _fwd_vmap_rule(
    axis_size, in_batched, *args, fn: jax.custom_batching.custom_vmap
):
  """`vmap` rule for Triton FlashAttention forward op."""
  q, k, v, bias, mask, k_start, k_end = args
  (
      q_batched,
      k_batched,
      v_batched,
      bias_batched,
      mask_batched,
      k_start_batched,
      k_end_batched,
  ) = in_batched

  if q_batched.base:
    q = _as_batched_array_view(q, axis_size)
  if k_batched.base:
    k = _as_batched_array_view(k, axis_size)
  if v_batched.base:
    v = _as_batched_array_view(v, axis_size)

  # Triton op requires `q`, `k`, `v` batch dims to be identical.
  if q_batched.base and k_batched.base and v_batched.base:
    if bias is not None and not bias_batched:
      bias = jax.lax.broadcast_to_rank(bias, bias.ndim + 1)
    if mask is not None and not mask_batched:
      mask = jax.lax.broadcast_to_rank(mask, mask.ndim + 1)
    if k_start is not None and not k_start_batched:
      k_start = jax.lax.broadcast_to_rank(k_start, k_start.ndim + 1)
    if k_end is not None and not k_end_batched:
      k_end = jax.lax.broadcast_to_rank(k_end, k_end.ndim + 1)
    out = fn(q, k, v, bias, mask, k_start, k_end)
    out_batched = True
    return out, out_batched

  # Fallback to sequential loop.
  q, k, v = map(jnp.asarray, (q, k, v))
  in_batched = [
      q_batched.base,
      k_batched.base,
      v_batched.base,
      bias_batched,
      mask_batched,
      k_start_batched,
      k_end_batched,
  ]

  def f(q, k, v, *args, **kwargs):
    q, k, v = map(array_view.ArrayView, (q, k, v))
    return fn.fun(q, k, v, *args, **kwargs)

  sequential_vmap = jax.custom_batching.sequential_vmap(f)
  return sequential_vmap.vmap_rule(axis_size, in_batched, q, k, v, *args[3:])


def _decompose_mask(mask, q, k, q_indices, k_indices):
  """Decomposes `mask` into a mask array, `is_causal`, `k_start` and `k_end`."""
  if mask is None:
    return None, False, None, None

  is_causal = False
  k_start = None
  k_end = None
  if q_indices is None and k_indices is None:
    mask, is_causal, k_start, k_end = mask.take("is_causal", "k_start", "k_end")
    if k_start is not None:
      k_start = jax.lax.broadcast_to_rank(k_start, 2)
    if k_end is not None:
      k_end = jax.lax.broadcast_to_rank(k_end, 2)
      if is_causal:  # Fold is_causal into k_end
        k_end = jnp.minimum(k_end, jnp.arange(1, q.shape[-3] + 1))
        is_causal = False

  q_len_or_indices = q.shape[-3] if q_indices is None else q_indices
  k_len_or_indices = k.shape[-3] if k_indices is None else k_indices
  return (
      mask.as_array(q_len_or_indices, k_len_or_indices),
      is_causal,
      k_start,
      k_end,
  )


@dataclasses.dataclass(frozen=True)
class TritonFlashAttention(base.DotProductAttention):
  """Triton FlashAttention implementation."""

  @jaxtyping.jaxtyped(typechecker=typeguard.typechecked)
  def _fwd(
      self,
      q: Float[array_view.ArrayView, "*B T H D"],
      k: Float[array_view.ArrayView, "*B t h D"],
      v: Float[array_view.ArrayView, "*B t h D"],
      bias: Float[Array, "*#B #H #T #t"] | None,
      *,
      q_k_dot_precision: precision_lib.DotPrecision,
      logits_dtype: jnp.dtype,
      logits_scale: float,
      mask: base.Mask | None,
      weights_v_dot_precision: precision_lib.DotPrecision,
      q_indices: Int[Array, "*#B #H T"] | None = None,
      k_indices: Int[Array, "*#B #H t"] | None = None,
  ) -> Float[Array, "*B T H D"]:
    if logits_dtype != jnp.float32:
      raise ValueError("`logits_dtype` must be float32.")

    kwargs = dict(
        logits_scale=logits_scale,
        q_k_dot_precision=q_k_dot_precision,
        weights_v_dot_precision=weights_v_dot_precision,
    )

    def attend_fwd(
        q,
        k,
        v,
        bias,
        mask_,
        q_indices,
        k_indices,
    ):

      mask, is_causal, k_start, k_end = _decompose_mask(
          mask_, q, k, q_indices, k_indices
      )

      fwd_closed_kwargs = dict(
          is_causal=is_causal,
          **kwargs,
      )

      fwd_closed = functools.partial(_fwd, **fwd_closed_kwargs)
      fwd_closed = jax.custom_batching.custom_vmap(fwd_closed)
      fwd_closed.def_vmap(functools.partial(_fwd_vmap_rule, fn=fwd_closed))

      return fwd_closed(q, k, v, bias, mask, k_start, k_end)

    return attend_fwd(
        q,
        k,
        v,
        bias,
        mask,
        q_indices,
        k_indices,
    )
