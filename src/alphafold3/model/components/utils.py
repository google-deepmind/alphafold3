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

"""Utility functions for training AlphaFold and similar models."""

from collections import abc
from collections.abc import Sequence
import contextlib
import numbers

from alphafold3.model import features
import haiku as hk
import jax.numpy as jnp
import numpy as np

VALID_DTYPES = [np.float32, np.float64, np.int8, np.int32, np.int64, bool]


def remove_invalidly_typed_feats(
    batch: features.BatchDict,
) -> features.BatchDict:
  """Remove features of types we don't want to send to the TPU e.g. strings."""
  return {
      k: v
      for k, v in batch.items()
      if hasattr(v, 'dtype') and v.dtype in VALID_DTYPES
  }


def bfloat16_getter(next_getter, value, context):
  """Ensures that a bfloat16 parameter is provided by casting if necessary."""
  if context.original_dtype == jnp.bfloat16:
    if value.dtype != jnp.bfloat16:
      value = value.astype(jnp.bfloat16)
  return next_getter(value)


@contextlib.contextmanager
def bfloat16_context():
  with hk.custom_getter(bfloat16_getter):
    yield


def mask_mean(
    mask: jnp.ndarray,
    value: jnp.ndarray,
    axis: int | Sequence[int] | None = None,
    keepdims: bool = False,
    eps: float = 1e-10,
) -> jnp.ndarray:
  """Masked mean."""

  mask_shape = mask.shape
  value_shape = value.shape

  if len(mask_shape) != len(value_shape):
    raise ValueError(f'Incompatible shapes: {mask_shape=}, {value_shape=}')

  if isinstance(axis, numbers.Integral):
    axis = [axis]
  elif axis is None:
    axis = list(range(len(mask_shape)))
  assert isinstance(axis, abc.Sequence)

  broadcast_factor = 1.0
  for axis_ in axis:
    value_size = value_shape[axis_]
    mask_size = mask_shape[axis_]
    if mask_size == 1:
      broadcast_factor *= value_size
    else:
      if mask_size != value_size:
        raise ValueError(f'Incompatible shapes: {mask_shape=}, {value_shape=}')

  return jnp.sum(mask * value, keepdims=keepdims, axis=axis) / (
      jnp.maximum(
          jnp.sum(mask, keepdims=keepdims, axis=axis) * broadcast_factor, eps
      )
  )
