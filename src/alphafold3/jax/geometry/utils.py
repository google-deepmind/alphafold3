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

"""Utils for geometry library."""

from collections.abc import Iterable
import numbers

import jax.numpy as jnp


def unstack(value: jnp.ndarray, axis: int = -1) -> list[jnp.ndarray]:
  return [
      jnp.squeeze(v, axis=axis)
      for v in jnp.split(value, value.shape[axis], axis=axis)
  ]


def angdiff(alpha: jnp.ndarray, beta: jnp.ndarray) -> jnp.ndarray:
  """Compute absolute difference between two angles."""
  d = alpha - beta
  d = (d + jnp.pi) % (2 * jnp.pi) - jnp.pi
  return d


def weighted_mean(
    *,
    weights: jnp.ndarray,
    value: jnp.ndarray,
    axis: int | Iterable[int] | None = None,
    eps: float = 1e-10,
) -> jnp.ndarray:
  """Computes weighted mean in a safe way that avoids NaNs.

  This is equivalent to jnp.average for the case eps=0.0, but adds a small
  constant to the denominator of the weighted average to avoid NaNs.
  'weights' should be broadcastable to the shape of value.

  Args:
    weights: Weights to weight value by.
    value: Values to average
    axis: Axes to average over.
    eps: Epsilon to add to the denominator.

  Returns:
    Weighted average.
  """

  weights = jnp.asarray(weights, dtype=value.dtype)
  weights = jnp.broadcast_to(weights, value.shape)

  weights_shape = weights.shape

  if isinstance(axis, numbers.Integral):
    axis = [axis]
  elif axis is None:
    axis = list(range(len(weights_shape)))

  return jnp.sum(weights * value, axis=axis) / (
      jnp.sum(weights, axis=axis) + eps
  )
