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

"""Rigid3Array Transformations represented by a Matrix and a Vector."""

from typing import Any, Final, Self, TypeAlias

from alphafold3.jax.geometry import rotation_matrix
from alphafold3.jax.geometry import struct_of_array
from alphafold3.jax.geometry import utils
from alphafold3.jax.geometry import vector
import jax
import jax.numpy as jnp


Float: TypeAlias = float | jnp.ndarray

VERSION: Final[str] = '0.1'


# Disabling name in pylint, since the relevant variable in math are typically
# referred to as X, Y in mathematical literature.
def _compute_covariance_matrix(
    row_values: vector.Vec3Array,
    col_values: vector.Vec3Array,
    weights: jnp.ndarray,
    epsilon=1e-6,
) -> jnp.ndarray:
  """Compute covariance matrix.

  The quantity computes is
  cov_xy = weighted_avg_i(row_values[i, x] col_values[j, y]).
  Here x and y run over the xyz coordinates.
  This is used to construct frames when aligning points.

  Args:
    row_values: Values used for rows of covariance matrix, shape [..., n_point]
    col_values: Values used for columns of covariance matrix, shape [...,
      n_point]
    weights: weights to weight points by, shape broacastable to [...]
    epsilon: small value to add to denominator to avoid Nan's when all weights
      are 0.

  Returns:
    Covariance Matrix as [..., 3, 3] array.
  """
  weights = jnp.asarray(weights)
  weights = jnp.broadcast_to(weights, row_values.shape)

  out = []

  normalized_weights = weights / (weights.sum(axis=-1, keepdims=True) + epsilon)

  weighted_average = lambda x: jnp.sum(normalized_weights * x, axis=-1)

  out.append(
      jnp.stack(
          (
              weighted_average(row_values.x * col_values.x),
              weighted_average(row_values.x * col_values.y),
              weighted_average(row_values.x * col_values.z),
          ),
          axis=-1,
      )
  )

  out.append(
      jnp.stack(
          (
              weighted_average(row_values.y * col_values.x),
              weighted_average(row_values.y * col_values.y),
              weighted_average(row_values.y * col_values.z),
          ),
          axis=-1,
      )
  )

  out.append(
      jnp.stack(
          (
              weighted_average(row_values.z * col_values.x),
              weighted_average(row_values.z * col_values.y),
              weighted_average(row_values.z * col_values.z),
          ),
          axis=-1,
      )
  )

  return jnp.stack(out, axis=-2)


@struct_of_array.StructOfArray(same_dtype=True)
class Rigid3Array:
  """Rigid Transformation, i.e. element of special euclidean group."""

  rotation: rotation_matrix.Rot3Array
  translation: vector.Vec3Array

  def __matmul__(self, other: Self) -> Self:
    new_rotation = self.rotation @ other.rotation
    new_translation = self.apply_to_point(other.translation)
    return Rigid3Array(new_rotation, new_translation)

  def inverse(self) -> Self:
    """Return Rigid3Array corresponding to inverse transform."""
    inv_rotation = self.rotation.inverse()
    inv_translation = inv_rotation.apply_to_point(-self.translation)
    return Rigid3Array(inv_rotation, inv_translation)

  def apply_to_point(self, point: vector.Vec3Array) -> vector.Vec3Array:
    """Apply Rigid3Array transform to point."""
    return self.rotation.apply_to_point(point) + self.translation

  def apply_inverse_to_point(self, point: vector.Vec3Array) -> vector.Vec3Array:
    """Apply inverse Rigid3Array transform to point."""
    new_point = point - self.translation
    return self.rotation.apply_inverse_to_point(new_point)

  def compose_rotation(self, other_rotation: rotation_matrix.Rot3Array) -> Self:
    rot = self.rotation @ other_rotation
    trans = jax.tree.map(
        lambda x: jnp.broadcast_to(x, rot.shape), self.translation
    )
    return Rigid3Array(rot, trans)

  @classmethod
  def identity(cls, shape: Any, dtype: jnp.dtype = jnp.float32) -> Self:
    """Return identity Rigid3Array of given shape."""
    return cls(
        rotation_matrix.Rot3Array.identity(shape, dtype=dtype),
        vector.Vec3Array.zeros(shape, dtype=dtype),
    )  # pytype: disable=wrong-arg-count  # trace-all-classes

  def scale_translation(self, factor: Float) -> Self:
    """Scale translation in Rigid3Array by 'factor'."""
    return Rigid3Array(self.rotation, self.translation * factor)

  def to_array(self):
    rot_array = self.rotation.to_array()
    vec_array = self.translation.to_array()
    return jnp.concatenate([rot_array, vec_array[..., None]], axis=-1)

  @classmethod
  def from_array(cls, array):
    rot = rotation_matrix.Rot3Array.from_array(array[..., :3])
    vec = vector.Vec3Array.from_array(array[..., -1])
    return cls(rot, vec)  # pytype: disable=wrong-arg-count  # trace-all-classes

  @classmethod
  def from_array4x4(cls, array: jnp.ndarray) -> Self:
    """Construct Rigid3Array from homogeneous 4x4 array."""
    if array.shape[-2:] != (4, 4):
      raise ValueError(f'array.shape({array.shape}) must be [..., 4, 4]')
    rotation = rotation_matrix.Rot3Array(
        *(array[..., 0, 0], array[..., 0, 1], array[..., 0, 2]),
        *(array[..., 1, 0], array[..., 1, 1], array[..., 1, 2]),
        *(array[..., 2, 0], array[..., 2, 1], array[..., 2, 2]),
    )
    translation = vector.Vec3Array(
        array[..., 0, 3], array[..., 1, 3], array[..., 2, 3]
    )
    return cls(rotation, translation)  # pytype: disable=wrong-arg-count  # trace-all-classes

  @classmethod
  def from_point_alignment(
      cls,
      points_to: vector.Vec3Array,
      points_from: vector.Vec3Array,
      weights: Float | None = None,
      epsilon: float = 1e-6,
  ) -> Self:
    """Constructs Rigid3Array by finding transform aligning points.

    This constructs the optimal Rigid Transform taking points_from to the
    arrangement closest to points_to.

    Args:
      points_to: Points to align to.
      points_from: Points to align from.
      weights: weights for points.
      epsilon: epsilon used to regularize covariance matrix.

    Returns:
      Rigid Transform.
    """
    if weights is None:
      weights = 1.0

    def compute_center(value):
      return utils.weighted_mean(value=value, weights=weights, axis=-1)

    points_to_center = jax.tree.map(compute_center, points_to)
    points_from_center = jax.tree.map(compute_center, points_from)
    centered_points_to = points_to - points_to_center[..., None]
    centered_points_from = points_from - points_from_center[..., None]
    cov_mat = _compute_covariance_matrix(
        centered_points_to,
        centered_points_from,
        weights=weights,
        epsilon=epsilon,
    )
    rots = rotation_matrix.Rot3Array.from_svd(
        jnp.reshape(cov_mat, cov_mat.shape[:-2] + (9,))
    )

    translations = points_to_center - rots.apply_to_point(points_from_center)

    return cls(rots, translations)  # pytype: disable=wrong-arg-count  # trace-all-classes

  def __getstate__(self):
    return (VERSION, (self.rotation, self.translation))

  def __setstate__(self, state):
    version, (rot, trans) = state
    del version
    object.__setattr__(self, 'rotation', rot)
    object.__setattr__(self, 'translation', trans)
