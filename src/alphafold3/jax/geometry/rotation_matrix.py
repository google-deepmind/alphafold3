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

"""Rot3Array Matrix Class."""

import dataclasses
from typing import Any, Final, Self

from alphafold3.jax.geometry import struct_of_array
from alphafold3.jax.geometry import utils
from alphafold3.jax.geometry import vector
import jax
import jax.numpy as jnp
import numpy as np


COMPONENTS: Final[tuple[str, ...]] = (
    *('xx', 'xy', 'xz'),
    *('yx', 'yy', 'yz'),
    *('zx', 'zy', 'zz'),
)
VERSION: Final[str] = '0.1'


def make_matrix_svd_factors() -> np.ndarray:
  """Generates factors for converting 3x3 matrix to symmetric 4x4 matrix."""
  factors = np.zeros((16, 9), dtype=np.float32)

  factors[0, [0, 4, 8]] = 1.0

  factors[[1, 4], 5] = 1.0
  factors[[1, 4], 7] = -1.0

  factors[[2, 8], 6] = 1.0
  factors[[2, 8], 2] = -1.0

  factors[[3, 12], 1] = 1.0
  factors[[3, 12], 3] = -1.0

  factors[5, 0] = 1.0
  factors[5, [4, 8]] = -1.0

  factors[[6, 9], 1] = 1.0
  factors[[6, 9], 3] = 1.0

  factors[[7, 13], 2] = 1.0
  factors[[7, 13], 6] = 1.0

  factors[10, 4] = 1.0
  factors[10, [0, 8]] = -1.0

  factors[[11, 14], 5] = 1.0
  factors[[11, 14], 7] = 1.0

  factors[15, 8] = 1.0
  factors[15, [0, 4]] = -1.0

  return factors


@jax.custom_jvp
def largest_evec(m):
  _, eigvecs = jnp.linalg.eigh(m)
  return eigvecs[..., -1]


def largest_evec_jvp(primals, tangents):
  """jvp for largest eigenvector."""
  (m,) = primals
  (t,) = tangents
  eigvals, eigvecs = jnp.linalg.eigh(m)
  large_eigvec = eigvecs[..., -1]
  large_eigval = eigvals[..., -1]
  other_eigvals = eigvals[..., :-1]
  other_eigvecs = eigvecs[..., :-1]
  other_ev_times_tangent = jnp.einsum(
      '...aj,...ab -> ...bj',
      other_eigvecs,
      t,
      precision=jax.lax.Precision.HIGHEST,
  )
  nominator = jnp.einsum(
      '...bj,...b -> ...j',
      other_ev_times_tangent,
      large_eigvec,
      precision=jax.lax.Precision.HIGHEST,
  )
  prefactor = nominator / jnp.maximum(
      large_eigval[..., None] - other_eigvals, 1e-6
  )
  grad = jnp.sum(prefactor[..., None, :] * other_eigvecs, axis=-1)
  return large_eigvec, grad


largest_evec.defjvp(largest_evec_jvp)


MATRIX_SVD_QUAT_FACTORS = make_matrix_svd_factors()


@struct_of_array.StructOfArray(same_dtype=True)
class Rot3Array:
  """Rot3Array Matrix in 3 dimensional Space implemented as struct of arrays."""

  xx: jnp.ndarray = dataclasses.field(metadata={'dtype': jnp.float32})
  xy: jnp.ndarray
  xz: jnp.ndarray
  yx: jnp.ndarray
  yy: jnp.ndarray
  yz: jnp.ndarray
  zx: jnp.ndarray
  zy: jnp.ndarray
  zz: jnp.ndarray

  __array_ufunc__ = None

  def inverse(self) -> Self:
    """Returns inverse of Rot3Array."""
    return Rot3Array(
        *(self.xx, self.yx, self.zx),
        *(self.xy, self.yy, self.zy),
        *(self.xz, self.yz, self.zz),
    )

  def apply_to_point(self, point: vector.Vec3Array) -> vector.Vec3Array:
    """Applies Rot3Array to point."""
    return vector.Vec3Array(
        self.xx * point.x + self.xy * point.y + self.xz * point.z,
        self.yx * point.x + self.yy * point.y + self.yz * point.z,
        self.zx * point.x + self.zy * point.y + self.zz * point.z,
    )

  def apply_inverse_to_point(self, point: vector.Vec3Array) -> vector.Vec3Array:
    """Applies inverse Rot3Array to point."""
    return self.inverse().apply_to_point(point)

  def __matmul__(self, other: Self) -> Self:
    """Composes two Rot3Arrays."""
    c0 = self.apply_to_point(vector.Vec3Array(other.xx, other.yx, other.zx))
    c1 = self.apply_to_point(vector.Vec3Array(other.xy, other.yy, other.zy))
    c2 = self.apply_to_point(vector.Vec3Array(other.xz, other.yz, other.zz))
    return Rot3Array(c0.x, c1.x, c2.x, c0.y, c1.y, c2.y, c0.z, c1.z, c2.z)

  @classmethod
  def identity(cls, shape: Any, dtype: jnp.dtype = jnp.float32) -> Self:
    """Returns identity of given shape."""
    ones = jnp.ones(shape, dtype=dtype)
    zeros = jnp.zeros(shape, dtype=dtype)
    return cls(ones, zeros, zeros, zeros, ones, zeros, zeros, zeros, ones)  # pytype: disable=wrong-arg-count  # trace-all-classes

  @classmethod
  def from_two_vectors(cls, e0: vector.Vec3Array, e1: vector.Vec3Array) -> Self:
    """Construct Rot3Array from two Vectors.

    Rot3Array is constructed such that in the corresponding frame 'e0' lies on
    the positive x-Axis and 'e1' lies in the xy plane with positive sign of y.

    Args:
      e0: Vector
      e1: Vector

    Returns:
      Rot3Array
    """
    # Normalize the unit vector for the x-axis, e0.
    e0 = e0.normalized()
    # make e1 perpendicular to e0.
    c = e1.dot(e0)
    e1 = (e1 - c * e0).normalized()
    # Compute e2 as cross product of e0 and e1.
    e2 = e0.cross(e1)
    return cls(e0.x, e1.x, e2.x, e0.y, e1.y, e2.y, e0.z, e1.z, e2.z)  # pytype: disable=wrong-arg-count  # trace-all-classes

  @classmethod
  def from_array(cls, array: jnp.ndarray) -> Self:
    """Construct Rot3Array Matrix from array of shape. [..., 3, 3]."""
    unstacked = utils.unstack(array, axis=-2)
    unstacked = sum([utils.unstack(x, axis=-1) for x in unstacked], [])
    return cls(*unstacked)

  def to_array(self) -> jnp.ndarray:
    """Convert Rot3Array to array of shape [..., 3, 3]."""
    return jnp.stack(
        [
            jnp.stack([self.xx, self.xy, self.xz], axis=-1),
            jnp.stack([self.yx, self.yy, self.yz], axis=-1),
            jnp.stack([self.zx, self.zy, self.zz], axis=-1),
        ],
        axis=-2,
    )

  @classmethod
  def from_quaternion(
      cls,
      w: jnp.ndarray,
      x: jnp.ndarray,
      y: jnp.ndarray,
      z: jnp.ndarray,
      normalize: bool = True,
      epsilon: float = 1e-6,
  ) -> Self:
    """Construct Rot3Array from components of quaternion."""
    if normalize:
      inv_norm = jax.lax.rsqrt(jnp.maximum(epsilon, w**2 + x**2 + y**2 + z**2))
      w *= inv_norm
      x *= inv_norm
      y *= inv_norm
      z *= inv_norm
    xx = 1 - 2 * (jnp.square(y) + jnp.square(z))
    xy = 2 * (x * y - w * z)
    xz = 2 * (x * z + w * y)
    yx = 2 * (x * y + w * z)
    yy = 1 - 2 * (jnp.square(x) + jnp.square(z))
    yz = 2 * (y * z - w * x)
    zx = 2 * (x * z - w * y)
    zy = 2 * (y * z + w * x)
    zz = 1 - 2 * (jnp.square(x) + jnp.square(y))
    return cls(xx, xy, xz, yx, yy, yz, zx, zy, zz)  # pytype: disable=wrong-arg-count  # trace-all-classes

  @classmethod
  def from_svd(cls, mat: jnp.ndarray, use_quat_formula: bool = True) -> Self:
    """Constructs Rot3Array from arbitrary array of shape [3 * 3] using SVD.

    The case when 'use_quat_formula' is False rephrases the problem of
    projecting the matrix to a rotation matrix as a problem of finding the
    largest eigenvector of a certain 4x4 matrix. This has the advantage of
    having fewer numerical issues.
    This approach follows:
    https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.65.971&rep=rep1&type=pdf
    In the other case we construct it via svd following
    https://arxiv.org/pdf/2006.14616.pdf
    In that case [∂L/∂M] is large if the two smallest singular values are close
    to each other, or if they are close to 0.

    Args:
      mat: Array of shape [..., 3 * 3]
      use_quat_formula: Whether to construct matrix via 4x4 eigenvalue problem.

    Returns:
      Rot3Array of shape [...]
    """
    assert mat.shape[-1] == 9
    if use_quat_formula:
      symmetric_4by4 = jnp.einsum(
          'ji, ...i -> ...j',
          MATRIX_SVD_QUAT_FACTORS,
          mat,
          precision=jax.lax.Precision.HIGHEST,
      )
      symmetric_4by4 = jnp.reshape(symmetric_4by4, mat.shape[:-1] + (4, 4))
      largest_eigvec = largest_evec(symmetric_4by4)
      return cls.from_quaternion(
          *utils.unstack(largest_eigvec, axis=-1)
      ).inverse()

    else:
      mat = jnp.reshape(mat, mat.shape[:-1] + (3, 3))
      u, _, v_t = jnp.linalg.svd(mat, full_matrices=False)
      det_uv_t = jnp.linalg.det(
          jnp.matmul(u, v_t, precision=jax.lax.Precision.HIGHEST)
      )
      ones = jnp.ones_like(det_uv_t)
      diag_array = jnp.stack([ones, ones, det_uv_t], axis=-1)
      # This is equivalent to making diag_array into a diagonal array and matrix
      # multiplying
      diag_times_v_t = diag_array[..., None] * v_t
      out = jnp.matmul(u, diag_times_v_t, precision=jax.lax.Precision.HIGHEST)
      return cls.from_array(out)

  @classmethod
  def random_uniform(cls, key, shape, dtype=jnp.float32) -> Self:
    """Samples uniform random Rot3Array according to Haar Measure."""
    quat_array = jax.random.normal(key, tuple(shape) + (4,), dtype=dtype)
    quats = utils.unstack(quat_array)
    return cls.from_quaternion(*quats)

  def __getstate__(self):
    return (VERSION, [np.asarray(getattr(self, field)) for field in COMPONENTS])

  def __setstate__(self, state):
    version, state = state
    del version
    for i, field in enumerate(COMPONENTS):
      object.__setattr__(self, field, state[i])
