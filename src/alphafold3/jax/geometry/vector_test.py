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

import pickle

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.jax.geometry import test_utils
from alphafold3.jax.geometry import vector
import jax
import jax.numpy as jnp
import numpy as np

jax.config.update('jax_threefry_partitionable', False)


def make_vector_array(shape, dtype=jnp.float32, random=False, key=42):
  shape = tuple(shape) + (3,)
  if random:
    key = jax.random.PRNGKey(key)
    return jax.random.normal(key, shape, dtype)
  else:
    return jnp.ones(shape, dtype=dtype)


class Vec3ArrayTest(parameterized.TestCase):

  @parameterized.parameters([{'shape': (12,)}])
  def test_from_array(self, shape=(12,)):
    vec_array = make_vector_array(shape)
    vec = vector.Vec3Array.from_array(vec_array)
    test_utils.assert_array_equal_to_vector(vec_array, vec)

  @parameterized.parameters([
      {'shape': (12,), 'key': 0},
      {'shape': (12,), 'key': slice(None, 3, None)},
      {'shape': (12,), 'key': ...},
  ])
  def test_indexing(self, shape=(12,), key=slice(None, 3, None)):
    vec_array = make_vector_array(shape)
    vec = vector.Vec3Array.from_array(vec_array)
    vec_array = vec_array[key]
    vec = vec[key]
    test_utils.assert_array_equal_to_vector(vec_array, vec)

  @parameterized.parameters([
      {'left_shape': (12,), 'right_shape': (12,)},
      {'left_shape': (12,), 'right_shape': (1,)},
      {'left_shape': (12,), 'right_shape': ()},
  ])
  def test_add(self, left_shape=(12,), right_shape=(12,)):
    vec1_array = make_vector_array(left_shape, random=True, key=16)
    vec1 = vector.Vec3Array.from_array(vec1_array)

    vec2_array = make_vector_array(right_shape, random=True, key=16)
    vec2 = vector.Vec3Array.from_array(vec2_array)

    out_vec = vec1 + vec2
    out_array = vec1_array + vec2_array

    test_utils.assert_array_equal_to_vector(out_array, out_vec)

  @parameterized.parameters(
      [{'shape': (12,), 'reverse': False}, {'shape': (12,), 'reverse': True}]
  )
  def test_mul(self, shape=(12,), reverse=False):
    vec_array = make_vector_array(shape, random=True, key=16)
    vec = vector.Vec3Array.from_array(vec_array)

    factor = jax.random.normal(jax.random.PRNGKey(42), [], dtype=jnp.float32)

    if reverse:
      out_vec = factor * vec
      out_array = factor * vec_array
    else:
      out_vec = vec * factor
      out_array = vec_array * factor

    test_utils.assert_array_equal_to_vector(out_array, out_vec)

  @parameterized.parameters([
      {'left_shape': (12,), 'right_shape': (12,)},
      {'left_shape': (12,), 'right_shape': (1,)},
      {'left_shape': (12,), 'right_shape': ()},
  ])
  def test_sub(self, left_shape=(12,), right_shape=(12,)):
    vec1_array = make_vector_array(left_shape, random=True, key=16)
    vec1 = vector.Vec3Array.from_array(vec1_array)

    vec2_array = make_vector_array(right_shape, random=True, key=16)
    vec2 = vector.Vec3Array.from_array(vec2_array)

    out_vec = vec1 - vec2
    out_vec_true = vec1 + (-1 * vec2)

    test_utils.assert_vectors_equal(out_vec_true, out_vec)

  @parameterized.parameters([{'shape': (12,)}])
  def test_div(self, shape=(12,)):
    vec_array = make_vector_array(shape, random=True, key=16)
    vec = vector.Vec3Array.from_array(vec_array)

    factor = jax.random.normal(jax.random.PRNGKey(42), [], dtype=jnp.float32)

    out_vec = vec / factor
    out_vec_true = vec * (1.0 / factor)

    test_utils.assert_vectors_close(out_vec_true, out_vec)

  @parameterized.parameters([{'shape': (12,)}])
  def test_neg(self, shape=(12,)):
    vec_array = make_vector_array(shape, random=True, key=16)
    vec = vector.Vec3Array.from_array(vec_array)

    out_vec = -vec
    out_vec_true = (-1.0) * vec

    test_utils.assert_vectors_equal(out_vec_true, out_vec)

  @parameterized.parameters([
      {'left_shape': (12,), 'right_shape': (12,)},
      {'left_shape': (12,), 'right_shape': (1,)},
      {'left_shape': (12,), 'right_shape': ()},
  ])
  def test_cross(self, left_shape=(12,), right_shape=(12,)):
    vec1_array = make_vector_array(left_shape, random=True, key=16)
    vec1 = vector.Vec3Array.from_array(vec1_array)

    vec2_array = make_vector_array(right_shape, random=True, key=16)
    vec2 = vector.Vec3Array.from_array(vec2_array)

    out_vec = vec1.cross(vec2)
    out_array = jnp.cross(vec1_array, vec2_array)

    test_utils.assert_array_close_to_vector(out_array, out_vec)

  @parameterized.parameters([
      {'left_shape': (12,), 'right_shape': (12,)},
      {'left_shape': (12,), 'right_shape': (1,)},
      {'left_shape': (12,), 'right_shape': ()},
  ])
  def test_dot(self, left_shape=(12,), right_shape=(12,)):
    vec1_array = make_vector_array(left_shape, random=True, key=16)
    vec1 = vector.Vec3Array.from_array(vec1_array)

    vec2_array = make_vector_array(right_shape, random=True, key=16)
    vec2 = vector.Vec3Array.from_array(vec2_array)

    out_vec = vec1.dot(vec2)
    out_array = jnp.sum(vec1_array * vec2_array, axis=-1)

    np.testing.assert_allclose(out_array, out_vec, rtol=1e-06)

  @parameterized.parameters([{'shape': (12,)}])
  def test_norm(self, shape=(12,)):
    vec_array = make_vector_array(shape, random=True, key=16)
    vec = vector.Vec3Array.from_array(vec_array)

    vec_norm = vec.norm()
    array_norm = jnp.linalg.norm(vec_array, axis=-1)
    np.testing.assert_allclose(array_norm, vec_norm, rtol=1e-6)

  @parameterized.parameters([{'shape': (12,)}])
  def test_norm2(self, shape=(12,)):
    vec_array = make_vector_array(shape, random=True, key=16)
    vec = vector.Vec3Array.from_array(vec_array)
    np.testing.assert_allclose(
        vec.norm2(), jnp.square(vec.norm(epsilon=1e-30)), rtol=1e-6
    )

  @parameterized.parameters([{'shape': (12,)}])
  def test_normalized(self, shape=(12,)):
    vec_array = make_vector_array(shape, random=True, key=16)
    vec = vector.Vec3Array.from_array(vec_array)
    vec_normalized = vec.normalized()

    # Normalized vector should have norm 1
    np.testing.assert_allclose(vec_normalized.norm(), 1.0, atol=1e-6)

    # Normalizing shouldn't change direction
    dot_product = vec.dot(vec_normalized)
    one = dot_product / vec.norm()
    np.testing.assert_allclose(one, 1.0, atol=1e-6)

  @parameterized.parameters([{'shape': (12,), 'dtype': jnp.float32}])
  def test_dtype(self, shape=(12,), dtype=jnp.float32):
    vec_array = make_vector_array(shape, dtype=dtype, random=True, key=16)
    vec = vector.Vec3Array.from_array(vec_array)
    self.assertEqual(vec.dtype, dtype)

  @parameterized.parameters([{'shape': (12,)}, {'shape': ()}])
  def test_shape(self, shape=(12,), dtype=jnp.float32):
    vec_array = make_vector_array(shape, dtype=dtype, random=True, key=16)
    vec = vector.Vec3Array.from_array(vec_array)
    self.assertEqual(vec.shape, shape)

  @parameterized.parameters([{'shape': (12,)}, {'shape': ()}])
  def test_len(self, shape=(12,), dtype=jnp.float32):
    vec_array = make_vector_array(shape, dtype=dtype, random=True, key=16)
    vec = vector.Vec3Array.from_array(vec_array)
    if shape:
      self.assertEqual(len(vec), len(np.zeros(shape)))
    else:
      with self.assertRaises(TypeError):
        len(vec)

  @parameterized.parameters([{'shape': (12,)}, {'shape': (12, 5)}])
  def test_zeros(self, shape=(12,), dtype=jnp.float32):
    zeros_array = jnp.zeros(shape + (3,), dtype)
    vec_zeros = vector.Vec3Array.zeros(shape, dtype)
    test_utils.assert_array_equal_to_vector(zeros_array, vec_zeros)

  @parameterized.parameters([{'shape': (12,)}, {'shape': (12, 5)}])
  def test_to_array(self, shape=(12,), dtype=jnp.float32):
    vec_array = make_vector_array(shape, dtype=dtype, random=True, key=16)
    vec = vector.Vec3Array.from_array(vec_array)
    vec_array_from_vec = vec.to_array()
    self.assertTrue(jnp.all(jnp.equal(vec_array_from_vec, vec_array)))

  @parameterized.parameters([{'shape': (12,)}, {'shape': (12, 5)}])
  def test_pickle(self, shape=(12,), dtype=jnp.float32):
    vec_array = make_vector_array(shape, dtype=dtype, random=True, key=16)
    vec = vector.Vec3Array.from_array(vec_array)
    for protocol in range(1, pickle.HIGHEST_PROTOCOL + 1):
      vec_after = pickle.loads(pickle.dumps(vec, protocol=protocol))
      test_utils.assert_vectors_equal(vec, vec_after)


if __name__ == '__main__':
  absltest.main()
