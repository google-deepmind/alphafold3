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

import itertools
import pickle

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.jax.geometry import rigid_matrix_vector
from alphafold3.jax.geometry import rotation_matrix
from alphafold3.jax.geometry import test_utils
from alphafold3.jax.geometry import vector
import jax
import jax.numpy as jnp

jax.config.update('jax_threefry_partitionable', False)


def make_rigid(shape, dtype=jnp.float32, random=False, key=42):
  if random:
    key = jax.random.PRNGKey(key)
    vec_key, rot_key = jax.random.split(key, 2)
    vec = vector.random_gaussian_vector(shape, vec_key, dtype=dtype)
    rot = rotation_matrix.Rot3Array.random_uniform(rot_key, shape, dtype=dtype)
    return rigid_matrix_vector.Rigid3Array(rot, vec), rot, vec  # pyrefly: ignore[bad-argument-count]
  else:
    vec = vector.Vec3Array.zeros(shape, dtype)
    rot = rotation_matrix.Rot3Array.identity(shape, dtype)  # pyrefly: ignore[bad-argument-type]
    return rigid_matrix_vector.Rigid3Array(rot, vec), rot, vec  # pyrefly: ignore[bad-argument-count]


class RigidMatrixVectorTest(parameterized.TestCase):

  @parameterized.parameters([
      {'shape': (12,), 'key': 0},
      {'shape': (12,), 'key': slice(None, 3, None)},
      {'shape': (12,), 'key': ...},
  ])
  def test_indexing(self, shape=(12,), key=slice(None, 3, None)):
    rigid, rot, vec = make_rigid(shape, random=True)
    rigid = rigid[key]
    rot = rot[key]
    vec = vec[key]
    test_utils.assert_rot_trans_equal_to_rigid(rot, vec, rigid)

  @parameterized.parameters([{'shape': (12,)}, {'shape': ()}])
  def test_identity(self, shape=(12,), dtype=jnp.float32):
    id_rigid = rigid_matrix_vector.Rigid3Array.identity(shape, dtype)  # pyrefly: ignore[bad-argument-type]
    id_rigid_direct, _, _ = make_rigid(shape, random=False, dtype=dtype)
    test_utils.assert_rigid_equal_to_rigid(id_rigid_direct, id_rigid)

  @parameterized.parameters([
      {'left_shape': (12,), 'right_shape': (12,)},
      {'left_shape': (12,), 'right_shape': (1,)},
      {'left_shape': (12,), 'right_shape': ()},
  ])
  def test_compose(self, left_shape=(12,), right_shape=(12,)):
    left_rigid, left_rot, left_vec = make_rigid(left_shape, random=True, key=18)
    right_rigid, right_rot, right_vec = make_rigid(
        right_shape, random=True, key=22
    )
    out_rigid = left_rigid @ right_rigid
    out_rots = left_rot @ right_rot
    out_vec = left_vec + left_rot.apply_to_point(right_vec)
    test_utils.assert_rot_trans_close_to_rigid(out_rots, out_vec, out_rigid)

  @parameterized.parameters([
      {'a_shape': (12,), 'b_shape': (12,), 'c_shape': (12,)},
      {'a_shape': (12,), 'b_shape': (1,), 'c_shape': (12,)},
      {'a_shape': (12,), 'b_shape': (), 'c_shape': ()},
  ])
  def test_compose_associative(
      self, a_shape=(12,), b_shape=(12,), c_shape=(12,)
  ):
    a, _, _ = make_rigid(a_shape, random=True, key=18)
    b, _, _ = make_rigid(b_shape, random=True, key=22)
    c, _, _ = make_rigid(c_shape, random=True, key=25)

    out_rigid1 = (a @ b) @ c
    out_rigid2 = a @ (b @ c)
    test_utils.assert_rigid_close_to_rigid(out_rigid1, out_rigid2)

  @parameterized.parameters([
      {'left_shape': (12,), 'right_shape': (12,), 'vec_shape': (12,)},
      {'left_shape': (12,), 'right_shape': (12,), 'vec_shape': ()},
      {'left_shape': (12,), 'right_shape': (), 'vec_shape': (12,)},
  ])
  def test_compose_apply(
      self, left_shape=(12,), right_shape=(12,), vec_shape=(12,)
  ):
    left_rigid, _, _ = make_rigid(left_shape, random=True, key=18)
    right_rigid, _, _ = make_rigid(right_shape, random=True, key=22)
    point = vector.random_gaussian_vector(
        shape=vec_shape, key=jax.random.PRNGKey(23)
    )
    out_point1 = (left_rigid @ right_rigid).apply_to_point(point)
    out_point2 = left_rigid.apply_to_point(right_rigid.apply_to_point(point))
    test_utils.assert_vectors_close(out_point1, out_point2)

  @parameterized.parameters([{'shape': (12,)}, {'shape': ()}])
  def test_inverse(self, shape=(12,), dtype=jnp.float32):
    id_rigid = rigid_matrix_vector.Rigid3Array.identity(shape, dtype)  # pyrefly: ignore[bad-argument-type]
    rigid, _, _ = make_rigid(shape, random=True, key=18, dtype=dtype)
    rigid_rigid_inv = rigid @ rigid.inverse()
    test_utils.assert_rigid_close_to_rigid(rigid_rigid_inv, id_rigid)
    rigid_inv_rigid = rigid.inverse() @ rigid
    test_utils.assert_rigid_close_to_rigid(rigid_inv_rigid, id_rigid)

  @parameterized.parameters([
      {'rigid_shape': (12,), 'vec_shape': (12,)},
      {'rigid_shape': (12,), 'vec_shape': (1,)},
      {'rigid_shape': (12,), 'vec_shape': ()},
  ])
  def test_apply_to_point(self, rigid_shape=(12,), vec_shape=(12,)):
    rigid, rot, vec = make_rigid(rigid_shape, random=True, key=18)
    point = vector.random_gaussian_vector(
        shape=vec_shape, key=jax.random.PRNGKey(23)
    )
    out_point_rigid = rigid.apply_to_point(point)
    out_point_rt = vec + rot.apply_to_point(point)
    test_utils.assert_vectors_close(out_point_rigid, out_point_rt)

  @parameterized.parameters([
      {'rigid_shape': (12,), 'vec_shape': (12,)},
      {'rigid_shape': (12,), 'vec_shape': (1,)},
      {'rigid_shape': (12,), 'vec_shape': ()},
  ])
  def test_apply_inverse_to_point(self, rigid_shape=(12,), vec_shape=(12,)):
    rigid, _, _ = make_rigid(rigid_shape, random=True, key=18)
    point = vector.random_gaussian_vector(
        shape=vec_shape, key=jax.random.PRNGKey(23)
    )
    out_point_rigid = rigid.apply_inverse_to_point(point)
    out_point_rigid_inv = rigid.inverse().apply_to_point(point)
    test_utils.assert_vectors_close(out_point_rigid, out_point_rigid_inv)

  @parameterized.parameters([
      {'shape': (12,), 'factor_shape': (12,)},
      {'shape': (12,), 'factor_shape': ()},
      {'shape': (), 'factor_shape': ()},
  ])
  def test_scale_translation(
      self, shape=(12,), factor_shape=(), dtype=jnp.float32
  ):
    factor = jax.random.normal(
        jax.random.PRNGKey(22), shape=factor_shape, dtype=dtype
    )
    rigid, rot, trans = make_rigid(shape, random=True, key=18, dtype=dtype)
    rigid = rigid.scale_translation(factor)
    test_utils.assert_vectors_close(rigid.translation, trans * factor)
    test_utils.assert_rotation_matrix_equal(rigid.rotation, rot)

  @parameterized.parameters([{'shape': (12,), 'dtype': jnp.float32}])
  def test_dtype(self, shape=(12,), dtype=jnp.float32):
    rigid, _, _ = make_rigid(shape, dtype=dtype, random=True, key=16)
    self.assertEqual(rigid.dtype, dtype)

  @parameterized.parameters([{'shape': (12,)}, {'shape': ()}])
  def test_shape(self, shape=(12,), dtype=jnp.float32):
    rigid, _, _ = make_rigid(shape, dtype=dtype, random=True, key=16)
    self.assertEqual(rigid.shape, shape)

  @parameterized.parameters([{'shape': (12,)}, {'shape': ()}])
  def test_flatten(self, shape=(12,), dtype=jnp.float32):
    rigid, _, _ = make_rigid(shape, dtype=dtype, random=True, key=16)
    rigid_flat, rigid_treedef = jax.tree_util.tree_flatten(rigid)
    unflattened_rigid = jax.tree_util.tree_unflatten(rigid_treedef, rigid_flat)
    test_utils.assert_rigid_equal_to_rigid(rigid, unflattened_rigid)

  @parameterized.parameters([{'shape': (12,)}])
  def test_pickle(self, shape=(12,), dtype=jnp.float32):
    rigid, _, _ = make_rigid(shape, dtype=dtype, random=True, key=16)
    for protocol in range(1, pickle.HIGHEST_PROTOCOL + 1):
      rigid_after = pickle.loads(pickle.dumps(rigid, protocol=protocol))
      test_utils.assert_rigid_equal_to_rigid(rigid, rigid_after)

  @parameterized.parameters(
      list(itertools.product([True, False], [True, False], [(), (8,)]))
  )
  def test_from_point_alignment(
      self, random_transformation, random_points, batch_shape=()
  ):

    if random_points:
      point_array = jax.random.normal(
          jax.random.PRNGKey(7), batch_shape + (3, 3)
      )
      points_from = vector.Vec3Array.from_array(point_array)
    else:
      points_from = vector.Vec3Array(
          x=jnp.array([1.0, 0.0, 0.0]),  # pyrefly: ignore[unexpected-keyword]
          y=jnp.array([0.0, 1.0, 0.0]),  # pyrefly: ignore[unexpected-keyword]
          z=jnp.array([0.0, 0.0, 1.0]),  # pyrefly: ignore[unexpected-keyword]
      )
      points_from = jax.tree.map(
          lambda x: jnp.broadcast_to(x, batch_shape + (3,)), points_from
      )

    if random_transformation:
      rot = rotation_matrix.Rot3Array.random_uniform(
          key=jax.random.PRNGKey(42), shape=batch_shape
      )
      translation_array = jax.random.normal(
          jax.random.PRNGKey(43), batch_shape + (3,)
      )
      translation = vector.Vec3Array.from_array(translation_array)
    else:
      rot = rotation_matrix.Rot3Array.from_array(
          jnp.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
      )
      rot = jax.tree.map(lambda x: jnp.broadcast_to(x, batch_shape), rot)
      translation = vector.Vec3Array.zeros([])
      translation = jax.tree.map(
          lambda x: jnp.broadcast_to(x, batch_shape), translation
      )
    rigid_to_apply = rigid_matrix_vector.Rigid3Array(rot, translation)  # pyrefly: ignore[bad-argument-count]

    points_to = rigid_to_apply[..., None].apply_to_point(points_from)  # pyrefly: ignore[bad-index]

    inferred_rigid = rigid_matrix_vector.Rigid3Array.from_point_alignment(
        points_to=points_to, points_from=points_from
    )

    test_utils.assert_rigid_close_to_rigid(
        rigid_to_apply, inferred_rigid, rotation_digits=4, translation_atol=2e-5
    )


if __name__ == '__main__':
  absltest.main()
