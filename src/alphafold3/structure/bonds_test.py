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

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.structure import bonds as bonds_table
import numpy as np


class BondsTest(parameterized.TestCase):

  @parameterized.parameters([
      dict(
          atom_key=[0, 1, 2, 3, 4, 5],
          exp_from_idxs=[4, 4, 0, 5],
          exp_dest_idxs=[5, 3, 1, 2],
      ),
      dict(
          atom_key=[5, 4, 3, 2, 1, 0],
          exp_from_idxs=[1, 1, 5, 0],
          exp_dest_idxs=[0, 2, 4, 3],
      ),
  ])
  def test_get_atom_indices(self, atom_key, exp_from_idxs, exp_dest_idxs):
    bonds = bonds_table.Bonds(
        key=np.array([2, 3, 4, 5]),
        type=np.array(['.', '.', '.', '.'], dtype=object),
        role=np.array(['.', '.', '.', '.'], dtype=object),
        from_atom_key=np.array([4, 4, 0, 5]),
        dest_atom_key=np.array([5, 3, 1, 2]),
    )
    from_atom_idxs, dest_atom_idxs = bonds.get_atom_indices(np.array(atom_key))
    np.testing.assert_array_equal(from_atom_idxs, exp_from_idxs)
    np.testing.assert_array_equal(dest_atom_idxs, exp_dest_idxs)

  @parameterized.parameters([
      dict(atom_key=[0, 1, 2, 3, 4]),
      dict(atom_key=[4, 3, 2, 1, 5]),
  ])
  def test_get_atom_indices_missing(self, atom_key):
    bonds = bonds_table.Bonds(
        key=np.array([1, 2, 3, 4]),
        type=np.array(['.', '.', '.', '.'], dtype=object),
        role=np.array(['.', '.', '.', '.'], dtype=object),
        from_atom_key=np.array([3, 3, 5, 4]),
        dest_atom_key=np.array([4, 2, 0, 1]),
    )
    with self.assertRaises(ValueError):
      bonds.get_atom_indices(np.array(atom_key))

  @parameterized.parameters([
      dict(atom_key=[4, 5], exp_bonds_key=[2]),
      dict(atom_key=[4, 2, 5], exp_bonds_key=[2, 5]),
      dict(atom_key=[4], exp_bonds_key=[]),
  ])
  def test_restrict_to_atoms(self, atom_key, exp_bonds_key):
    bonds = bonds_table.Bonds(
        key=np.array([2, 3, 4, 5]),
        type=np.array(['.', '.', '.', '.'], dtype=object),
        role=np.array(['.', '.', '.', '.'], dtype=object),
        from_atom_key=np.array([4, 4, 0, 5]),
        dest_atom_key=np.array([5, 3, 1, 2]),
    )
    restr_bonds = bonds.restrict_to_atoms(np.array(atom_key))
    np.testing.assert_array_equal(restr_bonds.key, np.array(exp_bonds_key))

  def test_concat_with_atom_keys_no_bonds(self):
    bonds_tables = [None, None, None]
    atom_key_arrays = [[1], [1, 2], [2]]
    atom_key_arrays = [np.array(a) for a in atom_key_arrays]
    bonds, atom_key = bonds_table.concat_with_atom_keys(
        bonds_tables, atom_key_arrays
    )
    self.assertIsNone(bonds)

    # We don't mind how atom keys are assigned, as long as they are unique.
    unique_atom_key = np.unique(atom_key)
    self.assertLen(atom_key, len(unique_atom_key))

  def test_concat_with_atom_keys(self):
    bonds_tables = [
        bonds_table.Bonds(
            key=np.array([1, 2]),
            type=np.array(['a', 'b'], dtype=object),
            role=np.array(['a', 'b'], dtype=object),
            from_atom_key=np.array([4, 0]),
            dest_atom_key=np.array([0, 3]),
        ),
        None,
        bonds_table.Bonds(
            key=np.array([1, 2, 3]),
            type=np.array(['a', 'b', 'c'], dtype=object),
            role=np.array(['a', 'b', 'c'], dtype=object),
            from_atom_key=np.array([1, 2, 4]),
            dest_atom_key=np.array([2, 4, 1]),
        ),
    ]
    atom_key_arrays = [[0, 4, 3], [0, 3], [1, 2, 4]]
    atom_key_arrays = [np.array(a) for a in atom_key_arrays]
    bonds, atom_key = bonds_table.concat_with_atom_keys(
        bonds_tables, atom_key_arrays
    )
    self.assertIsNotNone(bonds)

    # We don't mind how bond keys are assigned, as long as they are unique.
    unique_bond_key = np.unique(bonds.key)
    self.assertLen(bonds.key, len(unique_bond_key))
    np.testing.assert_array_equal(
        bonds.type, np.array(['a', 'b', 'a', 'b', 'c'])
    )
    np.testing.assert_array_equal(
        bonds.role, np.array(['a', 'b', 'a', 'b', 'c'])
    )

    # In order to be agnostic to how keys are assigned, we instead check
    # the indices of the bonded atoms.
    from_atom_idxs, dest_atom_idxs = bonds.get_atom_indices(atom_key)
    np.testing.assert_array_equal(from_atom_idxs, [1, 0, 5, 6, 7])
    np.testing.assert_array_equal(dest_atom_idxs, [0, 2, 6, 7, 5])

    # We don't mind how atom keys are assigned, as long as they are unique.
    unique_atom_key = np.unique(atom_key)
    self.assertLen(atom_key, len(unique_atom_key))

  def test_concat_with_atom_keys_all_empty(self):
    bonds_tables = [None, None, None]
    atom_key_arrays = [[], [], []]
    atom_key_arrays = [np.array(a, dtype=np.int64) for a in atom_key_arrays]
    bonds, atom_key = bonds_table.concat_with_atom_keys(
        bonds_tables, atom_key_arrays
    )
    self.assertIsNone(bonds)
    self.assertEqual(0, atom_key.size)

  def test_empty_falsiness(self):
    empty_bonds = bonds_table.Bonds.make_empty()
    self.assertEmpty(empty_bonds)
    self.assertFalse(empty_bonds)  # Implicitly False, since len is 0.


if __name__ == '__main__':
  absltest.main()
