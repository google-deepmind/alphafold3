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
      dict(atom_key=[], exp_bonds_key=[]),
      dict(atom_key=[0, 1, 2, 3, 4, 5], exp_bonds_key=[2, 3, 4, 5]),
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

  def test_restrict_to_atoms_int32_query_with_large_bond_key(self):
    bonds = bonds_table.Bonds(
        key=np.array([2, 3], dtype=np.int64),
        type=np.array(['.', '.'], dtype=object),
        role=np.array(['.', '.'], dtype=object),
        from_atom_key=np.array([4, 4], dtype=np.int64),
        dest_atom_key=np.array([5, 10**12], dtype=np.int64),
    )
    restr_bonds = bonds.restrict_to_atoms(np.array([4, 5], dtype=np.int32))
    np.testing.assert_array_equal(
        restr_bonds.key, np.array([2], dtype=np.int64), strict=True
    )

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

  @parameterized.named_parameters(
      dict(
          testcase_name='non_empty_atoms',
          atom_key=[0, 1, 2],
          atom_ends=[0, 2, 3],
          expected_len=3,
      ),
      dict(
          testcase_name='empty_atoms_and_ranges',
          atom_key=[],
          atom_ends=[],
          expected_len=0,
      ),
  )
  def test_split_by_atom_ranges_empty_bonds_returns_self(
      self, atom_key, atom_ends, expected_len
  ):
    empty = bonds_table.Bonds.make_empty()
    split = empty.split_by_atom_ranges(
        atom_key=np.array(atom_key, dtype=np.int64),
        atom_ends=np.array(atom_ends, dtype=np.int64),
    )
    self.assertLen(split, expected_len)
    for part in split:
      self.assertIs(part, empty)

  def test_split_by_atom_ranges_full_cover_returns_self(self):
    bonds = bonds_table.Bonds(
        key=np.array([10, 11], dtype=np.int64),
        type=np.array(['covale', 'disulf'], dtype=object),
        role=np.array(['r0', 'r1'], dtype=object),
        from_atom_key=np.array([100, 101], dtype=np.int64),
        dest_atom_key=np.array([101, 102], dtype=np.int64),
    )
    full = bonds.split_by_atom_ranges(
        atom_key=np.array([100, 101, 102], dtype=np.int64),
        atom_ends=np.array([0, 3], dtype=np.int64),
    )
    self.assertLen(full, 2)
    self.assertEmpty(full[0])
    self.assertIs(full[1], bonds)

  @parameterized.named_parameters(
      dict(
          testcase_name='drops_inter_range_bond',
          atom_key=[100, 101, 102, 103],
          from_atom_key=[100, 101, 102],
          dest_atom_key=[101, 102, 103],
          atom_ends=np.array([2, 4]),
          expected_keys=[[10], [12]],
      ),
      dict(
          testcase_name='multi_bond_contiguous_slice',
          atom_key=[100, 101, 102, 103],
          from_atom_key=[100, 100, 102],
          dest_atom_key=[101, 101, 103],
          atom_ends=np.array([2, 4]),
          expected_keys=[[10, 11], [12]],
      ),
      dict(
          testcase_name='unsorted_atom_keys_and_self_loops',
          atom_key=[103, 101, 100, 102],
          from_atom_key=[103, 100, 102],
          dest_atom_key=[101, 100, 102],
          atom_ends=np.array([2, 3]),
          expected_keys=[[10], [11]],
      ),
      dict(
          testcase_name='reordered_runs_with_empty_range',
          atom_key=[100, 101, 102, 103],
          from_atom_key=[102, 100, 101],
          dest_atom_key=[103, 101, 102],
          atom_ends=np.array([2, 2, 4]),
          expected_keys=[[11], [], [10]],
      ),
      dict(
          testcase_name='interleaved_runs_preserve_order_int32',
          atom_key=[100, 101, 102, 103],
          from_atom_key=[102, 100, 102],
          dest_atom_key=[103, 101, 103],
          atom_ends=np.array([2, 4], dtype=np.int32),
          expected_keys=[[11], [10, 12]],
      ),
      dict(
          testcase_name='uint64_atom_ends',
          atom_key=[100, 101, 102, 103],
          from_atom_key=[102, 100, 102],
          dest_atom_key=[103, 101, 103],
          atom_ends=np.array([2, 4], dtype=np.uint64),
          expected_keys=[[11], [10, 12]],
      ),
      dict(
          testcase_name='prefix_ranges_drop_trailing_bonds',
          atom_key=[100, 101, 102, 103],
          from_atom_key=[102, 100, 102],
          dest_atom_key=[103, 101, 103],
          atom_ends=np.array([2], dtype=np.int32),
          expected_keys=[[11]],
      ),
      dict(
          testcase_name='all_inter_range_bonds',
          atom_key=[100, 101, 102, 103],
          from_atom_key=[100, 101, 100],
          dest_atom_key=[102, 103, 103],
          atom_ends=np.array([2, 4]),
          expected_keys=[[], []],
      ),
      dict(
          testcase_name='empty_atom_ends',
          atom_key=[100, 101, 102, 103],
          from_atom_key=[100, 101, 102],
          dest_atom_key=[101, 102, 103],
          atom_ends=np.empty((0,), dtype=np.int64),
          expected_keys=[],
      ),
  )
  def test_split_by_atom_ranges(
      self,
      atom_key,
      from_atom_key,
      dest_atom_key,
      atom_ends,
      expected_keys,
  ):
    from_arr = np.array(from_atom_key, dtype=np.int64)
    dest_arr = np.array(dest_atom_key, dtype=np.int64)
    bonds = bonds_table.Bonds(
        key=np.array([10, 11, 12], dtype=np.int64),
        type=np.array(['a', 'b', 'c'], dtype=object),
        role=np.array(['r0', 'r1', 'r2'], dtype=object),
        from_atom_key=from_arr,
        dest_atom_key=dest_arr,
    )
    split = bonds.split_by_atom_ranges(
        atom_key=np.array(atom_key, dtype=np.int64),
        atom_ends=atom_ends,
    )
    self.assertLen(split, len(expected_keys))
    for part, exp_keys in zip(split, expected_keys, strict=True):
      rows = [k - 10 for k in exp_keys]
      np.testing.assert_array_equal(
          part.key, np.array(exp_keys, dtype=np.int64), strict=True
      )
      np.testing.assert_array_equal(
          part.type,
          np.array([['a', 'b', 'c'][r] for r in rows], dtype=object),
          strict=True,
      )
      np.testing.assert_array_equal(
          part.role,
          np.array([['r0', 'r1', 'r2'][r] for r in rows], dtype=object),
          strict=True,
      )
      np.testing.assert_array_equal(
          part.from_atom_key, from_arr[rows], strict=True
      )
      np.testing.assert_array_equal(
          part.dest_atom_key, dest_arr[rows], strict=True
      )
      self.assertFalse(part.key.flags.writeable)

  @parameterized.named_parameters(
      dict(
          testcase_name='missing_bond_atom_key',
          atom_key=np.array([0, 1], dtype=np.int64),
          atom_ends=np.array([2]),
          empty_bonds=False,
          expected_error='Bond atom key missing from atom_key',
      ),
      dict(
          testcase_name='empty_atom_key_with_nonempty_bonds',
          atom_key=np.empty((0,), dtype=np.int64),
          atom_ends=np.empty((0,), dtype=np.int64),
          empty_bonds=False,
          expected_error='Bond atom key missing from atom_key',
      ),
      dict(
          testcase_name='two_dim_atom_key_empty_bonds',
          atom_key=np.zeros((2, 2), dtype=np.int64),
          atom_ends=np.array([2, 4]),
          empty_bonds=True,
          expected_error='atom_key must be a 1D signed-integer array',
      ),
      dict(
          testcase_name='uint64_atom_key',
          atom_key=np.array([0, 1, 2, 3], dtype=np.uint64),
          atom_ends=np.array([2, 4]),
          empty_bonds=False,
          expected_error='atom_key must be a 1D signed-integer array',
      ),
      dict(
          testcase_name='decreasing_atom_ends',
          atom_key=np.array([0, 1, 2, 3], dtype=np.int64),
          atom_ends=np.array([3, 2]),
          empty_bonds=False,
          expected_error='atom_ends must be a 1D non-decreasing integer array',
      ),
      dict(
          testcase_name='atom_end_past_num_atoms_empty_bonds',
          atom_key=np.array([0, 1, 2, 3], dtype=np.int64),
          atom_ends=np.array([2, 5]),
          empty_bonds=True,
          expected_error='atom_ends must be a 1D non-decreasing integer array',
      ),
      dict(
          testcase_name='negative_atom_end',
          atom_key=np.array([0, 1, 2, 3], dtype=np.int64),
          atom_ends=np.array([-1, 4]),
          empty_bonds=False,
          expected_error='atom_ends must be a 1D non-decreasing integer array',
      ),
      dict(
          testcase_name='float_atom_ends',
          atom_key=np.array([0, 1, 2, 3], dtype=np.int64),
          atom_ends=np.array([2.0, 4.0]),
          empty_bonds=False,
          expected_error='atom_ends must be a 1D non-decreasing integer array',
      ),
      dict(
          testcase_name='bool_atom_ends',
          atom_key=np.array([0, 1, 2, 3], dtype=np.int64),
          atom_ends=np.array([False, True]),
          empty_bonds=False,
          expected_error='atom_ends must be a 1D non-decreasing integer array',
      ),
      dict(
          testcase_name='zero_dim_atom_ends',
          atom_key=np.array([0, 1, 2, 3], dtype=np.int64),
          atom_ends=np.array(2),
          empty_bonds=False,
          expected_error='atom_ends must be a 1D non-decreasing integer array',
      ),
      dict(
          testcase_name='two_dim_atom_ends',
          atom_key=np.array([0, 1, 2, 3], dtype=np.int64),
          atom_ends=np.array([[2, 4]]),
          empty_bonds=False,
          expected_error='atom_ends must be a 1D non-decreasing integer array',
      ),
  )
  def test_split_by_atom_ranges_invalid_inputs_raise(
      self, atom_key, atom_ends, empty_bonds, expected_error
  ):
    if empty_bonds:
      bonds = bonds_table.Bonds.make_empty()
    else:
      bonds = bonds_table.Bonds(
          key=np.array([1], dtype=np.int64),
          type=np.array(['covale'], dtype=object),
          role=np.array(['?'], dtype=object),
          from_atom_key=np.array([0], dtype=np.int64),
          dest_atom_key=np.array([3], dtype=np.int64),
      )
    with self.assertRaisesRegex(ValueError, expected_error):
      bonds.split_by_atom_ranges(atom_key=atom_key, atom_ends=atom_ends)


if __name__ == '__main__':
  absltest.main()
