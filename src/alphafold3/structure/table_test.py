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

from collections.abc import Collection, Mapping
import dataclasses
from typing import Any, ClassVar, Final, Self

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.structure import table
import numpy as np


@dataclasses.dataclass(frozen=True)
class OneColTable(table.Table):
  data: np.ndarray


@dataclasses.dataclass(frozen=True)
class TwoColTable(table.Table):
  data: np.ndarray
  ints: np.ndarray


@dataclasses.dataclass(frozen=True)
class MultiDimColTable(table.Table):
  multi_data: np.ndarray


_NAN: Final[float] = float('nan')


class TableTest(parameterized.TestCase):

  def test_post_init_validation(self):
    with self.assertRaises(ValueError):
      _ = OneColTable(
          key=np.array([1, 2, 3], dtype=np.int64),
          data=np.array(['A'], dtype=object),
      )

  @parameterized.named_parameters(
      dict(
          testcase_name='numeric_column_str_predicate',
          col='ints',
          predicate='A',
      ),
      dict(
          testcase_name='character_column_numeric_predicate',
          col='data',
          predicate=(1, 2, 3),
      ),
  )
  def test_make_filter_mask_mismatched_value_type_raises_value_error(
      self, col, predicate
  ):
    tbl = TwoColTable(
        key=np.array([0, 1, 2], dtype=np.int64),
        data=np.array(['A', 'B', 'C'], dtype='<U1'),
        ints=np.array([3, 5, 7], dtype=np.int32),
    )
    with self.assertRaisesRegex(ValueError, f'Cannot filter column {col}'):
      tbl.make_filter_mask(**{col: predicate})

  def test_single_table(self):
    test_table = OneColTable(
        key=np.array([1, 2, 3, 5, 6, 8]),
        data=np.array(['A', 'B', 'C', 'D', 'E', 'F'], dtype=object),
    )

    with self.subTest('columns'):
      self.assertSequenceEqual(test_table.columns, ['key', 'data'])
    with self.subTest('size'):
      self.assertEqual(test_table.size, 6)
    with self.subTest('len'):
      self.assertLen(test_table, 6)
    with self.subTest('get_column_by_name'):
      np.testing.assert_array_equal(test_table['data'], test_table.data)

    mask = np.array([False, True, False, False, True, True])
    with self.subTest('get_column_with_mask'):
      np.testing.assert_array_equal(
          test_table['data', mask], np.array(['B', 'E', 'F'])
      )
    with self.subTest('apply_mask'):
      masked_test_table = test_table[mask]
      np.testing.assert_array_equal(
          masked_test_table.data, np.array(['B', 'E', 'F'])  # pyrefly: ignore[missing-attribute]
      )

    key_arr = np.array([1, 5, 6])
    with self.subTest('slice_into_column'):
      np.testing.assert_array_equal(
          test_table['data', key_arr], np.array(['A', 'D', 'E'])
      )
    with self.subTest('slice_table'):
      sliced_table = test_table[key_arr]
      np.testing.assert_array_equal(
          sliced_table.data, np.array(['A', 'D', 'E'])  # pyrefly: ignore[missing-attribute]
      )
    with self.subTest('apply_index'):
      masked_test_table = test_table.apply_index(
          np.array([0, 3, 4], dtype=np.int64)
      )
      np.testing.assert_array_equal(
          masked_test_table.data, np.array(['A', 'D', 'E'])
      )

    with self.subTest('get_row_by_key'):
      row = test_table.get_row_by_key(6)
      self.assertDictEqual(row, {'key': 6, 'data': 'E'})
      self.assertDictEqual(row, {'key': 6, 'data': 'E'})

    with self.subTest('get_row_by_index'):
      row = test_table.get_row_by_index(4)
      self.assertDictEqual(row, {'key': 6, 'data': 'E'})

    with self.subTest('get_row_by_key_renamed'):
      row = test_table.get_row_by_key(6, column_name_map={'data_col': 'data'})
      self.assertDictEqual(row, {'data_col': 'E'})

    with self.subTest('get_row_by_index_renamed'):
      row = test_table.get_row_by_index(4, column_name_map={'data_col': 'data'})
      self.assertDictEqual(row, {'data_col': 'E'})

  def test_tables_with_foreign_keys(self):
    @dataclasses.dataclass(frozen=True)
    class TableA(table.Table):
      data: np.ndarray

    @dataclasses.dataclass(frozen=True)
    class TableB(table.Table):
      table_a_key: np.ndarray
      data: np.ndarray

    table_a = TableA(
        key=np.array([1, 2, 3, 5, 6, 8]),
        data=np.array(['A', 'B', 'C', 'D', 'E', 'F'], dtype=object),
    )

    table_b = TableB(
        key=np.array([0, 1, 2, 3]),
        table_a_key=np.array([2, 2, 5, 8]),
        data=np.array([2.3, 5.2, 8.1, 3.0], dtype=np.float32),
    )

    with self.subTest('iterrows'):
      expected = [
          {'A_data': 'B', 'B_data': np.float32(2.3)},
          {'A_data': 'B', 'B_data': np.float32(5.2)},
          {'A_data': 'D', 'B_data': np.float32(8.1)},
          {'A_data': 'F', 'B_data': np.float32(3.0)},
      ]
      rows = list(
          table_b.iterrows(
              column_name_map={'B_data': 'data'},
              table_a_key=table_a.with_column_names({'A_data': 'data'}),
          )
      )
      self.assertSequenceEqual(rows, expected)

    with self.subTest('iterrows_with_keys'):
      expected = [
          {'A_data': 'B', 'B_data': np.float32(5.2)},
          {'A_data': 'D', 'B_data': np.float32(8.1)},
      ]
      rows = list(
          table_b.iterrows(
              row_keys=np.array([1, 2]),
              column_name_map={'B_data': 'data'},
              table_a_key=table_a.with_column_names({'A_data': 'data'}),
          )
      )
      self.assertSequenceEqual(rows, expected)

  def test_filter(self):
    test_table = TwoColTable(
        key=np.array([1, 2, 3, 5, 6, 8]),
        data=np.array(['A', 'B', 'C', 'D', 'E', 'F'], dtype=object),
        ints=np.array([3, 5, 1, 7, -1, 5]),
    )

    with self.subTest('filter_single_value'):
      filtered = test_table.filter(data='C')
      self.assertEqual(filtered.size, 1)
      np.testing.assert_array_equal(filtered.key, [3])
      np.testing.assert_array_equal(filtered.data, ['C'])
      np.testing.assert_array_equal(filtered.ints, [1])

    with self.subTest('filter_multiple_values'):
      filtered = test_table.filter(data=('B', 'E', 'F'))
      self.assertEqual(filtered.size, 3)
      np.testing.assert_array_equal(filtered.key, [2, 6, 8])
      np.testing.assert_array_equal(filtered.data, ['B', 'E', 'F'])
      np.testing.assert_array_equal(filtered.ints, [5, -1, 5])

    with self.subTest('filter_multiple_cols'):
      filtered = test_table.filter(data='B', ints=5)
      self.assertEqual(filtered.size, 1)
      np.testing.assert_array_equal(filtered.key, [2])
      np.testing.assert_array_equal(filtered.data, ['B'])
      np.testing.assert_array_equal(filtered.ints, [5])

    with self.subTest('filter_multiple_cols_and_vals'):
      filtered = test_table.filter(data=('B', 'C', 'D'), ints=(5, 6, 7))
      self.assertEqual(filtered.size, 2)
      np.testing.assert_array_equal(filtered.key, [2, 5])
      np.testing.assert_array_equal(filtered.data, ['B', 'D'])
      np.testing.assert_array_equal(filtered.ints, [5, 7])

    with self.subTest('filter_with_lambda'):
      filtered = test_table.filter(ints=lambda x: x < 4)
      self.assertEqual(filtered.size, 3)
      np.testing.assert_array_equal(filtered.key, [1, 3, 6])
      np.testing.assert_array_equal(filtered.data, ['A', 'C', 'E'])
      np.testing.assert_array_equal(filtered.ints, [3, 1, -1])

    with self.subTest('filter_with_mask'):
      mask = np.array([False, True, False, False, True, True])
      filtered = test_table.filter(mask)
      self.assertEqual(filtered.size, 3)
      np.testing.assert_array_equal(filtered.key, [2, 6, 8])
      np.testing.assert_array_equal(filtered.data, ['B', 'E', 'F'])
      np.testing.assert_array_equal(filtered.ints, [5, -1, 5])

    with self.subTest('filter_with_mask_and_values'):
      mask = np.array([False, True, False, False, True, True])
      filtered = test_table.filter(mask, data=('D', 'E', 'F'))
      self.assertEqual(filtered.size, 2)
      np.testing.assert_array_equal(filtered.key, [6, 8])
      np.testing.assert_array_equal(filtered.data, ['E', 'F'])
      np.testing.assert_array_equal(filtered.ints, [-1, 5])

    with self.subTest('filter_single_inverted'):
      filtered = test_table.filter(data='C', invert=True)
      self.assertEqual(filtered.size, 5)
      np.testing.assert_array_equal(filtered.key, [1, 2, 5, 6, 8])
      np.testing.assert_array_equal(filtered.data, ['A', 'B', 'D', 'E', 'F'])
      np.testing.assert_array_equal(filtered.ints, [3, 5, 7, -1, 5])

    with self.subTest('filter_multiple_inverted'):
      filtered = test_table.filter(
          data=('B', 'C', 'D'), ints=(5, 6, 7), invert=True
      )
      self.assertEqual(filtered.size, 4)
      np.testing.assert_array_equal(filtered.key, [1, 3, 6, 8])
      np.testing.assert_array_equal(filtered.data, ['A', 'C', 'E', 'F'])
      np.testing.assert_array_equal(filtered.ints, [3, 1, -1, 5])

    with self.subTest('filter_apply_by_element'):
      filtered = test_table.filter(
          data=lambda x: x.lower() == 'c', apply_per_element=True
      )
      self.assertEqual(filtered.size, 1)
      np.testing.assert_array_equal(filtered.key, [3])
      np.testing.assert_array_equal(filtered.data, ['C'])
      np.testing.assert_array_equal(filtered.ints, [1])

    with self.subTest('filter_all'):
      filtered = test_table.filter(data=lambda _: False, apply_per_element=True)
      self.assertEqual(filtered.size, 0)
      np.testing.assert_array_equal(filtered.key, [])
      np.testing.assert_array_equal(filtered.data, [])
      np.testing.assert_array_equal(filtered.ints, [])

  def test_copy_methods(self):
    test_table = TwoColTable(
        key=np.array([1, 2, 3, 5, 6, 8], dtype=np.int64),
        data=np.array(['A', 'B', 'C', 'D', 'E', 'F'], dtype=object),
        ints=np.array([3, 5, 1, 7, -1, 5], dtype=np.int64),
    )

    with self.subTest('copy'):
      copied = test_table.copy_and_update()
      self.assertIsNot(copied, test_table)
      for col in copied.columns:
        self.assertIs(
            copied.get_column(col), test_table.get_column(col), msg=col
        )

    with self.subTest('copy_and_update_invalid_column_name'):
      with self.assertRaisesRegex(ValueError, 'Invalid column names: '):
        _ = test_table.copy_and_update(bad_col_name=test_table.key)

    with self.subTest('copy_and_update_invalid_length'):
      with self.assertRaisesRegex(ValueError, 'All columns should have length'):
        _ = test_table.copy_and_update(ints=np.array([1, 2, 3], dtype=np.int64))

    with self.subTest('copy_and_update'):
      new_data = np.array(['U', 'V', 'W', 'X', 'Y', 'Z'], dtype=object)
      new_ints = np.array([0, 0, 1, 1, 2, 2], dtype=np.int64)
      updated = test_table.copy_and_update(data=new_data, ints=new_ints)
      self.assertIsNot(updated, test_table)
      self.assertIs(updated.data, new_data)
      self.assertIs(updated.ints, new_ints)
      self.assertIs(updated.key, test_table.key)

    with self.subTest('copy_and_remap'):
      remapped = test_table.copy_and_remap(
          key={2: 10, 5: 15, 15: 5},
          data={'A': 'B', 'B': 'A', 'Q': 'Z'},
      )
      self.assertIsNot(remapped, test_table)
      self.assertIs(remapped.ints, test_table.ints)
      np.testing.assert_array_equal(remapped.key, [1, 10, 3, 15, 6, 8])
      np.testing.assert_array_equal(
          remapped.data, ['B', 'A', 'C', 'D', 'E', 'F']
      )

    with self.subTest('copy_and_remap_invalid_column_name'):
      with self.assertRaisesRegex(ValueError, 'Invalid column names: '):
        _ = test_table.copy_and_remap(bad_col_name={})

  def test_copy_and_remap_empty_table(self):
    test_table = OneColTable(
        key=np.array([], dtype=np.int64),
        data=np.array([], dtype=object),
    )

    remapped = test_table.copy_and_remap(key={1: 5}, data={'a': 'A', 'b': 'B'})
    with self.subTest('copies_table'):
      self.assertEqual(remapped.size, 0)

    with self.subTest('invalid_column_name'):
      with self.assertRaisesRegex(ValueError, 'Invalid column names: '):
        _ = test_table.copy_and_remap(bad_col_name={})

  def test_make_filter_mask_float32_column(self):
    @dataclasses.dataclass(frozen=True, kw_only=True)
    class FloatColTable(table.Table):
      vals: np.ndarray

    tbl = FloatColTable(
        key=np.array([0, 1, 2], dtype=np.int64),
        vals=np.array([0.3, 0.5, 0.7], dtype=np.float32),
    )
    np.testing.assert_array_equal(
        tbl.make_filter_mask(vals=(0.3, 0.7)),
        np.array([True, False, True]),
        strict=True,
    )

  def test_make_filter_mask_bool_column_multi_value(self):
    @dataclasses.dataclass(frozen=True, kw_only=True)
    class BoolColTable(table.Table):
      flags: np.ndarray

    tbl = BoolColTable(
        key=np.array([0, 1], dtype=np.int64),
        flags=np.array([True, False], dtype=bool),
    )
    np.testing.assert_array_equal(
        tbl.make_filter_mask(flags=(True, False)),
        np.array([True, True]),
        strict=True,
    )
    np.testing.assert_array_equal(
        tbl.make_filter_mask(flags=tuple(str(i) for i in range(12))),
        np.array([False, False]),
        strict=True,
    )

  @parameterized.named_parameters(
      dict(testcase_name='object_int', dtype=object, predicate=1),
      dict(testcase_name='object_bool', dtype=object, predicate=True),
      dict(testcase_name='object_nan', dtype=object, predicate=_NAN),
      dict(
          testcase_name='object_nan_tuple',
          dtype=object,
          predicate=(_NAN, 'B'),
      ),
      dict(testcase_name='object_bytes', dtype=object, predicate=b'A'),
      dict(
          testcase_name='string_dtype_int',
          dtype=np.dtypes.StringDType(),
          predicate=(1, 2),
      ),
      dict(
          testcase_name='string_dtype_nan',
          dtype=np.dtypes.StringDType(),
          predicate=_NAN,
      ),
  )
  def test_make_filter_mask_string_or_object_invalid_predicate_raises(
      self, dtype, predicate
  ):
    tbl = OneColTable(
        key=np.array([0, 1], dtype=np.int64),
        data=np.array(['A', 'B'], dtype=dtype),
    )
    with self.assertRaisesRegex(ValueError, 'Cannot filter column data'):
      tbl.make_filter_mask(data=predicate)

  def test_make_filter_mask_object_column_multi_value_with_none(self):
    tbl = OneColTable(
        key=np.arange(4, dtype=np.int64),
        data=np.array(['A', 'B', None, _NAN], dtype=object),
    )
    np.testing.assert_array_equal(
        tbl.make_filter_mask(data=('B', None)),
        np.array([False, True, True, False]),
        strict=True,
    )

  @parameterized.named_parameters(
      dict(
          testcase_name='string',
          predicate='B',
          expected=[False, True, False],
      ),
      dict(
          testcase_name='none',
          predicate=None,
          expected=[False, False, True],
      ),
  )
  def test_make_filter_mask_object_column_single_value(
      self, predicate, expected
  ):
    tbl = OneColTable(
        key=np.arange(3, dtype=np.int64),
        data=np.array(['A', 'B', None], dtype=object),
    )
    np.testing.assert_array_equal(
        tbl.make_filter_mask(data=predicate), np.array(expected), strict=True
    )

  def test_make_filter_mask_empty_iterable_returns_all_false(self):
    tbl = OneColTable(
        key=np.array([0, 1, 2], dtype=np.int64),
        data=np.array(['A', 'B', 'C'], dtype=object),
    )
    np.testing.assert_array_equal(
        tbl.make_filter_mask(data=()),
        np.array([False, False, False]),
        strict=True,
    )

  def test_make_filter_mask_normalizes_callable_mask_shapes(self):
    tbl = OneColTable(
        key=np.array([0, 1, 2], dtype=np.int64),
        data=np.array(['A', 'B', 'C'], dtype=object),
    )
    read_only = np.array([True, False, True], dtype=bool)
    read_only.flags.writeable = False
    for pred, expected in (
        (lambda _: read_only, [True, False, True]),
        (lambda _: True, [True, True, True]),
        (lambda _: np.array([1, 0, 1]), [True, False, True]),
    ):
      np.testing.assert_array_equal(
          tbl.make_filter_mask(data=pred),
          np.array(expected, dtype=bool),
          strict=True,
      )
    with self.assertRaisesRegex(ValueError, 'broadcast'):
      tbl.make_filter_mask(data=lambda _: np.array([True, False]))

  def test_make_filter_mask_per_element_truthy_callable(self):
    tbl = TwoColTable(
        key=np.array([0, 1, 2], dtype=np.int64),
        data=np.array(['A', 'B', 'C'], dtype=object),
        ints=np.array([3, 5, 7], dtype=np.int32),
    )
    truthy_fn: Any = lambda x: x - 3
    np.testing.assert_array_equal(
        tbl.make_filter_mask(ints=truthy_fn, apply_per_element=True),
        np.array([False, True, True]),
        strict=True,
    )

  def test_make_filter_mask_multichar_string_and_generator_predicates(self):
    tbl = OneColTable(
        key=np.array([0, 1, 2], dtype=np.int64),
        data=np.array(['AB', 'A', 'B'], dtype=object),
    )
    for pred, expected in (
        ('AB', [True, False, False]),
        ((c for c in ('A', 'C')), [False, True, False]),
    ):
      np.testing.assert_array_equal(
          tbl.make_filter_mask(data=pred),
          np.array(expected, dtype=bool),
          strict=True,
      )

  def test_make_filter_mask_does_not_mutate_input_mask(self):
    full_data = np.array(['A', 'X', 'B', 'Y', 'C', 'Z'], dtype=object)
    tbl = TwoColTable(
        key=np.array([0, 1, 2], dtype=np.int64),
        data=full_data[::2],
        ints=np.array([3, 5, 7], dtype=np.int32),
    )
    input_mask = np.array([True, True, False], dtype=bool)
    input_mask.flags.writeable = False
    np.testing.assert_array_equal(
        tbl.make_filter_mask(
            mask=input_mask,
            data=(c for c in ('B', 'C')),
            ints=(5, 9),
        ),
        np.array([False, True, False]),
        strict=True,
    )

  def test_make_filter_mask_does_not_mutate_callable_masks(self):
    tbl = TwoColTable(
        key=np.array([0, 1, 2], dtype=np.int64),
        data=np.array(['A', 'B', 'C'], dtype=object),
        ints=np.array([3, 5, 7], dtype=np.int32),
    )
    read_only_mask = np.array([True, True, False], dtype=bool)
    read_only_mask.flags.writeable = False
    writeable_mask = np.array([False, True, True], dtype=bool)
    np.testing.assert_array_equal(
        tbl.make_filter_mask(
            data=lambda _: read_only_mask, ints=lambda _: writeable_mask
        ),
        np.array([False, True, False]),
        strict=True,
    )
    np.testing.assert_array_equal(
        writeable_mask, np.array([False, True, True]), strict=True
    )

  def test_make_filter_mask_multi_value_int_out_of_dtype_range(self):
    tbl = TwoColTable(
        key=np.array([0, 1, 2], dtype=np.int64),
        data=np.array(['A', 'B', 'C'], dtype=object),
        ints=np.array([3, 5, 7], dtype=np.int32),
    )
    np.testing.assert_array_equal(
        tbl.make_filter_mask(ints=(5, 10**12)),
        np.array([False, True, False]),
        strict=True,
    )

  @parameterized.named_parameters(
      dict(
          testcase_name='contiguous_keys',
          table_keys=[0, 1, 2],
          query_keys=np.array([2, 0, 1, 2], dtype=np.int32),
          expected=['C', 'A', 'B', 'C'],
      ),
      dict(
          testcase_name='sparse_unsorted_keys',
          table_keys=[100_000, 10, 5_000],
          query_keys=np.array([5_000, 100_000, 10], dtype=np.int64),
          expected=['C', 'A', 'B'],
      ),
      dict(
          testcase_name='sparse_empty_query',
          table_keys=[10, 5_000, 100_000],
          query_keys=np.array([], dtype=np.int64),
          expected=[],
      ),
  )
  def test_apply_array_to_column_avoids_dense_index_by_key(
      self, table_keys, query_keys, expected
  ):
    tbl = OneColTable(
        key=np.array(table_keys, dtype=np.int64),
        data=np.array(['A', 'B', 'C'][: len(table_keys)], dtype=object),
    )
    np.testing.assert_array_equal(
        tbl.apply_array_to_column('data', query_keys),
        np.array(expected, dtype=object),
        strict=True,
    )
    self.assertNotIn('index_by_key', tbl.__dict__)

  @parameterized.named_parameters(
      dict(
          testcase_name='float',
          bad_query=np.array([0.0, 1.0], dtype=np.float64),
      ),
      dict(testcase_name='uint', bad_query=np.array([2**63], dtype=np.uint64)),
  )
  def test_apply_array_to_column_non_signed_int_query_raises_index_error(
      self, bad_query
  ):
    tbl = OneColTable(
        key=np.array([0, 1], dtype=np.int64),
        data=np.array(['A', 'B'], dtype=object),
    )
    with self.assertRaisesRegex(IndexError, 'Keys'):
      tbl.apply_array_to_column('data', bad_query)

  @parameterized.named_parameters(
      dict(testcase_name='empty_table', table_keys=[], query_keys=[0]),
      dict(testcase_name='negative_key', table_keys=[0, 1, 2], query_keys=[-1]),
      dict(
          testcase_name='missing_sparse_key',
          table_keys=[10, 5_000, 100_000],
          query_keys=[5_001],
      ),
  )
  def test_apply_array_to_column_invalid_key_raises_index_error(
      self, table_keys, query_keys
  ):
    tbl = OneColTable(
        key=np.array(table_keys, dtype=np.int64),
        data=np.array(['A', 'B', 'C'][: len(table_keys)], dtype=object),
    )
    with self.assertRaisesRegex(IndexError, 'Keys'):
      tbl.apply_array_to_column('data', np.array(query_keys, dtype=np.int64))

  def test_iterrows_row_keys_avoids_dense_index_by_key(self):
    tbl = OneColTable(
        key=np.array([100_000, 10, 5_000], dtype=np.int64),
        data=np.array(['A', 'B', 'C'], dtype=object),
    )
    rows = list(
        tbl.iterrows(row_keys=np.array([5_000, 100_000], dtype=np.int64))
    )
    self.assertSequenceEqual(
        rows, [{'key': 5_000, 'data': 'C'}, {'key': 100_000, 'data': 'A'}]
    )
    self.assertNotIn('index_by_key', tbl.__dict__)

  def test_sliced_multidim_table_columns_remain_immutable(self):
    tbl = MultiDimColTable(
        key=np.array([10, 20, 30], dtype=np.int64),
        multi_data=np.array(
            [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float32
        ),
    )
    for sliced in (
        tbl.apply_array(np.array([True, False, True])),
        tbl.apply_array(np.array([10, 30], dtype=np.int64)),
        tbl.apply_index(np.array([0, 2], dtype=np.int64)),
    ):
      self.assertFalse(sliced.key.flags.writeable)
      self.assertFalse(sliced.multi_data.flags.writeable)
      np.testing.assert_array_equal(
          sliced.multi_data,
          np.array([[1.0, 3.0], [4.0, 6.0]], dtype=np.float32),
          strict=True,
      )

  @parameterized.named_parameters(
      dict(
          testcase_name='empty_spans_preserve_multidim_shape',
          total_size=5,
          starts=[2, 4],
          ends=[2, 4],
          expected_indices=[],
      ),
      dict(
          testcase_name='adjacent_spans_merge_to_full_table',
          total_size=5,
          starts=[0, 2],
          ends=[2, 5],
          expected_indices=[0, 1, 2, 3, 4],
      ),
      dict(
          testcase_name='single_interior_span',
          total_size=6,
          starts=[1],
          ends=[4],
          expected_indices=[1, 2, 3],
      ),
      dict(
          testcase_name='disjoint_spans_concat',
          total_size=6,
          starts=[0, 3],
          ends=[2, 5],
          expected_indices=[0, 1, 3, 4],
      ),
      dict(
          testcase_name='many_disjoint_spans_concat',
          total_size=70,
          starts=list(range(0, 70, 2)),
          ends=list(range(1, 71, 2)),
          expected_indices=list(range(0, 70, 2)),
      ),
  )
  def test_slice_spans(self, total_size, starts, ends, expected_indices):
    tbl = MultiDimColTable(
        key=np.arange(total_size, dtype=np.int64) * 10,
        multi_data=np.arange(2 * total_size, dtype=np.float32).reshape(
            2, total_size
        ),
    )
    sliced = tbl.slice_spans(
        np.array(starts, dtype=np.int64), np.array(ends, dtype=np.int64)
    )
    expected = tbl.apply_index(np.array(expected_indices, dtype=np.int64))
    np.testing.assert_array_equal(sliced.key, expected.key, strict=True)
    np.testing.assert_array_equal(
        sliced.multi_data, expected.multi_data, strict=True
    )
    self.assertFalse(sliced.key.flags.writeable)
    self.assertFalse(sliced.multi_data.flags.writeable)

  def test_slice_spans_adjacent_full_cover_returns_self(self):
    tbl = MultiDimColTable(
        key=np.arange(5, dtype=np.int64),
        multi_data=np.zeros((1, 5), dtype=np.float32),
    )
    self.assertIs(
        tbl.slice_spans(
            np.array([0, 2], dtype=np.int64), np.array([2, 5], dtype=np.int64)
        ),
        tbl,
    )

  def test_slice_spans_single_span_copies_columns(self):
    tbl = MultiDimColTable(
        key=np.array([0, 1, 2], dtype=np.int64),
        multi_data=np.array([[1.0, 2.0, 3.0]], dtype=np.float32),
    )
    sliced = tbl.slice_spans(
        np.array([1], dtype=np.int64), np.array([2], dtype=np.int64)
    )
    self.assertIsNone(sliced.multi_data.base)

  def test_slice_range_full_cover_returns_self(self):
    empty = MultiDimColTable(
        key=np.empty((0,), dtype=np.int64),
        multi_data=np.empty((1, 0), dtype=np.float32),
    )
    self.assertIs(empty.slice_range(start=0, end=0), empty)

    tbl = MultiDimColTable(
        key=np.arange(3, dtype=np.int64),
        multi_data=np.array([[1.0, 2.0, 3.0]], dtype=np.float32),
    )
    self.assertIs(tbl.slice_range(start=0, end=3), tbl)

  @parameterized.named_parameters(
      dict(testcase_name='interior', start=1, end=3),
      dict(testcase_name='empty_subrange', start=2, end=2),
  )
  def test_slice_range_copies_columns_and_preserves_multidim_shape(
      self, start, end
  ):
    tbl = MultiDimColTable(
        key=np.array([10, 20, 30, 40], dtype=np.int64),
        multi_data=np.arange(8, dtype=np.float32).reshape(2, 4),
    )
    sliced = tbl.slice_range(start=start, end=end)
    np.testing.assert_array_equal(sliced.key, tbl.key[start:end], strict=True)
    np.testing.assert_array_equal(
        sliced.multi_data, tbl.multi_data[:, start:end], strict=True
    )
    self.assertFalse(sliced.key.flags.writeable)
    self.assertFalse(sliced.multi_data.flags.writeable)
    self.assertIsNone(sliced.multi_data.base)

  @parameterized.named_parameters(
      dict(testcase_name='negative_start', start=-1, end=2),
      dict(testcase_name='start_after_end', start=3, end=2),
      dict(testcase_name='end_past_size', start=0, end=5),
  )
  def test_slice_range_invalid_bounds_raises(self, start, end):
    tbl = OneColTable(
        key=np.array([0, 1, 2], dtype=np.int64),
        data=np.array(['A', 'B', 'C'], dtype=object),
    )
    with self.assertRaises(ValueError):
      tbl.slice_range(start=start, end=end)

  @parameterized.named_parameters(
      dict(testcase_name='scattered_rows', kept_rows=slice(None, None, 3)),
      dict(testcase_name='long_run', kept_rows=slice(5, 30)),
  )
  def test_apply_array_bool_mask(self, kept_rows):
    tbl = MultiDimColTable(
        key=np.arange(40, dtype=np.int64) * 10,
        multi_data=np.arange(80, dtype=np.float32).reshape(2, 40),
    )
    mask = np.zeros(40, dtype=bool)
    mask[kept_rows] = True
    sliced = tbl.apply_array(mask)
    np.testing.assert_array_equal(sliced.key, tbl.key[mask], strict=True)
    np.testing.assert_array_equal(
        sliced.multi_data, tbl.multi_data[:, mask], strict=True
    )

  def test_apply_array_wrong_length_bool_mask_raises_index_error(self):
    tbl = OneColTable(
        key=np.array([0, 1, 2], dtype=np.int64),
        data=np.array(['A', 'B', 'C'], dtype=object),
    )
    with self.assertRaises(IndexError):
      tbl.apply_array(np.array([True, False]))

  @parameterized.named_parameters(
      dict(
          testcase_name='str_dtype',
          data=np.array(['X', 'Y', 'Z'], dtype=str),
          expected=np.array(['Z', 'Y'], dtype='<U1'),
      ),
      dict(
          testcase_name='strided_object',
          data=np.array(['A', 'skip', 'B', 'skip', 'C'], dtype=object)[::2],
          expected=np.array(['C', 'B'], dtype=object),
      ),
  )
  def test_apply_array_to_column_fallback_columns(self, data, expected):
    tbl = OneColTable(
        key=np.array([100_000, 10, 5_000], dtype=np.int64), data=data
    )
    gathered = tbl.apply_array_to_column(
        'data', np.array([5_000, 10], dtype=np.int64)
    )
    np.testing.assert_array_equal(gathered, expected, strict=True)
    self.assertNotIn('index_by_key', tbl.__dict__)


@dataclasses.dataclass(frozen=True)
class Residues(table.Table):
  id: np.ndarray


@dataclasses.dataclass(frozen=True)
class Atoms(table.Table):
  name: np.ndarray
  coords: np.ndarray  # 3 x N
  res_key: np.ndarray


@dataclasses.dataclass(frozen=True)
class Bonds(table.Table):
  from_atom_key: np.ndarray
  dest_atom_key: np.ndarray


@dataclasses.dataclass(frozen=True)
class Polymer(table.Database):
  residues: Residues
  atoms: Atoms
  bonds: Bonds

  tables: ClassVar[Collection[str]] = ('residues', 'atoms', 'bonds')  # pyrefly: ignore[bad-override]
  foreign_keys: ClassVar[Mapping[str, Collection[tuple[str, str]]]] = {  # pyrefly: ignore[bad-override]
      'atoms': (('res_key', 'residues'),),
      'bonds': (
          ('from_atom_key', 'atoms'),
          ('dest_atom_key', 'atoms'),
      ),
  }

  def get_table(self, table_name: str) -> table.Table:
    return getattr(self, table_name)

  def copy_and_update(self, **changes) -> Self:
    return dataclasses.replace(self, **changes)


class DatabaseTest(absltest.TestCase):

  def test_concat_fails_on_empty(self):
    with self.assertRaisesRegex(ValueError, 'Need at least one'):
      table.concat_databases([])

  def test_concat_fails_on_different_types(self):
    @dataclasses.dataclass(frozen=True)
    class OtherDB(table.Database):
      atoms: Atoms
      tables: ClassVar[Collection[str]] = ('atoms',)  # pyrefly: ignore[bad-override]
      foreign_keys: ClassVar[Mapping[str, Collection[tuple[str, str]]]] = {}  # pyrefly: ignore[bad-override]
      get_table = getattr  # pyrefly: ignore[bad-override]
      copy_and_update = dataclasses.replace

    with self.assertRaisesRegex(
        ValueError, 'All `dbs` must be of the same type.'
    ):
      table.concat_databases([
          Polymer(
              residues=Residues(
                  key=np.array([], dtype=np.int64),
                  id=np.array([], dtype=object),
              ),
              atoms=Atoms(
                  key=np.array([], dtype=np.int64),
                  name=np.array([], dtype=object),
                  coords=np.array([[], [], []], dtype=np.float32),
                  res_key=np.array([], dtype=np.int64),
              ),
              bonds=Bonds(
                  key=np.array([], dtype=np.int64),
                  from_atom_key=np.array([], dtype=np.int64),
                  dest_atom_key=np.array([], dtype=np.int64),
              ),
          ),
          OtherDB(
              atoms=Atoms(
                  key=np.array([], dtype=np.int64),
                  name=np.array([], dtype=object),
                  coords=np.array([[], [], []], dtype=np.float32),
                  res_key=np.array([], dtype=np.int64),
              ),
          ),
      ])

  def test_concat(self):
    p1 = Polymer(
        residues=Residues(
            key=np.array([1, 2], dtype=np.int64),
            id=np.array([1, 2], dtype=object),
        ),
        atoms=Atoms(
            key=np.array([1, 2, 3, 4], dtype=np.int64),
            name=np.array(['a', 'b', 'c', 'd'], dtype=object),
            coords=np.array(
                [
                    [0.0, 1.0, 2.0, 3.0],
                    [0.1, 1.1, 2.1, 3.1],
                    [0.2, 1.2, 2.2, 3.2],
                ],
                dtype=np.float32,
            ),
            res_key=np.array([1, 2, 1, 2], dtype=np.int64),
        ),
        bonds=Bonds(
            key=np.array([1, 2, 3], dtype=np.int64),
            from_atom_key=np.array([1, 2, 1], dtype=np.int64),
            dest_atom_key=np.array([3, 4, 2], dtype=np.int64),
        ),
    )
    p2 = Polymer(
        residues=Residues(
            key=np.array([1, 2], dtype=np.int64),
            id=np.array([5, 6], dtype=object),
        ),
        atoms=Atoms(
            key=np.array([1, 2, 3], dtype=np.int64),
            name=np.array(['x', 'y', 'z'], dtype=object),
            coords=np.array(
                [
                    [10.0, 11.0, 12.0],
                    [10.1, 11.1, 12.1],
                    [10.2, 11.2, 12.2],
                ],
                dtype=np.float32,
            ),
            res_key=np.array([2, 2, 1], dtype=np.int64),
        ),
        bonds=Bonds(
            key=np.array([1, 2], dtype=np.int64),
            from_atom_key=np.array([2, 3], dtype=np.int64),
            dest_atom_key=np.array([3, 1], dtype=np.int64),
        ),
    )
    empty = Polymer(
        residues=Residues(
            key=np.array([1, 2], dtype=np.int64),
            id=np.array([8, 9], dtype=object),
        ),
        atoms=Atoms(
            key=np.array([], dtype=np.int64),
            name=np.array([], dtype=object),
            coords=np.array([[], [], []], dtype=np.float32),
            res_key=np.array([], dtype=np.int64),
        ),
        bonds=Bonds(
            key=np.array([], dtype=np.int64),
            from_atom_key=np.array([], dtype=np.int64),
            dest_atom_key=np.array([], dtype=np.int64),
        ),
    )
    concatted = table.concat_databases([p1, empty, p2])
    with self.subTest('plain_columns_concatenate'):
      np.testing.assert_array_equal(
          concatted.residues.id,
          [1, 2, 8, 9, 5, 6],
      )
      np.testing.assert_array_equal(
          concatted.atoms.name, ['a', 'b', 'c', 'd', 'x', 'y', 'z']
      )
    with self.subTest('foreign_keys_concatenate'):
      np.testing.assert_array_equal(
          concatted.residues.apply_array(concatted.atoms.res_key).id,
          [1, 2, 1, 2, 6, 6, 5],
      )
      np.testing.assert_array_almost_equal(
          concatted.atoms.coords,
          [
              [0.0, 1.0, 2.0, 3.0, 10.0, 11.0, 12.0],
              [0.1, 1.1, 2.1, 3.1, 10.1, 11.1, 12.1],
              [0.2, 1.2, 2.2, 3.2, 10.2, 11.2, 12.2],
          ],
      )
      np.testing.assert_array_equal(
          concatted.atoms.apply_array(concatted.bonds.from_atom_key).name,
          ['a', 'b', 'a', 'y', 'z'],
      )
      np.testing.assert_array_equal(
          concatted.atoms.apply_array(concatted.bonds.dest_atom_key).name,
          ['c', 'd', 'b', 'z', 'x'],
      )
    with self.subTest('key_is_unique'):
      for table_name in Polymer.tables:
        self.assertLen(
            set(concatted.get_table(table_name).key),
            len(concatted.get_table(table_name).key),
        )


if __name__ == '__main__':
  absltest.main()
