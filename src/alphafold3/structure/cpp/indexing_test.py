# Copyright 2026 DeepMind Technologies Limited
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

from concurrent import futures
import sys

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.cpp import indexing
import numpy as np


class IndexingTest(parameterized.TestCase):

  @parameterized.named_parameters(
      dict(
          testcase_name='identity',
          table_keys=[0, 1, 2],
          keys=[2, 0, 1],
          expected=[2, 0, 1],
      ),
      dict(
          testcase_name='shifted_contiguous',
          table_keys=[7, 8, 9],
          keys=[9, 7, 8],
          expected=[2, 0, 1],
      ),
      dict(
          testcase_name='shifted_contiguous_chunked',
          table_keys=list(range(10, 310)),
          keys=list(range(309, 9, -1)),
          expected=list(range(299, -1, -1)),
      ),
      dict(
          testcase_name='singleton_large',
          table_keys=[10_000_000],
          keys=[10_000_000, 10_000_000],
          expected=[0, 0],
      ),
      dict(
          testcase_name='sorted_binary_search_with_runs',
          table_keys=list(range(0, 200, 2)),
          keys=[198, 198, 0, 64],
          expected=[99, 99, 0, 32],
      ),
      dict(
          testcase_name='shifted_lut',
          table_keys=[10, 14, 19],
          keys=[19, 10, 14, 10],
          expected=[2, 0, 1, 0],
      ),
      dict(
          testcase_name='sparse_hash_sorted',
          table_keys=[10, 5_000_000, 10_000_000],
          keys=[10_000_000, 10, 5_000_000],
          expected=[2, 0, 1],
      ),
      dict(
          testcase_name='sparse_hash_unsorted_with_runs',
          table_keys=[10_000_000, 10, 5_000_000],
          keys=[5_000_000, 5_000_000, 10_000_000, 10],
          expected=[2, 2, 0, 1],
      ),
      dict(
          testcase_name='duplicate_table_keys_last_wins',
          table_keys=[0, 0, 2],
          keys=[0, 2],
          expected=[1, 2],
      ),
      dict(
          testcase_name='duplicate_sparse_last_wins',
          table_keys=[10_000_000, 10_000_000, 0],
          keys=[10_000_000, 0],
          expected=[1, 2],
      ),
      dict(
          testcase_name='negative_table_keys',
          table_keys=[-100, -50, 25],
          keys=[25, -100, -50],
          expected=[2, 0, 1],
      ),
      dict(
          testcase_name='extreme_int64_bounds',
          table_keys=[-(2**63), 0, 2**63 - 1],
          keys=[2**63 - 1, -(2**63), 0],
          expected=[2, 0, 1],
      ),
      dict(
          testcase_name='empty_query_non_empty_table',
          table_keys=[1, 2],
          keys=[],
          expected=[],
      ),
      dict(
          testcase_name='empty_query_empty_table',
          table_keys=[],
          keys=[],
          expected=[],
      ),
  )
  def test_index_of(self, table_keys, keys, expected):
    np.testing.assert_array_equal(
        indexing.index_of(
            np.array(keys, dtype=np.int64),
            np.array(table_keys, dtype=np.int64),
        ),
        np.array(expected, dtype=np.int64),
        strict=True,
    )

  @parameterized.named_parameters(
      dict(
          testcase_name='int16_query_int32_table',
          keys=np.array([30, 10], dtype=np.int16),
          table_keys=np.array([10, 20, 30], dtype=np.int32),
          expected=np.array([2, 0], dtype=np.int64),
      ),
      dict(
          testcase_name='negative_strided',
          keys=np.array([30, 99, 10, 99, 20], dtype=np.int64)[::-2],
          table_keys=np.array([30, 99, 20, 99, 10], dtype=np.int64)[::-2],
          expected=np.array([1, 0, 2], dtype=np.int64),
      ),
  )
  def test_index_of_dtypes_and_strided_layouts(
      self, keys, table_keys, expected
  ):
    np.testing.assert_array_equal(
        indexing.index_of(keys, table_keys), expected, strict=True
    )

  @parameterized.parameters(np.int64, np.longlong)
  def test_index_of_identity_borrows_read_only_int64(self, dtype):
    read_only_keys = np.arange(300, dtype=dtype)[::-1].copy()
    read_only_keys.flags.writeable = False
    borrowed = indexing.index_of(read_only_keys, np.arange(300, dtype=dtype))
    self.assertIs(borrowed, read_only_keys)

  def test_index_of_identity_copies_writeable_int64(self):
    writeable_keys = np.array([2, 0, 1], dtype=np.int64)
    copied = indexing.index_of(
        writeable_keys, np.array([0, 1, 2], dtype=np.int64)
    )
    self.assertIsNot(copied, writeable_keys)
    self.assertTrue(copied.flags.writeable)
    np.testing.assert_array_equal(copied, writeable_keys, strict=True)

  @parameterized.parameters(0, 10)
  def test_index_of_concurrent_threads(self, start):
    table_keys = np.arange(start, start + 300, dtype=np.int64)
    keys = table_keys[::-1].copy()
    keys.flags.writeable = False
    expected = np.arange(299, -1, -1, dtype=np.int64)
    with futures.ThreadPoolExecutor(max_workers=8) as pool:
      results = list(
          pool.map(lambda _: indexing.index_of(keys, table_keys), range(32))
      )
    for result in results:
      np.testing.assert_array_equal(result, expected, strict=True)

  @parameterized.named_parameters(
      dict(
          testcase_name='empty_table',
          keys=np.array([0], dtype=np.int64),
          table_keys=np.array([], dtype=np.int64),
      ),
      dict(
          testcase_name='missing_identity_negative',
          keys=np.array([-1], dtype=np.int64),
          table_keys=np.array([0, 1, 2], dtype=np.int64),
          read_only=True,
      ),
      dict(
          testcase_name='missing_identity_chunked',
          keys=np.insert(np.zeros(299, dtype=np.int64), 100, 300),
          table_keys=np.arange(300, dtype=np.int64),
          read_only=True,
      ),
      dict(
          testcase_name='missing_identity_prefix_spoof',
          keys=np.array([150], dtype=np.int64),
          table_keys=np.concatenate([
              np.arange(150, dtype=np.int64),
              [0],
              np.arange(151, 300, dtype=np.int64),
          ]),
          read_only=True,
      ),
      dict(
          testcase_name='missing_shifted_contiguous_chunked',
          keys=np.insert(np.full(299, 7, dtype=np.int64), 100, 6),
          table_keys=np.arange(7, 307, dtype=np.int64),
      ),
      dict(
          testcase_name='missing_binary_search_gap',
          keys=np.array([3], dtype=np.int64),
          table_keys=np.arange(0, 200, 2, dtype=np.int64),
      ),
      dict(
          testcase_name='missing_lut_gap',
          keys=np.array([12], dtype=np.int64),
          table_keys=np.array([10, 14, 19], dtype=np.int64),
      ),
      dict(
          testcase_name='missing_duplicate_gap',
          keys=np.array([1], dtype=np.int64),
          table_keys=np.array([0, 0, 2], dtype=np.int64),
      ),
      dict(
          testcase_name='missing_sparse_hash_in_span',
          keys=np.array([5], dtype=np.int64),
          table_keys=np.array([0, 10_000_000], dtype=np.int64),
      ),
      dict(
          testcase_name='missing_sparse_hash_out_of_span',
          keys=np.array([-1], dtype=np.int64),
          table_keys=np.array([0, 10_000_000], dtype=np.int64),
      ),
      dict(
          testcase_name='uint64_keys',
          keys=np.array([0], dtype=np.uint64),
          table_keys=np.array([0, 1], dtype=np.int64),
      ),
      dict(
          testcase_name='float_empty_keys',
          keys=np.array([], dtype=np.float64),
          table_keys=np.array([0, 1], dtype=np.int64),
      ),
      dict(
          testcase_name='uint64_table_keys',
          keys=np.array([0], dtype=np.int64),
          table_keys=np.array([0, 1], dtype=np.uint64),
      ),
      dict(
          testcase_name='0d_keys',
          keys=np.array(0, dtype=np.int64),
          table_keys=np.array([0, 1], dtype=np.int64),
      ),
      dict(
          testcase_name='2d_keys',
          keys=np.array([[0, 1]], dtype=np.int64),
          table_keys=np.array([0, 1], dtype=np.int64),
      ),
      dict(
          testcase_name='2d_table_keys',
          keys=np.array([0], dtype=np.int64),
          table_keys=np.array([[0, 1]], dtype=np.int64),
      ),
  )
  def test_index_of_invalid_raises_index_error(
      self, keys, table_keys, read_only=False
  ):
    if read_only:
      keys = keys.copy()
      keys.flags.writeable = False
    with self.assertRaisesRegex(IndexError, 'Keys'):
      indexing.index_of(keys, table_keys)

  @parameterized.named_parameters(
      dict(
          testcase_name='list_keys',
          keys=[0],
          table_keys=np.array([0], dtype=np.int64),
      ),
      dict(
          testcase_name='list_table_keys',
          keys=np.array([0], dtype=np.int64),
          table_keys=[0],
      ),
  )
  def test_index_of_non_ndarray_raises_type_error(self, keys, table_keys):
    with self.assertRaisesRegex(TypeError, 'incompatible function arguments'):
      indexing.index_of(keys, table_keys)

  @parameterized.named_parameters(
      dict(testcase_name='all_false', mask=[False, False, False]),
      dict(testcase_name='sparse', mask=[False] * 17 + [True] + [False] * 9),
      dict(testcase_name='dense', mask=[True, False, True, True, False]),
  )
  def test_nonzero_indices(self, mask):
    arr = np.array(mask, dtype=bool)
    np.testing.assert_array_equal(
        indexing.nonzero_indices(arr),
        np.flatnonzero(arr),
        strict=True,
    )

  @parameterized.parameters(0, 5)
  def test_nonzero_indices_all_true_returns_none(self, size):
    self.assertIsNone(indexing.nonzero_indices(np.ones(size, dtype=bool)))

  @parameterized.product(
      case=(
          (list(range(300)), list(range(299, -1, -1))),
          ([7, 8, 9], [9, 7, 8]),
          (list(range(0, 200, 2)), [198, 198, 0, 64]),
          ([10, 14, 19], [19, 10, 14, 10]),
          ([10_000_000, 10, 5_000_000], [5_000_000, 5_000_000, 10_000_000, 10]),
          ([0, 0, 2], [0, 2]),
          ([1, 2], []),
          ([], []),
      ),
      dtype=(np.int32, object),
  )
  def test_take_by_key_regimes(self, case, dtype):
    table_keys = np.array(case[0], dtype=np.int64)
    keys = np.array(case[1], dtype=np.int64)
    values = (
        np.arange(table_keys.size, dtype=np.int32)
        if dtype is np.int32
        else np.array([f'V{i}' for i in range(table_keys.size)], dtype=object)
    )
    expected = values[indexing.index_of(keys, table_keys)]
    np.testing.assert_array_equal(
        indexing.take_by_key(values, keys, table_keys), expected, strict=True
    )

  @parameterized.named_parameters(
      (
          'pod_1b_bool',
          np.array([True, False, True], dtype=bool),
          np.array([30, 20, 10], dtype=np.int64),
          np.array([10, 20, 30], dtype=np.int64),
      ),
      (
          'pod_2b_int16',
          np.array([-5, 12, 99], dtype=np.int16),
          np.array([30, 10, 20, 30], dtype=np.int32),
          np.array([10, 20, 30], dtype=np.int64),
      ),
      (
          'pod_2b_float16',
          np.array([1.5, -2.5, 3.25], dtype=np.float16),
          np.array([20, 10], dtype=np.int16),
          np.array([10, 20, 30], dtype=np.int32),
      ),
      (
          'pod_4b_float32',
          np.array([1.5, -2.5, 3.25], dtype=np.float32),
          np.array([20, 10], dtype=np.int16),
          np.array([10, 20, 30], dtype=np.int32),
      ),
      (
          'pod_8b_float64',
          np.array([1.5, -2.5], dtype=np.float64),
          np.array([1, 0, 1], dtype=np.longlong),
          np.array([0, 1], dtype=np.longlong),
      ),
      (
          'neg_strided_keys',
          np.array([10, 20, 30], dtype=np.int32),
          np.array([30, 99, 10, 99, 20], dtype=np.int64)[::-2],
          np.array([30, 99, 20, 99, 10], dtype=np.int64)[::-2],
      ),
  )
  def test_take_by_key_dtypes(self, values, keys, table_keys):
    np.testing.assert_array_equal(
        indexing.take_by_key(values, keys, table_keys),
        values[indexing.index_of(keys, table_keys)],
        strict=True,
    )

  def test_take_by_key_object_refcounts_conserved(self):
    item_a = object()
    item_b = object()
    values = np.array([item_a, item_b], dtype=object)
    table_keys = np.array([10, 20], dtype=np.int64)
    base_a = sys.getrefcount(item_a)
    base_b = sys.getrefcount(item_b)

    gathered = indexing.take_by_key(
        values, np.array([10, 10, 20], dtype=np.int64), table_keys
    )
    self.assertEqual(sys.getrefcount(item_a), base_a + 2)
    self.assertEqual(sys.getrefcount(item_b), base_b + 1)
    del gathered
    self.assertEqual(sys.getrefcount(item_a), base_a)
    self.assertEqual(sys.getrefcount(item_b), base_b)

    with self.assertRaisesRegex(IndexError, 'Keys missing'):
      indexing.take_by_key(
          values, np.array([10, 99, 20], dtype=np.int64), table_keys
      )
    self.assertEqual(sys.getrefcount(item_a), base_a)
    self.assertEqual(sys.getrefcount(item_b), base_b)

  @parameterized.parameters(np.int32, object)
  def test_take_by_key_concurrent_threads(self, dtype):
    table_keys = np.arange(10, 310, dtype=np.int64)
    keys = table_keys[::-1].copy()
    values = (
        np.arange(table_keys.size, dtype=np.int32)
        if dtype is np.int32
        else np.array([f'V{i}' for i in range(table_keys.size)], dtype=object)
    )
    expected = values[::-1]
    with futures.ThreadPoolExecutor(max_workers=8) as pool:
      results = list(
          pool.map(
              lambda _: indexing.take_by_key(values, keys, table_keys),
              range(32),
          )
      )
    for result in results:
      np.testing.assert_array_equal(result, expected, strict=True)

  @parameterized.named_parameters(
      (
          'empty_table',
          np.array([], dtype=np.int32),
          np.array([0], dtype=np.int64),
          np.array([], dtype=np.int64),
          IndexError,
      ),
      (
          'missing_key',
          np.array(['A', 'B'], dtype=object),
          np.array([0, 5], dtype=np.int64),
          np.array([0, 1], dtype=np.int64),
          IndexError,
      ),
      (
          '2d_keys',
          np.array([1, 2], dtype=np.int32),
          np.array([[0, 1]], dtype=np.int64),
          np.array([0, 1], dtype=np.int64),
          IndexError,
      ),
      (
          '2d_values',
          np.array([[1, 2]], dtype=np.int32),
          np.array([0], dtype=np.int64),
          np.array([0, 1], dtype=np.int64),
          ValueError,
      ),
      (
          'length_mismatch_empty_keys',
          np.array([1], dtype=np.int32),
          np.array([], dtype=np.int64),
          np.array([0, 1], dtype=np.int64),
          ValueError,
      ),
      (
          'strided_values',
          np.array([1, 2, 3, 4], dtype=np.int32)[::2],
          np.array([0], dtype=np.int64),
          np.array([0, 1], dtype=np.int64),
          ValueError,
      ),
      (
          'unaligned_values',
          np.frombuffer(
              b'\x00' + np.array([1, 2], dtype=np.int32).tobytes(),
              dtype=np.int32,
              count=2,
              offset=1,
          ),
          np.array([0], dtype=np.int64),
          np.array([0, 1], dtype=np.int64),
          ValueError,
      ),
      (
          'list_values',
          [1, 2],
          np.array([0], dtype=np.int64),
          np.array([0, 1], dtype=np.int64),
          TypeError,
      ),
      (
          'byte_swapped_values',
          np.array([1, 2], dtype='>i4'),
          np.array([0], dtype=np.int64),
          np.array([0, 1], dtype=np.int64),
          TypeError,
      ),
      (
          'unicode_values_empty_keys',
          np.array(['A', 'B'], dtype='<U2'),
          np.array([], dtype=np.int64),
          np.array([0, 1], dtype=np.int64),
          TypeError,
      ),
      (
          'string_dtype_values',
          np.array(['A', 'B'], dtype=np.dtype('T')),
          np.array([0], dtype=np.int64),
          np.array([0, 1], dtype=np.int64),
          TypeError,
      ),
      (
          'complex128_values',
          np.array([1.0 + 2.0j, 3.0], dtype=np.complex128),
          np.array([0], dtype=np.int64),
          np.array([0, 1], dtype=np.int64),
          TypeError,
      ),
  )
  def test_take_by_key_invalid_errors(self, values, keys, table_keys, exc):
    with self.assertRaises(exc):
      indexing.take_by_key(values, keys, table_keys)


if __name__ == '__main__':
  absltest.main()
