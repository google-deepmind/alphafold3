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


if __name__ == '__main__':
  absltest.main()
