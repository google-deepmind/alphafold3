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

from concurrent import futures
import sys

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.cpp import membership
import numpy as np


def _prevent_intern(s: str) -> str:
  """Prevents Python from interning the given string."""
  return ''.join(s)


_INT_DTYPES = (
    np.int64,
    np.int32,
    np.int16,
    np.int8,
    np.uint32,
    np.uint16,
    np.uint8,
)


class MembershipTest(parameterized.TestCase):

  @parameterized.parameters(
      ([], set(), False, []),
      ([], {1, 2}, False, []),
      ([1, 2, 3], set(), False, [False, False, False]),
      ([1, 2, 3], {2}, False, [False, True, False]),
      ([1, 2, 3], {1, 2}, False, [True, True, False]),
      ([1, 2, 3], {1, 2, 3}, False, [True, True, True]),
      ([1, 2, 2], {1, 2, 3}, False, [True, True, True]),
      ([2, 1, 3], {1, 2, 3}, False, [True, True, True]),
      ([], set(), True, []),
      ([], {1, 2}, True, []),
      ([1, 2, 3], set(), True, [True, True, True]),
      ([1, 2, 3], {2}, True, [True, False, True]),
      ([1, 2, 3], {1, 2}, True, [False, False, True]),
      ([1, 2, 3], {1, 2, 3}, True, [False, False, False]),
      ([1, 2, 2], {1, 2, 3}, True, [False, False, False]),
      ([2, 1, 3], {1, 2, 3}, True, [False, False, False]),
      ([2, 1, 3], {1, 2, 3}, True, [False, False, False], np.int32),
      ([], set(), True, [], np.int32),
  )
  def test_isin(self, array, test_elements, invert, exp_result, dtype=np.int64):
    result = membership.isin(
        np.array(array, dtype=dtype), test_elements, invert=invert
    )
    self.assertSequenceEqual(result.tolist(), exp_result)

  @parameterized.parameters(*_INT_DTYPES)
  def test_isin_multi_dim(self, dtype):
    result = membership.isin(
        np.array([[1, 2, 0], [3, 4, 0]], dtype=dtype), {1, 4}
    )
    self.assertSequenceEqual(
        result.tolist(), [[True, False, False], [False, True, False]]
    )

  @parameterized.parameters(
      *((dtype, 42, 1) for dtype in _INT_DTYPES),
      (object, 'A', 'B'),
      (np.dtype('T'), 'A', 'B'),
  )
  def test_isin_zero_dim(self, dtype, found, missing):
    result = membership.isin(np.array(found, dtype=dtype), {missing, found})
    self.assertEqual(result.tolist(), True)

  @parameterized.parameters(
      ([], set(), False, []),
      ([], {'A', 'B'}, False, []),
      ([None, 'A'], {'A', 'B', 'C'}, False, [False, True]),
      (['A', 'B', 'C'], set(), False, [False, False, False]),
      (['A', 'B', 'C'], {'B'}, False, [False, True, False]),
      (['A', 'B', 'C'], {'A', 'B'}, False, [True, True, False]),
      (['A', 'B', 'C'], {'A', 'B', 'C'}, False, [True, True, True]),
      (['A', 'B', 'B'], {'A', 'B', 'C'}, False, [True, True, True]),
      (['B', 'A', 'C'], {'A', 'B', 'C'}, False, [True, True, True]),
      (['AB::AB'], {_prevent_intern('AB::AB')}, False, [True]),
      ([], set(), True, []),
      ([], {'A', 'B'}, True, []),
      ([None, 'A'], {'A', 'B', 'C'}, True, [True, False]),
      (['A', 'B', 'C'], set(), True, [True, True, True]),
      (['A', 'B', 'C'], {'B'}, True, [True, False, True]),
      (['A', 'B', 'C'], frozenset(['A', 'B']), True, [False, False, True]),
      (['A', 'B', 'C'], {'A', 'B'}, True, [False, False, True]),
      (['A', 'B', 'C'], {'A', 'B', 'C'}, True, [False, False, False]),
      (['A', 'B', 'B'], {'A', 'B', 'C'}, True, [False, False, False]),
      (['B', 'A', 'C'], {'A', 'B', 'C'}, True, [False, False, False]),
      (['AB::AB'], {_prevent_intern('AB::AB')}, True, [False]),
  )
  def test_isin_string_object_array(
      self, array, test_elements, invert, exp_result
  ):
    result = membership.isin(
        np.array(array, dtype=object), test_elements, invert=invert
    )
    self.assertSequenceEqual(result.tolist(), exp_result)

  @parameterized.product(
      (
          dict(array=[], test_set=set(), isin_gt=[]),
          dict(array=[], test_set={'A', 'B'}, isin_gt=[]),
          dict(
              array=[None, 'A'],
              test_set={'A', 'B', 'C'},
              isin_gt=[False, True],
          ),
          dict(
              array=['A', 'B', 'C'],
              test_set=set(),
              isin_gt=[False, False, False],
          ),
          dict(
              array=['A', 'B', 'C'],
              test_set={'B'},
              isin_gt=[False, True, False],
          ),
          dict(
              array=['A', 'B', 'C'],
              test_set=frozenset(['A', 'B']),
              isin_gt=[True, True, False],
          ),
          dict(
              array=['A', 'B', 'C'],
              test_set={'A', 'B'},
              isin_gt=[True, True, False],
          ),
          dict(
              array=['A', 'B', 'C'],
              test_set={'A', 'B', 'C'},
              isin_gt=[True, True, True],
          ),
          dict(
              array=['A', 'B', 'B'],
              test_set={'A', 'B', 'C'},
              isin_gt=[True, True, True],
          ),
          dict(
              array=['B', 'A', 'C'],
              test_set={'A', 'B', 'C'},
              isin_gt=[True, True, True],
          ),
          dict(
              array=np.full((9, 9), 'A')[3:6, 3:6],
              test_set={'A', 'B'},
              isin_gt=[[True] * 3] * 3,
          ),
          dict(
              array=['AB::AB'],
              test_set={_prevent_intern('AB::AB')},
              isin_gt=[True],
          ),
      ),
      dtype=(object, str, np.dtype('T'), np.dtypes.StringDType(na_object=None)),
      invert=(False, True),
  )
  def test_isin_string_array(self, array, test_set, invert, isin_gt, dtype):
    exp_results = np.array(isin_gt, dtype=bool)
    if invert:
      exp_results = ~exp_results
    result = membership.isin(
        np.array(array, dtype=dtype), test_set, invert=invert
    )
    np.testing.assert_array_equal(result.tolist(), exp_results)

  @parameterized.product(
      targets=({'None', 'A'}, {None, 'A'}),
      invert=(False, True),
  )
  def test_isin_string_array_null_string_na_object(self, targets, invert):
    arr = np.array([None, 'A'], dtype=np.dtypes.StringDType(na_object=None))
    result = membership.isin(arr, targets, invert=invert)
    expected = np.array([invert, not invert], dtype=bool)
    np.testing.assert_array_equal(result, expected, strict=True)

  @parameterized.product(
      targets_and_expected=(
          ({'N/A'}, [True, False]),
          ({'A'}, [False, True]),
          ({'N/A', 'A'}, [True, True]),
          ({'B'}, [False, False]),
      ),
      invert=(False, True),
  )
  def test_isin_string_array_string_na_object(
      self, targets_and_expected, invert
  ):
    targets, base_expected = targets_and_expected
    arr = np.array(['N/A', 'A'], dtype=np.dtypes.StringDType(na_object='N/A'))
    result = membership.isin(arr, targets, invert=invert)
    expected = np.array(base_expected, dtype=bool)
    if invert:
      expected = ~expected
    np.testing.assert_array_equal(result, expected, strict=True)

  @parameterized.parameters(object, str, np.dtype('T'))
  def test_isin_string_array_multi_dim(self, dtype):
    result = membership.isin(
        np.array([['A', 'B', 'Z'], ['C', 'D', 'Z']], dtype=dtype),
        {'A', 'D'},
    )
    self.assertSequenceEqual(
        result.tolist(), [[True, False, False], [False, True, False]]
    )

  @parameterized.named_parameters(
      dict(testcase_name='scalar', test_elements=1),
      dict(testcase_name='numpy_scalar', test_elements=np.int64(1)),
      dict(testcase_name='list', test_elements=[1]),
      dict(testcase_name='tuple', test_elements=(1,)),
  )
  def test_isin_test_elements_invalid_type(self, test_elements):
    with self.assertRaises(TypeError):
      membership.isin(np.array([1]), test_elements=test_elements)

  def test_isin_float64_dtype_array(self):
    result = membership.isin(
        np.array([1.0, 2.0, 3.0], dtype=np.float64), {1, 2}
    )
    self.assertSequenceEqual(result.tolist(), [True, True, False])

  @parameterized.product(
      layout_and_f_contig=(
          ('F', True),
          ('transpose', True),
          ('strided', False),
          ('neg_strided', False),
      ),
      dtype=(np.dtypes.StringDType(), object, str),
      test_elements=({'A', 'E', 'F'}, set()),
      invert=(False, True),
  )
  def test_isin_string_array_non_c_contiguous_layouts(
      self, layout_and_f_contig, dtype, test_elements, invert
  ):
    layout, expect_f_contiguous = layout_and_f_contig
    base = np.array(
        [['A', 'B', 'C', 'D'], ['E', 'F', 'G', 'H'], ['I', 'J', 'K', 'L']],
        dtype=dtype,
    )
    if layout == 'F':
      arr = np.asfortranarray(base)
    elif layout == 'transpose':
      arr = base.T
    elif layout == 'strided':
      arr = base[::2, ::2]
    elif layout == 'neg_strided':
      arr = base[::-1, ::-1]
    else:
      raise ValueError(f'Unknown layout: {layout}')

    result = membership.isin(arr, test_elements, invert=invert)
    expected = np.isin(arr.astype(object), list(test_elements), invert=invert)
    np.testing.assert_array_equal(result, expected)
    if dtype == np.dtypes.StringDType():
      self.assertEqual(result.flags.f_contiguous, expect_f_contiguous)

  @parameterized.named_parameters(
      dict(
          testcase_name='non_empty',
          arr_str=['A', 'B', 'C'],
          test_str={'A', 'B'},
      ),
      dict(
          testcase_name='empty_test_elements',
          arr_str=['A', 'B', 'C'],
          test_str=set(),
      ),
      dict(
          testcase_name='zero_size',
          arr_str=[],
          test_str={'A', 'B'},
      ),
  )
  def test_isin_string_array_refcounts(self, arr_str, test_str):
    # The returned array must not carry an extra reference from the NumPy
    # iterator operand it is borrowed from as this would leak every array
    # returned by isin(). No need to check elements as they are packed
    # strings in an arena with no per-element refcount.
    str_array = np.array(arr_str, dtype=np.dtypes.StringDType())
    res_str = membership.isin(str_array, test_str)
    # Compare against a freshly allocated array rather than hardcoding a
    # constant: two fresh locals always agree, without depending on the
    # interpreter's borrowed-reference semantics.
    fresh_numpy_array = np.empty(0, dtype=bool)
    self.assertEqual(
        sys.getrefcount(res_str),
        sys.getrefcount(fresh_numpy_array),
    )

  def test_isin_string_array_zero_size_non_contiguous(self):
    arr = np.empty((5, 0), dtype=np.dtypes.StringDType())[::2, :]
    res = membership.isin(arr, {'A', 'B'})
    self.assertEqual(res.shape, (3, 0))

  @parameterized.product(
      layout_fn=(
          np.asfortranarray,
          np.transpose,
          lambda x: x[::2, ::2],
          lambda x: x[::-1, ::-1],
      ),
      dtype=(*_INT_DTYPES, np.longlong),
      test_elements=({1, 5, 6}, set()),
      invert=(False, True),
  )
  def test_isin_int_array_non_c_contiguous_layouts(
      self, layout_fn, dtype, test_elements, invert
  ):
    arr = layout_fn(np.arange(1, 13, dtype=dtype).reshape(3, 4))
    self.assertFalse(arr.flags.c_contiguous)
    np.testing.assert_array_equal(
        membership.isin(arr, test_elements, invert=invert),
        np.isin(arr, list(test_elements), invert=invert),
        strict=True,
    )

  @parameterized.parameters(*_INT_DTYPES, np.longlong)
  def test_isin_non_contiguous_int_rejects_incompatible_test_elements(
      self, dtype
  ):
    with self.assertRaisesRegex(RuntimeError, 'Unable to cast'):
      membership.isin(np.arange(4, dtype=dtype)[::2], {'A'})

  @parameterized.parameters(*_INT_DTYPES)
  def test_isin_int_out_of_range_test_elements(self, dtype):
    arr = np.array([1, 2, 3], dtype=dtype)
    np.testing.assert_array_equal(
        membership.isin(arr, {-1, 2, 10**12, 10**30}),
        [False, True, False],
        strict=True,
    )

  def test_isin_uint64_above_int64_max(self):
    arr = np.array([1, 2**63 + 5, 3], dtype=np.uint64)
    np.testing.assert_array_equal(
        membership.isin(arr, {2**63 + 5, 99}),
        [False, True, False],
        strict=True,
    )

  def test_isin_object_array_mutated_during_eq(self):
    arr = np.empty(2, dtype=object)

    class MutatingKey:

      def __hash__(self):
        return 1

      def __eq__(self, other):
        arr[0] = None
        return other == 1

    arr[0] = MutatingKey()
    arr[1] = 1
    np.testing.assert_array_equal(
        membership.isin(arr, {1}),
        np.array([True, True], dtype=bool),
        strict=True,
    )

  @parameterized.product(
      case=(
          ([], [1, 2]),
          ([1, 2, 3], []),
          ([], []),
          ([0, 1, 2], [1]),
          ([-1, 0, 1, 2, 3], [0, 0, 1, 2, 2]),
          ([-1, 0, 50, 51, 998, 999], list(range(0, 1000, 2))),
          ([10, 12, 15, 25], [25, 10, 15, 15]),
          ([-(10**9), 0, 42, 10**9], [10**9, -(10**9), 42, 42]),
          (
              [np.iinfo(np.int64).min, -1, 0, 1, np.iinfo(np.int64).max],
              [np.iinfo(np.int64).max, np.iinfo(np.int64).min],
          ),
      ),
      invert=(False, True),
  )
  def test_isin_int_array_test_elements(self, case, invert):
    arr = np.array(case[0], dtype=np.int64)
    test_arr = np.array(case[1], dtype=np.int64)
    np.testing.assert_array_equal(
        membership.isin(arr, test_arr, invert=invert),
        np.isin(arr, test_arr, invert=invert),
        strict=True,
    )

  @parameterized.product(
      array_dtype=(np.longlong, np.int64, np.int32, np.int8),
      test_dtype=(np.longlong, np.int64, np.int32, np.int8),
      strided=(False, True),
      invert=(False, True),
  )
  def test_isin_int_array_cross_dtypes_and_layouts(
      self, array_dtype, test_dtype, strided, invert
  ):
    vals = (
        [-(2**40), -5, 0, 3, 7, 2**40]
        if array_dtype in (np.int64, np.longlong)
        else [-128, -5, 0, 3, 7, 127]
    )
    arr = np.array(vals, dtype=array_dtype).reshape(2, 3)
    test_arr = np.array([-128, 3, 3, 127], dtype=test_dtype).reshape(2, 2)
    if strided:
      arr = arr[::-1, ::2]
      test_arr = test_arr.T[:, ::-1]
    np.testing.assert_array_equal(
        membership.isin(arr, test_arr, invert=invert),
        np.isin(arr, test_arr, invert=invert),
        strict=True,
    )

  @parameterized.named_parameters(
      dict(
          testcase_name='int_array_str_test',
          arr=np.array([1, 2]),
          test_arr=np.array(['1']),
      ),
      dict(
          testcase_name='int_array_empty_untyped_test',
          arr=np.array([1, 2], dtype=np.int64),
          test_arr=np.array([]),
      ),
      dict(
          testcase_name='int_array_uint_test',
          arr=np.array([1, 2], dtype=np.int64),
          test_arr=np.array([1], dtype=np.uint64),
      ),
      dict(
          testcase_name='int_array_float_test',
          arr=np.array([1, 2], dtype=np.int64),
          test_arr=np.array([1.0], dtype=np.float64),
      ),
      dict(
          testcase_name='int_array_bool_test',
          arr=np.array([0, 1]),
          test_arr=np.array([True]),
      ),
      dict(
          testcase_name='empty_int_array_str_test',
          arr=np.array([], dtype=np.int64),
          test_arr=np.array(['1']),
      ),
      dict(
          testcase_name='string_dtype_array_test',
          arr=np.array(['A'], dtype=np.dtype('T')),
          test_arr=np.array(['A'], dtype=np.dtype('T')),
      ),
      dict(
          testcase_name='string_dtype_empty_untyped_test',
          arr=np.array(['A'], dtype=np.dtype('T')),
          test_arr=np.array([]),
      ),
      dict(
          testcase_name='object_array_test',
          arr=np.array(['A'], dtype=object),
          test_arr=np.array(['A'], dtype=object),
      ),
      dict(
          testcase_name='float_array_int_test',
          arr=np.array([1.0], dtype=np.float64),
          test_arr=np.array([1], dtype=np.int64),
      ),
      dict(
          testcase_name='list_array_int_test',
          arr=[1, 2],
          test_arr=np.array([1], dtype=np.int64),
      ),
  )
  def test_isin_array_incompatible_dtypes(self, arr, test_arr):
    with self.assertRaises(TypeError):
      membership.isin(arr, test_arr)

  @parameterized.parameters(np.int64, np.int32, object, np.dtype('T'))
  def test_isin_concurrent_threads(self, dtype):
    if dtype in (np.int64, np.int32):
      arr = np.arange(300, dtype=dtype)
      targets = {10, 50, 150, 299}
      expected = np.isin(arr, np.array(sorted(targets), dtype=dtype))
    else:
      arr = np.array([f'A{i % 20}' for i in range(300)], dtype=dtype)
      targets = {'A1', 'A7', 'A19'}
      expected = np.isin(arr.astype(object), list(targets))
    with futures.ThreadPoolExecutor(max_workers=8) as pool:
      results = list(
          pool.map(lambda _: membership.isin(arr, targets), range(32))
      )
    for result in results:
      np.testing.assert_array_equal(result, expected, strict=True)


if __name__ == '__main__':
  absltest.main()
