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

  @parameterized.parameters(False, True)
  def test_isin_string_array_null_string_na_object(self, invert):
    arr = np.array([None, 'A'], dtype=np.dtypes.StringDType(na_object=None))
    result = membership.isin(arr, {'None', 'A'}, invert=invert)
    expected = [invert, not invert]
    self.assertSequenceEqual(result.tolist(), expected)

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
      dict(testcase_name='list', test_elements=[1]),
      dict(testcase_name='numpy_array', test_elements=np.array([1])),
  )
  def test_isin_test_elements_not_set(self, test_elements):
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


if __name__ == '__main__':
  absltest.main()
