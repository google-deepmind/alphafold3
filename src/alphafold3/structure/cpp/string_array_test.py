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
from alphafold3.cpp import string_array
import numpy as np

_TEST_SUCCESS_EXAMPLES = [
    dict(
        values=list('Hello'),
        mapping=dict(H='W', l='o'),
        expected=list('Weooo'),
    ),
    dict(
        values=list('Hello'),
        mapping=dict(H='W', l='o'),
        default_value='a',
        expected=list('Waooa'),
    ),
    dict(
        values=list('AllMissing'),
        mapping=dict(H='W', z='o'),
        expected=list('AllMissing'),
    ),
    dict(
        values=[1, dict(a='b'), 3],
        mapping={3: 30},
        expected=[1, dict(a='b'), 30],
    ),
    dict(
        values=[1, dict(a='b'), 3],
        mapping={3: 30},
        default_value='a',
        expected=np.array(['a', 'a', 30], dtype=object),
    ),
    dict(
        values=[],
        mapping={'Empty': 'Values'},
        expected=[],
    ),
    dict(
        values=[[['A'], ['B']], [['A'], ['B']]],
        mapping={'B': 'C'},
        expected=[[['A'], ['C']], [['A'], ['C']]],
    ),
    dict(
        values=list('Empty Dict'),
        mapping={},
        expected=list('Empty Dict'),
    ),
    dict(
        values=list('Empty Dict Has Default'),
        mapping={},
        default_value='a',
        expected=['a'] * 22,
    ),
    dict(
        values=list('None default'),
        mapping={},
        default_value=None,
        expected=[None] * 12,
    ),
]


def _prevent_intern(s: str) -> str:
  """Prevents Python from interning the given string."""
  return ''.join(s)


class StringArrayTest(parameterized.TestCase):

  @parameterized.parameters(_TEST_SUCCESS_EXAMPLES)
  def test_copy(self, values, mapping, expected, **kwargs):
    start_values = np.asarray(values, dtype=object)
    actual_values = string_array.remap(start_values, mapping, **kwargs)
    np.testing.assert_array_equal(
        actual_values, np.array(expected, dtype=object), strict=True
    )
    # Start left unchanged.
    np.testing.assert_array_equal(
        start_values, np.array(values, dtype=object), strict=True
    )

  @parameterized.parameters(_TEST_SUCCESS_EXAMPLES)
  def test_inplace(self, values, mapping, expected, **kwargs):
    start_values = np.asarray(values, dtype=object)
    string_array.remap(start_values, mapping, inplace=True, **kwargs)
    np.testing.assert_array_equal(
        start_values, np.array(expected, dtype=object), strict=True
    )

  def test_not_dict(self):
    start_values = np.array(['asd', 'fgh', 'jkl'], dtype=object)
    with self.assertRaises(TypeError):
      string_array.remap(start_values, ('1', '2'), inplace=True)

  def test_not_numpy_array(self):
    start_values = ['asd', 'fgh', 'jkl']
    with self.assertRaises(TypeError):
      string_array.remap(start_values, {}, inplace=True)

  def test_not_numpy_array_objects(self):
    start_values = np.array(['asd', 'fgh', 'jkl'])
    with self.assertRaises(TypeError):
      string_array.remap(start_values, {}, inplace=True)

  @parameterized.parameters(
      (np.array([1.07, 2.12345], dtype=np.float32), 0, ['1', '2']),
      (np.array([1.07, 2.12345], dtype=np.float32), 1, ['1.1', '2.1']),
      (np.array([1.07, 2.12345], dtype=np.float32), 2, ['1.07', '2.12']),
      (np.array([1.07, 2.12345], dtype=np.float32), 3, ['1.070', '2.123']),
      (np.array([1.07, 2.12345], dtype=np.float32), 4, ['1.0700', '2.1235']),
      (
          np.array([1.07, 2.12345, -55.5625, 55.5625], dtype=np.float32),
          3,
          ['1.070', '2.123', '%.3f' % -55.5625, '%.3f' % 55.5625],
      ),
  )
  def test_format_float_array(self, values, num_decimals, formatted_array):
    self.assertEqual(
        string_array.format_float_array(values, num_decimals),
        formatted_array,
    )

  @parameterized.product(num_decimals=list(range(1, 4)))
  def test_format_float_array_coords(self, num_decimals):
    rng = np.random.default_rng(42)
    values = rng.uniform(-100, 100, size=40_000).astype(np.float32)
    actuals = string_array.format_float_array(
        values, num_decimal_places=num_decimals
    )
    for i, (value, actual) in enumerate(zip(values, actuals)):
      self.assertEqual(
          actual,
          (expected := '%.*f' % (num_decimals, value)),
          msg=f'Mismatch at index {i}: {actual} != {expected}',
      )

  @parameterized.parameters(
      ([], set(), False, []),
      ([], set(['A', 'B']), False, []),
      ([None, 'A'], set(['A', 'B', 'C']), False, [False, True]),
      (['A', 'B', 'C'], set([]), False, [False, False, False]),
      (['A', 'B', 'C'], set(['B']), False, [False, True, False]),
      (['A', 'B', 'C'], set(['A', 'B']), False, [True, True, False]),
      (['A', 'B', 'C'], set(['A', 'B', 'C']), False, [True, True, True]),
      (['A', 'B', 'B'], set(['A', 'B', 'C']), False, [True, True, True]),
      (['B', 'A', 'C'], set(['A', 'B', 'C']), False, [True, True, True]),
      (['AB::AB'], set([_prevent_intern('AB::AB')]), False, [True]),
      ([], set(), True, []),
      ([], set(['A', 'B']), True, []),
      ([None, 'A'], set(['A', 'B', 'C']), True, [True, False]),
      (['A', 'B', 'C'], set([]), True, [True, True, True]),
      (['A', 'B', 'C'], set(['B']), True, [True, False, True]),
      (['A', 'B', 'C'], set(['A', 'B']), True, [False, False, True]),
      (['A', 'B', 'C'], set(['A', 'B', 'C']), True, [False, False, False]),
      (['A', 'B', 'B'], set(['A', 'B', 'C']), True, [False, False, False]),
      (['B', 'A', 'C'], set(['A', 'B', 'C']), True, [False, False, False]),
      (['AB::AB'], set([_prevent_intern('AB::AB')]), True, [False]),
  )
  def test_isin(self, array, test_elements, invert, exp_result):
    result = string_array.isin(
        np.array(array, dtype=object), test_elements, invert=invert
    )
    self.assertSequenceEqual(result.tolist(), exp_result)

  def test_isin_multi_dim(self):
    result = string_array.isin(
        np.array([['A', 'B', 'Z'], ['C', 'D', 'Z']], dtype=object),
        set(['A', 'D']),
    )
    self.assertSequenceEqual(
        result.tolist(), [[True, False, False], [False, True, False]]
    )

  def test_remap_multiple_returns_expected(self):
    arrays = (
        np.array(['a', 'b', 'c'], dtype=object),
        np.array(['d', 'e', 'f'], dtype=object),
    )
    mapping = {('a', 'd'): 1, ('b', 'e'): 2, ('c', 'f'): 3}
    np.testing.assert_array_equal(
        string_array.remap_multiple(arrays, mapping),
        np.array([1, 2, 3], dtype=np.int64),
        strict=True,
    )

  def test_remap_multiple_empty(self):
    np.testing.assert_array_equal(
        string_array.remap_multiple(arrays=[], mapping={}),
        np.array([], dtype=np.int64),
        strict=True,
    )

  def test_remap_multiple_raises_error_on_inconsistent_input(self):
    arrays = (
        np.array(['a', 'b', 'c'], dtype=object),
        np.array(['d', 'e'], dtype=object),
    )
    mapping = {('a', 'd'): 1, ('b', 'e'): 2, ('c', 'f'): 3}

    with self.assertRaises(ValueError):
      string_array.remap_multiple(arrays, mapping)

  def test_remap_multiple_raises_error_on_non_int_mapping_value(self):
    arrays = (
        np.array(['a', 'b', 'c'], dtype=object),
        np.array(['d', 'e', 'f'], dtype=object),
    )
    mapping = {('a', 'd'): 1, ('b', 'e'): 2, ('c', 'f'): 3.6}

    with self.assertRaises(TypeError):
      string_array.remap_multiple(arrays, mapping)

  def test_remap_multiple_raises_error_on_non_existing_key(self):
    arrays = (
        np.array(['a', 'b', 'c'], dtype=object),
        np.array(['d', 'e', 'f'], dtype=object),
    )
    mapping = {('a', 'd'): 1, ('b', 'e'): 2}

    with self.assertRaises(KeyError):
      string_array.remap_multiple(arrays, mapping)


if __name__ == '__main__':
  absltest.main()
