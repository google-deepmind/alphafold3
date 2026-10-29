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

from collections.abc import Mapping, Sequence
import pickle

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.common import resources
from alphafold3.common.testing import data
from alphafold3.cpp import cif_dict
import numpy as np


_EXP_CIF_STRING_MULTILINE_FIRST = """data_test
#
loop_
_mmCIF.a
_mmCIF.z
;L1
L2
;
z1 
a2 z2 
#
"""

_EXP_CIF_STRING_MULTILINE_MIDDLE = """data_test
#
loop_
_mmCIF.a
_mmCIF.k
_mmCIF.z
a1 
;L1
L2
;
z1 
a2 k2 z2 
#
"""

_EXP_CIF_STRING_MULTILINE_LAST = """data_test
#
loop_
_mmCIF.a
_mmCIF.z
a1 
;L1
L2
;
a2 z2 
#
"""

_EXP_CIF_STRING_QUOTES_IN_VALUES = """data_tst
#
loop_
_k.a
'A "A" A' 
"B 'B' B" 
;C "C" 'C' C
;

#
"""

_TEST_DATA = data.Data(
    resources.ROOT / 'test_data/miniature_databases/pdb_mmcif'
)


class CifDictTest(parameterized.TestCase):

  def setUp(self):
    cif_dict_data = _TEST_DATA.load('2pv7.cif')
    self.cif_2pv7 = cif_dict.from_string(cif_dict_data)
    super().setUp()

  def assertSameMapping(
      self,
      lhs: Mapping[str, Sequence[str] | str] | cif_dict.CifDict,
      rhs: Mapping[str, Sequence[str] | str] | cif_dict.CifDict,
  ):
    for k, lhs_v in lhs.items():
      self.assertIn(k, rhs)
      self.assertEqual(lhs_v, rhs[k])
    for k in rhs:
      self.assertIn(k, lhs)

  def test_tokenize_comment(self):
    tokenized = cif_dict.tokenize('# Comment')
    self.assertListEqual(tokenized, [])

  def test_tokenize_multiline(self):
    tokenized = cif_dict.tokenize(';A\nB \n C\n;')
    self.assertListEqual(tokenized, ['A\nB\n C'])

  def test_tokenize_multiline_containing_empty_line(self):
    tokenized = cif_dict.tokenize(';A\n\n;')
    self.assertListEqual(tokenized, ['A\n'])

  def test_split_line(self):
    split_line = cif_dict.split_line('A "B" ""C" "D E" \'F\' # Comment')
    self.assertListEqual(split_line, ['A', 'B', '"C', 'D E', 'F'])

  def test_split_line_quote_in_word(self):
    tokenized = cif_dict.split_line('T1" T2\' T"3 T4')
    self.assertListEqual(tokenized, ['T1"', "T2'", 'T"3', 'T4'])

  def test_split_line_quote_in_quoted_word(self):
    tokenized = cif_dict.split_line('""A" \'\'B\' "C"" \'D\'\' ""E "')
    self.assertListEqual(tokenized, ['"A', "'B", 'C"', "D'", '"E '])

  def test_split_line_quote_at_the_end(self):
    with self.assertRaisesRegex(
        ValueError, r'Line ended with quote open: T "T2'
    ):
      cif_dict.split_line('T "T2')

  @parameterized.named_parameters(
      (
          'first_column',
          {
              'data_': ['test'],
              '_mmCIF.a': ['L1\nL2', 'a2'],
              '_mmCIF.z': ['z1', 'z2'],
          },
          _EXP_CIF_STRING_MULTILINE_FIRST,
      ),
      (
          'middle_column',
          {
              'data_': ['test'],
              '_mmCIF.a': ['a1', 'a2'],
              '_mmCIF.k': ['L1\nL2', 'k2'],
              '_mmCIF.z': ['z1', 'z2'],
          },
          _EXP_CIF_STRING_MULTILINE_MIDDLE,
      ),
      (
          'last_column',
          {
              'data_': ['test'],
              '_mmCIF.a': ['a1', 'a2'],
              '_mmCIF.z': ['L1\nL2', 'z2'],
          },
          _EXP_CIF_STRING_MULTILINE_LAST,
      ),
  )
  def test_to_string_loop_with_multiline(self, cif_in, exp_cif):
    cif = cif_dict.CifDict(cif_in)
    self.assertEqual(cif.to_string(), exp_cif)

  @parameterized.parameters(['', '  ', '\t \n  \t'])
  def test_from_string_empty(self, empty_string: str):
    with self.assertRaisesRegex(ValueError, 'must not be empty'):
      cif_dict.from_string(empty_string)

  def test_to_string_quotes_in_value(self):
    raw = {'data_': ['tst'], '_k.a': ['A "A" A', "B 'B' B", 'C "C" \'C\' C']}
    cif = cif_dict.CifDict(raw)
    roundtrip_cif = cif_dict.from_string(cif.to_string())
    self.assertSameMapping(roundtrip_cif, raw)
    self.assertEqual(cif.to_string(), _EXP_CIF_STRING_QUOTES_IN_VALUES)

  def test_from_string_just_name(self):
    cif = cif_dict.from_string('data_name')
    self.assertEqual(cif.to_string(), 'data_name\n#\n')
    self.assertSameMapping(cif, dict(data_=['name']))

  @parameterized.named_parameters(
      ('double_quote_key', '"_key.a"', ['A', 'B']),
      ('single_quote_key', "'_key.a'", ['A', 'B']),
      ('double_quote_key_middle', '_"key".a', ['A', 'B']),
      ('single_quote_key_middle', "_'key'.a", ['A', 'B']),
      ('double_quote_value', '"_key.a"', ['"A"', 'B']),
      ('single_quote_value', "'_key.a'", ["'A'", 'B']),
      ('double_quote_value_middle', '_key.a', ['A"A"A', 'B']),
      ('single_quote_value_middle', '_key.a', ["A'A'A", 'B']),
  )
  def test_from_string_paired_quotes(self, key, values):
    cif_values = '\n'.join(values)
    cif = cif_dict.from_string(f'data_TST\nloop_\n{key}\n{cif_values}')
    expected = {
        key.strip('\'"'): [v.strip('\'"') for v in values],
        'data_': ['TST'],
    }
    self.assertSameMapping(cif, expected)

  @parameterized.named_parameters(
      ('double_quote_key_middle', '_key."a', ['A', 'B']),
      ('double_quote_key_end', '_key.a"', ['A', 'B']),
      ('single_quote_key_middle', "_key.'a", ['A', 'B']),
      ('single_quote_key_end', "_key.a'", ['A', 'B']),
      ('double_quote_value_middle', '_key.a', ['A"A', 'B']),
      ('double_quote_value_end', '_key.a', ['A"', 'B']),
      ('single_quote_value_middle', '_key.a', ["A'A", 'B']),
      ('single_quote_value_end', '_key.a', ["A'", 'B']),
  )
  def test_from_string_unpaired_quotes(self, key, values):
    cif_values = '\n'.join(values)
    cif = cif_dict.from_string(f'data_TST\nloop_\n{key}\n{cif_values}')
    expected = {key: values, 'data_': ['TST']}
    self.assertSameMapping(cif, expected)

  @parameterized.named_parameters(
      ('at_end', 'data_TST\nloop_\n_k.a\n_k.b\na11\nb11\na11\n'),
      ('before_loop', 'data_TST\nloop_\n_k.a\n_k.b\na11\nloop_\n_k.c\nc'),
      ('before_key', 'data_TST\nloop_\n_k.a\n_k.b\na11\n_k.c c11'),
  )
  def test_from_string_bad_loop_columns(self, cif):
    with self.assertRaisesRegex(ValueError, 'The number of values'):
      cif_dict.from_string(cif)

  @parameterized.parameters(
      ('data_TST\n_k v1 v2', 'v2'), ('data_TST\n_k v1 v2 v3', 'v2')
  )
  def test_from_string_non_loop_key_multiple_values(self, cif, bad_key):
    with self.assertRaisesRegex(ValueError, f"Key '{bad_key}' does not"):
      cif_dict.from_string(cif)

  def test_crlf_newlines(self):
    cif = cif_dict.from_string('data_CRLF\r\n_k1 v1\r\n_k2 v2\r\n')
    self.assertEqual(cif['data_'], ['CRLF'])
    self.assertEqual(cif['_k1'], ['v1'])
    self.assertEqual(cif['_k2'], ['v2'])

  def test_crlf_newlines_multiline(self):
    cif = cif_dict.from_string('data_CRLF\r\n_k\r\n;L1\r\nL2\r\n;')
    self.assertEqual(cif['data_'], ['CRLF'])
    self.assertEqual(cif['_k'], ['L1\nL2'])

  def test_value_length(self):
    cif = cif_dict.CifDict({
        'data_': ['name'],
        '_count_2': ['1', '2'],
        '_count_4': ['1', '2', '3', '4'],
    })
    self.assertEqual(cif.value_length('_count_2'), 2)
    self.assertEqual(cif.value_length('_count_4'), 4)
    self.assertEqual(cif.value_length('_count_missing'), 0)

  @parameterized.product(
      cif_in=(
          {'data_': ['name']},
          {'data_': ['name'], '_counts': ['1', '2', '3', '4']},
          {
              'data_': ['name'],
              '_counts.x': ['x11\nx12', 'x2', 'x\t3', 'x4'],
              '_counts.y': ['y1', 'y2', 'y3', 'y4'],
              '_counts.z': ['_loop', 'z2', 'z 3', 'z4'],
          },
      )
  )
  def test_round_trip(self, cif_in):
    convert_to_cif_dict = cif_dict.CifDict(cif_in)
    convert_to_string = convert_to_cif_dict.to_string()
    cif_out = cif_dict.from_string(convert_to_string)
    self.assertSameMapping(cif_out, cif_in)

  def test_keys(self):
    cif = cif_dict.CifDict({
        'data_': ['name'],
        '_counts.x': ['x1', 'x2', 'x3', 'x4'],
        '_counts.y': ['y1', 'y2', 'y3', 'y4'],
        '_counts.z': ['z1', 'z2', 'z3', 'z4'],
    })
    self.assertSameElements(
        cif.keys(), ['data_', '_counts.x', '_counts.y', '_counts.z']
    )
    self.assertSameElements(
        cif, ['data_', '_counts.x', '_counts.y', '_counts.z']
    )

  def test_values(self):
    cif = cif_dict.CifDict({
        'data_': ['name'],
        '_counts.x': ['x1', 'x2'],
        '_counts.y': ['y1', 'y2'],
        '_counts.z': ['z1', 'z2'],
    })
    self.assertSameElements(
        cif.values(), [['name'], ['x1', 'x2'], ['y1', 'y2'], ['z1', 'z2']]
    )

  def test_items(self):
    cif = cif_dict.CifDict({
        'data_': ['name'],
        '_counts.x': ['x1', 'x2'],
        '_counts.y': ['y1', 'y2'],
        '_counts.z': ['z1', 'z2'],
    })
    self.assertSameElements(
        cif.items(),
        [
            ('data_', ['name']),
            ('_counts.x', ['x1', 'x2']),
            ('_counts.y', ['y1', 'y2']),
            ('_counts.z', ['z1', 'z2']),
        ],
    )

  def test_to_dict(self):
    raw = {
        'data_': ['name'],
        '_counts.x': ['x1', 'x2'],
        '_counts.y': ['y1', 'y2'],
        '_counts.z': ['z1', 'z2'],
    }
    cif = cif_dict.CifDict(raw)
    result = cif.to_dict()
    self.assertIsInstance(result, dict)
    self.assertDictEqual(result, raw)

  def test_len(self):
    cif = cif_dict.CifDict({
        'data_': ['name'],
        '_counts.x': ['x1', 'x2'],
        '_counts.y': ['y1', 'y2'],
        '_counts.z': ['z1', 'z2'],
    })
    self.assertLen(cif, 4)

  def test_contains(self):
    cif = cif_dict.CifDict({
        'data_': ['name'],
        '_counts.x': ['x1', 'x2'],
    })
    self.assertIn('data_', cif)
    self.assertIn('_counts.x', cif)
    self.assertIn(b'_counts.x', cif)
    self.assertNotIn('_counts.y', cif)
    with self.assertRaises(TypeError):
      _ = ['name'] in cif

  def test_extract_loop_as_list(self):
    cif = cif_dict.CifDict({
        'data_': ['name'],
        '_counts.x': ['x1', 'x2'],
        '_counts.y': ['y1', 'y2'],
        '_counts.z': ['z1', 'z2'],
    })
    result = cif.extract_loop_as_list('_counts.')
    self.assertLen(result, 2)
    self.assertSameMapping(
        result[0], {'_counts.x': 'x1', '_counts.y': 'y1', '_counts.z': 'z1'}
    )
    self.assertSameMapping(
        result[1], {'_counts.x': 'x2', '_counts.y': 'y2', '_counts.z': 'z2'}
    )

  def test_extract_loop_as_dict(self):
    cif = cif_dict.CifDict({
        'data_': ['name'],
        '_a.ix': ['r1', 'r2', 'r3'],
        '_a.1': ['a.1.1', 'a.1.2', 'a.1.3'],
        '_a.2': ['a.2.1', 'a.2.2', 'a.2.3'],
    })
    result = cif.extract_loop_as_dict('_a.', index='_a.ix')
    self.assertDictEqual(
        result,
        {
            'r1': {
                '_a.ix': 'r1',
                '_a.1': 'a.1.1',
                '_a.2': 'a.2.1',
            },
            'r2': {
                '_a.ix': 'r2',
                '_a.1': 'a.1.2',
                '_a.2': 'a.2.2',
            },
            'r3': {
                '_a.ix': 'r3',
                '_a.1': 'a.1.3',
                '_a.2': 'a.2.3',
            },
        },
    )

  def test_pickle(self):
    cif = cif_dict.CifDict({
        'data_': ['name'],
        '_counts.x': ['x1', 'x2', 'x3', 'x4'],
        '_counts.y': ['y1', 'y2', 'y3', 'y4'],
        '_counts.z': ['z1', 'z2', 'z3', 'z4'],
    })
    dumped = pickle.dumps(cif, protocol=4)
    loaded = pickle.loads(dumped)
    self.assertSameMapping(cif, loaded)

  @parameterized.parameters(
      np.int8,
      np.int16,
      np.int32,
      np.int64,
  )
  def test_get_array_as_signed_ints(self, dtype):
    values = ['-1', '1', '2', '3', '4']
    cif = cif_dict.CifDict({
        'data_': ['name'],
        '_counts.x': values,
    })
    np.testing.assert_array_equal(
        cif.get_array('_counts.x', dtype=dtype),
        np.array(values, dtype=object).astype(dtype),
        strict=True,
    )

  @parameterized.parameters(
      np.uint8,
      np.uint16,
      np.uint32,
      np.uint64,
  )
  def test_get_array_as_unsigned_ints(self, dtype):
    values = ['0', '1', '2', '3', '4']
    cif = cif_dict.CifDict({
        'data_': ['name'],
        '_counts.x': values,
    })
    np.testing.assert_array_equal(
        cif.get_array('_counts.x', dtype=dtype),
        np.array(values, dtype=object).astype(dtype),
        strict=True,
    )

  def test_get_array_as_bool(self):
    cif = cif_dict.CifDict({
        'data_': ['name'],
        '_counts.is_true': ['y', 'y', 'no', 'yes', 'no'],
    })
    np.testing.assert_array_equal(
        cif.get_array('_counts.is_true', dtype=bool),
        np.array([True, True, False, True, False]),
        strict=True,
    )

  def test_get_array_as_strings(self):
    values = ['0', 'aaaa', 'bbb', 'aaaa', 'aaaa']
    cif = cif_dict.CifDict({
        'data_': ['name'],
        '_counts.name': values,
    })
    result = cif.get_array('_counts.name', dtype=object)
    np.testing.assert_array_equal(
        result, np.array(values, dtype=object), strict=True
    )
    with self.subTest('Reuses existing strings'):
      self.assertIs(result[1], result[3])
      self.assertIs(result[1], result[4])

  def test_get_array_default(self):
    values = ['0', 'aaaa', 'bbb', 'aaaa', 'aaaa']
    cif = cif_dict.CifDict({
        'data_': ['name'],
        '_counts.name': values,
    })
    result = cif.get_array('_counts.name')
    np.testing.assert_array_equal(
        result, np.array(values, dtype=object), strict=True
    )
    self.assertEqual(result.dtype, np.dtype(object))

  @parameterized.parameters(
      (np.int8, '-128', '128'),
      (np.int16, '-32768', '32768'),
      (np.int32, '-2147483648', '2147483648'),
      (np.uint8, '255', '-1'),
      (np.uint16, '65535', '-1'),
      (np.uint32, '4294967295', '-1'),
      (bool, 'yes', 'YeS!'),
      (np.float32, '0.5', '1e0.1'),
      (np.float64, '1.79e+308', '1e0.1'),
  )
  def test_get_array_invalid_value(self, dtype, valid_value, invalid_value):
    cif = cif_dict.CifDict({
        'data_': ['name'],
        '_counts.valid': [valid_value, valid_value, valid_value],
        '_counts.invalidfirst': [invalid_value, valid_value, valid_value],
        '_counts.invalidlast': [valid_value, valid_value, invalid_value],
        '_counts.invalidmiddle': [valid_value, invalid_value, valid_value],
    })

    np.testing.assert_array_equal(
        cif.get_array('_counts.valid', dtype=dtype),
        np.array([valid_value, valid_value, valid_value], dtype=object).astype(
            dtype
        ),
        strict=True,
    )

    with self.assertRaisesWithLiteralMatch(ValueError, invalid_value):
      cif.get_array('_counts.invalidfirst', dtype=dtype)

    with self.assertRaisesWithLiteralMatch(ValueError, invalid_value):
      cif.get_array('_counts.invalidlast', dtype=dtype)

    with self.assertRaisesWithLiteralMatch(ValueError, invalid_value):
      cif.get_array('_counts.invalidmiddle', dtype=dtype)

  @parameterized.parameters(np.float32, np.float64)
  def test_get_array_missing_values(self, dtype):
    cif = cif_dict.CifDict({
        'data_': ['name'],
        '_counts.sequence_id': ['1', '2', '.', '4'],
    })
    np.testing.assert_array_equal(
        np.isnan(cif.get_array('_counts.sequence_id', dtype=dtype)),
        [False, False, True, False],
        strict=True,
    )

  def test_get_array_bad_dtype(self):
    cif = cif_dict.CifDict({
        'data_': ['name'],
        '_counts.dates': ['2007-07-13', '2006-01-13', '2010-08-13'],
    })
    with self.assertRaisesWithLiteralMatch(
        TypeError, "Unsupported dtype 'datetime64'"
    ):
      cif.get_array('_counts.dates', dtype='datetime64')

    with self.assertRaisesWithLiteralMatch(TypeError, "Invalid dtype 'blah'"):
      cif.get_array('_counts.dates', dtype='blah')

  @parameterized.parameters(
      ('_atom_site.group_PDB', object),
      ('_atom_site.id', np.int32),
      ('_atom_site.type_symbol', object),
      ('_atom_site.label_atom_id', object),
      ('_atom_site.label_alt_id', object),
      ('_atom_site.label_comp_id', object),
      ('_atom_site.label_asym_id', object),
      ('_atom_site.label_entity_id', object),
      ('_atom_site.label_seq_id', object),
      ('_atom_site.pdbx_PDB_ins_code', object),
      ('_atom_site.Cartn_x', np.float32),
      ('_atom_site.Cartn_y', np.float32),
      ('_atom_site.Cartn_z', np.float32),
      ('_atom_site.occupancy', np.float32),
      ('_atom_site.B_iso_or_equiv', object),
      ('_atom_site.pdbx_formal_charge', object),
      ('_atom_site.auth_seq_id', object),
      ('_atom_site.auth_comp_id', object),
      ('_atom_site.auth_asym_id', object),
      ('_atom_site.auth_atom_id', object),
      ('_atom_site.pdbx_PDB_model_num', np.int32),
  )
  def test_get_array(self, key, dtype):
    array = self.cif_2pv7.get_array(key, dtype=dtype)
    self.assertEqual(array.shape, (4983,))
    self.assertEqual(array.dtype, np.dtype(dtype))

  @parameterized.product(
      key=['_atom_site.id'],
      dtype=[
          np.int16,
          np.uint16,
          np.int32,
          np.uint32,
          np.int64,
          np.uint64,
          np.float32,
          np.float64,
          object,
      ],
      slice_=[
          slice(None),
          slice(10),
          slice(5, 15),
          slice(-10, None),
          slice(5, 25, 2),
      ],
  )
  def test_get_array_slice(self, key, dtype, slice_):
    array1 = self.cif_2pv7.get_array(key, dtype=dtype)
    self.assertEqual(array1.shape, (4983,))
    array2 = self.cif_2pv7.get_array(key, dtype=dtype, gather=slice_)
    np.testing.assert_array_equal(array1[slice_], array2, strict=True)

  @parameterized.product(
      key=['_atom_site.id'],
      dtype=[
          np.int16,
          np.uint16,
          np.int32,
          np.uint32,
          np.int64,
          np.uint64,
          np.float32,
          np.float64,
          object,
      ],
      lookup_indexes=[
          20,
          [1, 2, 3, 50],
          [[10, 20, 30], [1, 2, 3]],
      ],
      lookup_dtype=[
          np.int16,
          np.uint16,
          np.int32,
          np.uint32,
          np.int64,
          np.uint64,
      ],
  )
  def test_get_array_gather(self, key, dtype, lookup_indexes, lookup_dtype):
    array1 = self.cif_2pv7.get_array(key, dtype=dtype)
    self.assertEqual(array1.shape, (4983,))
    look_up = np.array(lookup_indexes, lookup_dtype)
    array2 = self.cif_2pv7.get_array(key, dtype=dtype, gather=look_up)
    self.assertEqual(array2.shape, look_up.shape)
    np.testing.assert_array_equal(array1[look_up], array2)

  @parameterized.named_parameters(
      dict(
          testcase_name='missing_key',
          key='missing_key',
          dtype=object,
          gather=None,
          exception=KeyError,
          msg="'missing_key'",
      ),
      dict(
          testcase_name='invalid_gather',
          key='_atom_site.id',
          dtype=object,
          gather='Tomato',
          exception=TypeError,
          msg="Invalid gather 'Tomato'",
      ),
      dict(
          testcase_name='invalid_slice',
          key='_atom_site.id',
          dtype=object,
          gather=slice('Tomato'),
          exception=TypeError,
          msg=None,  # Not owned error message.
      ),
      dict(
          testcase_name='out_of_range_index',
          key='_atom_site.id',
          dtype=object,
          gather=np.array([2, -5, 10], dtype=np.int64),
          exception=IndexError,
          msg='index -5 is out of bounds for column with size 4983',
      ),
      dict(
          testcase_name='out_of_range_index_scalar',
          key='_atom_site.id',
          dtype=object,
          gather=np.array(-5, dtype=np.int64),
          exception=IndexError,
          msg='index -5 is out of bounds for column with size 4983',
      ),
      dict(
          testcase_name='invalid_index_dtype',
          key='_atom_site.id',
          dtype=object,
          gather=np.array('2002-10-27T04:30', dtype='datetime64'),
          exception=TypeError,
          msg='Unsupported NumPy array type.',
      ),
      dict(
          testcase_name='unsupported_dtype',
          key='_atom_site.id',
          dtype='datetime64',
          gather=np.array([2, 5, 10], dtype=np.int64),
          exception=TypeError,
          msg="Unsupported dtype 'datetime64'",
      ),
      dict(
          testcase_name='invalid_dtype',
          key='_atom_site.id',
          dtype=27,
          gather=np.array([2, 5, 10], dtype=np.int64),
          exception=TypeError,
          msg='Invalid dtype 27',
      ),
      dict(
          testcase_name='invalid_convert',
          key='_atom_site.id',
          dtype=np.uint8,
          gather=np.array([20, 1000, 30], dtype=np.int32),
          exception=ValueError,
          msg='1001',
      ),
      dict(
          testcase_name='invalid_convert_scalar',
          key='_atom_site.id',
          dtype=np.uint8,
          gather=np.array(1000, dtype=np.int32),
          exception=ValueError,
          msg='1001',
      ),
      dict(
          testcase_name='invalid_convert_multidim',
          key='_atom_site.id',
          dtype=np.uint8,
          gather=np.array([[[1, 1000]]], dtype=np.int32),
          exception=ValueError,
          msg='1001',
      ),
      dict(
          testcase_name='invalid_convert_slice',
          key='_atom_site.id',
          dtype=np.uint8,
          gather=slice(200, 500),
          exception=ValueError,
          msg='256',
      ),
  )
  def test_get_array_failure(self, key, dtype, gather, exception, msg):
    if msg is not None:
      with self.assertRaisesWithLiteralMatch(exception, msg):
        self.cif_2pv7.get_array(key, dtype=dtype, gather=gather)
    else:
      with self.assertRaises(exception):
        self.cif_2pv7.get_array(key, dtype=dtype, gather=gather)


if __name__ == '__main__':
  absltest.main()
