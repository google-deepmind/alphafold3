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

from collections.abc import Sequence
import pickle
import sys

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.common import resources
from alphafold3.common.testing import data
from alphafold3.constants import chemical_components
from alphafold3.structure import mmcif
import numpy as np


_PARSERS_TESTDATA = data.Data(
    resources.ROOT / 'test_data/miniature_databases/pdb_mmcif'
)


def _read_test_cif_file(name):
  return _PARSERS_TESTDATA.load(f'{name}.cif')


class MmcifTest(parameterized.TestCase):

  @parameterized.parameters(
      (1, 'A'),
      (26, 'Z'),
      (27, 'AA'),
      (28, 'BA'),
      (29, 'CA'),
      (26**2 + 26, 'ZZ'),
      (26**2 + 27, 'AAA'),
      (26**2 + 28, 'BAA'),
      (26**2 + 26 + 27, 'ABA'),
      (26**3 + 26**2 + 27, 'AAAA'),
      (26**3 + 26**2 + 27 + 26 + 26**2, 'ABBA'),
  )
  def test_int_id_to_str_id(self, num, exp_int_id_to_str_id):
    self.assertEqual(mmcif.int_id_to_str_id(num), exp_int_id_to_str_id)

  @parameterized.parameters((0,), (-1), (-100))
  def test_int_id_to_str_id_negative(self, num):
    with self.assertRaises(ValueError):
      mmcif.int_id_to_str_id(num)

  @parameterized.parameters(
      ('A', 1),
      ('Z', 26),
      ('AA', 27),
      ('BA', 28),
      ('CA', 29),
      ('ZZ', 26**2 + 26),
      ('AAA', 26**2 + 27),
      ('BAA', 26**2 + 28),
      ('ABA', 26**2 + 26 + 27),
      ('AAAA', 26**3 + 26**2 + 27),
      ('ABBA', 26**3 + 26**2 + 27 + 26 + 26**2),
  )
  def test_str_id_to_int_id(self, str_id, exp_str_id_to_int_id):
    self.assertEqual(mmcif.str_id_to_int_id(str_id), exp_str_id_to_int_id)

  @parameterized.parameters('', ' ', ' A', 'A ', '1', 'a', 'aA')
  def test_str_id_to_int_id_failures(self, str_id):
    with self.assertRaises(ValueError):
      mmcif.str_id_to_int_id(str_id)

  def test_int_to_str_to_int_id_roundtrip(self):
    for i in range(1, 1000):
      str_id = mmcif.int_id_to_str_id(i)
      int_id = mmcif.str_id_to_int_id(str_id)
      self.assertEqual(int_id, i)

  def test_broken_file_quote_open(self):
    cif_string = _read_test_cif_file('3d1c_with_unclosed_quote')

    with self.assertRaisesRegex(ValueError, 'Line ended with quote open: .*'):
      mmcif.from_string(cif_string)

  def test_broken_file_no_data_header(self):
    cif_string = '# \n_entry.id   1WIB\n'

    with self.assertRaisesRegex(
        ValueError, 'The CIF file does not start with the data_ field.*'
    ):
      mmcif.from_string(cif_string)

  @parameterized.named_parameters(('empty', ''), ('blank', ' '))
  def test_empty_mmcif_file(self, cif_string):
    with self.assertRaisesRegex(ValueError, 'The CIF file must not be empty.*'):
      mmcif.from_string(cif_string)

  @parameterized.parameters(
      ('2pv7', '1'), ('3b9f', '1'), ('1wib', '15'), ('openmm_bonds', '1')
  )
  def test_get_bond_atom_indices(self, pdb_id, model_id):
    cif_string = _read_test_cif_file(pdb_id)
    mmcif_parsed = mmcif.from_string(cif_string)
    ptnr1_atoms, ptnr2_atoms = mmcif.get_bond_atom_indices(
        mmcif_parsed, model_id
    )
    self.assertLen(ptnr1_atoms, len(ptnr2_atoms))
    struct_conn = mmcif_parsed.extract_loop_as_list(prefix='_struct_conn.')
    self.assertLen(struct_conn, len(ptnr1_atoms))

    label_asym_id = mmcif_parsed['_atom_site.label_asym_id']
    label_seq_id = mmcif_parsed['_atom_site.label_seq_id']
    label_comp_id = mmcif_parsed['_atom_site.label_comp_id']
    label_atom_id = mmcif_parsed['_atom_site.label_atom_id']
    for struct_conn_row, ptnr1_atom_i, ptnr2_atom_i in zip(
        struct_conn, ptnr1_atoms, ptnr2_atoms, strict=True
    ):
      self.assertEqual(
          struct_conn_row['_struct_conn.ptnr1_label_asym_id'],
          label_asym_id[ptnr1_atom_i],
      )
      self.assertEqual(
          struct_conn_row['_struct_conn.ptnr2_label_asym_id'],
          label_asym_id[ptnr2_atom_i],
      )
      self.assertEqual(
          struct_conn_row['_struct_conn.ptnr1_label_seq_id'],
          label_seq_id[ptnr1_atom_i],
      )
      self.assertEqual(
          struct_conn_row['_struct_conn.ptnr2_label_seq_id'],
          label_seq_id[ptnr2_atom_i],
      )
      self.assertEqual(
          struct_conn_row['_struct_conn.ptnr1_label_comp_id'],
          label_comp_id[ptnr1_atom_i],
      )
      self.assertEqual(
          struct_conn_row['_struct_conn.ptnr2_label_comp_id'],
          label_comp_id[ptnr2_atom_i],
      )
      self.assertEqual(
          struct_conn_row['_struct_conn.ptnr1_label_atom_id'],
          label_atom_id[ptnr1_atom_i],
      )
      self.assertEqual(
          struct_conn_row['_struct_conn.ptnr2_label_atom_id'],
          label_atom_id[ptnr2_atom_i],
      )

  @parameterized.named_parameters(
      ('has_type_symbol', '2pv7'),
      ('missing_type_symbol', '5w6l_missing_type_symbol'),
  )
  def test_ret_of_infer_type_symbol(self, cif_file_name):
    cif_string = _read_test_cif_file(cif_file_name)
    mmcif_parsed = mmcif.from_string(cif_string)

    type_symbol = mmcif.get_or_infer_type_symbol(mmcif_parsed)

    exp_type_symbol = []
    for res_name, atom_name in zip(
        mmcif_parsed['_atom_site.label_comp_id'],
        mmcif_parsed['_atom_site.label_atom_id'],
        strict=True,
    ):
      element = chemical_components.type_symbol(
          ccd=chemical_components.Ccd(),
          res_name=res_name,
          atom_name=atom_name,
      )
      exp_type_symbol.append(element or '?')
    self.assertSequenceEqual(type_symbol, exp_type_symbol)

  def test_mmcif_pickleable(self):
    raw_data = dict(a=['1'], b=['2'], c=['3', '4'])
    mmcif_parsed = mmcif.Mmcif(raw_data)
    for protocol in range(pickle.HIGHEST_PROTOCOL, 1, -1):
      pickled = pickle.dumps(mmcif_parsed, protocol=protocol)  # Test pickling.
      unpickled = pickle.loads(pickled)  # Test unpickling.
      self.assertEqual(set(mmcif_parsed.keys()), set(unpickled.keys()))
      for k in mmcif_parsed.keys():
        self.assertEqual(mmcif_parsed[k], unpickled[k])

  def test_mmcif_extract_loop_as_list(self):
    raw_data = {
        'a.1': ['a.1.1', 'a.1.2', 'a.1.3'],
        'a.2': ['a.2.1', 'a.2.2', 'a.2.3'],
        'b.1': ['b.1.1'],
    }
    mmcif_parsed = mmcif.Mmcif(raw_data)
    with_prefix_a = mmcif_parsed.extract_loop_as_list(prefix='a')
    self.assertLen(with_prefix_a, 3)
    self.assertCountEqual(with_prefix_a[0], {'a.1': 'a.1.1', 'a.2': 'a.2.1'})
    self.assertCountEqual(with_prefix_a[1], {'a.1': 'a.1.2', 'a.2': 'a.2.2'})
    self.assertCountEqual(with_prefix_a[2], {'a.1': 'a.1.3', 'a.2': 'a.2.3'})
    with_prefix_b = mmcif_parsed.extract_loop_as_list(prefix='b')
    self.assertLen(with_prefix_b, 1)
    self.assertCountEqual(with_prefix_b[0], {'b.1': 'b.1.1'})
    with_prefix_not_in_mmcif = mmcif_parsed.extract_loop_as_list(prefix='c')
    self.assertEmpty(with_prefix_not_in_mmcif)

  def test_mmcif_extract_loop_as_list_missing_dot(self):
    raw_data = {'data_': ['TEST'], 'a.1': ['a', 'a'], 'a_other.1': ['a']}
    mmcif_parsed = mmcif.Mmcif(raw_data)

    with self.assertRaisesRegex(
        ValueError,
        "TEST: Columns do not have the same number of rows for prefix: 'a'",
    ):
      mmcif_parsed.extract_loop_as_list(prefix='a')

  def test_mmcif_extract_loop_as_dict(self):
    raw_data = {
        'a.ix': ['1', '2', '3'],
        'a.1': ['a.1.1', 'a.1.2', 'a.1.3'],
        'a.2': ['a.2.1', 'a.2.2', 'a.2.3'],
    }
    mmcif_parsed = mmcif.Mmcif(raw_data)
    indexed_loop = mmcif_parsed.extract_loop_as_dict(prefix='a', index='a.ix')
    self.assertLen(indexed_loop, 3)
    self.assertCountEqual(
        indexed_loop['1'], {'a.ix': '1', 'a.1': 'a.1.1', 'a.2': 'a.2.1'}
    )
    self.assertCountEqual(
        indexed_loop['2'], {'a.ix': '2', 'a.1': 'a.1.2', 'a.2': 'a.2.2'}
    )
    self.assertCountEqual(
        indexed_loop['3'], {'a.ix': '3', 'a.1': 'a.1.3', 'a.2': 'a.2.3'}
    )

  def test_mmcif_extract_loop_as_dict_bad_index(self):
    raw_data = {'data_': ['TEST'], 'a.index': ['1', '2', '3']}
    mmcif_parsed = mmcif.Mmcif(raw_data)
    with self.assertRaisesRegex(
        ValueError, "TEST: The loop index 'BAD' .* loop prefix 'a'"
    ):
      mmcif_parsed.extract_loop_as_dict(prefix='a', index='BAD')
    with self.assertRaisesRegex(
        ValueError, "TEST: The index column 'a.BAD' .* prefix 'a'"
    ):
      mmcif_parsed.extract_loop_as_dict(prefix='a', index='a.BAD')

  def test_get_internal_to_author_chain_id_map(self):
    cif_string = _read_test_cif_file('2pv7')
    mmcif_parsed = mmcif.from_string(cif_string)
    mapping = mmcif.get_internal_to_author_chain_id_map(mmcif_parsed)
    # There are 2 author chains: A and B. There are 8 internal chains A-H.
    # Internal chains A, B contain normal residues of the author A-B chains.
    # Internal chains C, D, E, F, G, H are HETATM-only chains.
    exp_mapping = {i: 'A' for i in 'ACDG'} | {i: 'B' for i in 'BEFH'}
    self.assertDictEqual(mapping, exp_mapping)

  def test_mmcif_get(self):
    raw_dict = dict(A=['a', 'b', 'c'], B=['A', 'B', 'C'])
    mmcif_parsed = mmcif.Mmcif(raw_dict)
    self.assertEqual(mmcif_parsed['A'], raw_dict['A'])
    self.assertEqual(mmcif_parsed['B'], raw_dict['B'])
    other = 'other'
    self.assertEqual(mmcif_parsed.get('A', other), raw_dict.get('A', other))
    self.assertEqual(mmcif_parsed.get('C', other), raw_dict.get('C', other))
    self.assertEqual(sys.getrefcount(mmcif_parsed.get('A', None)), 1)
    self.assertEqual(sys.getrefcount(mmcif_parsed.get('D', {})), 1)
    self.assertEqual(
        sys.getrefcount(mmcif_parsed.get('D', other)), sys.getrefcount(other)
    )

    self.assertEqual(mmcif_parsed.get('A'), raw_dict.get('A'))
    self.assertIs(mmcif_parsed.get('C'), raw_dict.get('C'))

  def test_mmcif_keys(self):
    raw_dict = dict(A=['a', 'b', 'c'], B=['A', 'B', 'C'])
    mmcif_parsed = mmcif.Mmcif(raw_dict)
    self.assertCountEqual(mmcif_parsed.keys(), raw_dict.keys())
    self.assertEqual(
        sys.getrefcount(mmcif_parsed.keys()), sys.getrefcount(raw_dict.keys())
    )

  def test_mmcif_key_error(self):
    mmcif_parsed = mmcif.Mmcif({})
    with self.assertRaisesRegex(KeyError, 'bad_key'):
      _ = mmcif_parsed['bad_key']

  def test_mmcif_values(self):
    raw_dict = dict(A=['a', 'b', 'c'], B=['A', 'B', 'C'])
    mmcif_parsed = mmcif.Mmcif(raw_dict)
    self.assertCountEqual(mmcif_parsed.values(), raw_dict.values())
    self.assertEqual(
        sys.getrefcount(mmcif_parsed.values()),
        sys.getrefcount(raw_dict.values()),
    )

  def test_mmcif_items(self):
    raw_dict = dict(A=['a', 'b', 'c'], B=['A', 'B', 'C'])
    mmcif_parsed = mmcif.Mmcif(raw_dict)
    self.assertCountEqual(mmcif_parsed.items(), raw_dict.items())
    self.assertEqual(
        sys.getrefcount(mmcif_parsed.items()), sys.getrefcount(raw_dict.items())
    )

  def test_mmcif_iter(self):
    raw_dict = dict(A=['a', 'b', 'c'], B=['A', 'B', 'C'])
    mmcif_parsed = mmcif.Mmcif(raw_dict)
    self.assertCountEqual(list(mmcif_parsed), list(raw_dict))

  def test_raises_key_error(self):
    raw_dict = dict(A=['a', 'b', 'c'], B=['A', 'B', 'C'])
    mmcif_parsed = mmcif.Mmcif(raw_dict)
    with self.assertRaises(KeyError):
      _ = mmcif_parsed['C']

  def test_mmcif_len(self):
    mmcif_parsed = mmcif.Mmcif({'a': ['1', '2', '3'], 'b': [], 'c': ['1', '2']})
    self.assertLen(mmcif_parsed, 3)

  def test_mmcif_update(self):
    raw_dict = {'a': ['1', '2', '3'], 'b': [], 'c': ['1', '2']}
    original_raw_dict = dict(raw_dict)
    mmcif_parsed = mmcif.Mmcif(raw_dict)

    other = {'a': ['1'], 'd': ['1']}
    updated = mmcif_parsed.copy_and_update(other)
    raw_dict.update(other)
    # Test that the Mmcif object was updated.
    self.assertCountEqual(mmcif_parsed.items(), original_raw_dict.items())
    self.assertCountEqual(updated.items(), raw_dict.items())

    # Test that removing the other dictionary doesn't modify the Mmcif object.
    del other
    self.assertCountEqual(mmcif_parsed.items(), original_raw_dict.items())
    self.assertCountEqual(updated.items(), raw_dict.items())

    # Test that removing the original raw dict doesn't modify the Mmcif object.
    copied_raw_dict = dict(raw_dict)
    del raw_dict
    self.assertCountEqual(mmcif_parsed.items(), original_raw_dict.items())
    self.assertCountEqual(updated.items(), copied_raw_dict.items())

  @parameterized.named_parameters(
      (
          'single_key_value',
          {'_key.a': ['value']},
          'data_TEST\n#\n_key.a value\n#\n',
      ),
      (
          'quote_key_value',
          {'_key.q': ['Never laugh at live dragons']},
          'data_TEST\n#\n_key.q "Never laugh at live dragons"\n#\n',
      ),
      (
          # Keys should be padded to the length of the longest.
          'multi_key_value',
          {'_key.a': ['value1'], '_key.bbb': ['value2']},
          'data_TEST\n#\n_key.a   value1\n_key.bbb value2\n#\n',
      ),
      (
          'plain_multiline_value',
          {'_key.m': ['ABCD\nEFGH']},
          'data_TEST\n#\n_key.m \n;ABCD\nEFGH\n;\n#\n',
      ),
      (
          'multiline_value_with_spaces',  # No quotes should be added.
          {'_key.m': ['AB CD\nEF GH']},
          'data_TEST\n#\n_key.m \n;AB CD\nEF GH\n;\n#\n',
      ),
      (
          'loop',
          {'_l.a': ['a1', 'aa2'], '_l.b': ['bb1', 'b2']},
          'data_TEST\n#\nloop_\n_l.a\n_l.b\na1  bb1 \naa2 b2  \n#\n',
      ),
      (
          # Non-alphabetical ordering for the _atom_site loop.
          'atom_site_loop',
          {
              '_atom_site.group_PDB': ['a1', 'a2'],
              '_atom_site.Cartn_x': ['b1', 'b2'],
              '_atom_site.aaa': ['c1', 'c2'],
          },
          (
              'data_TEST\n#\nloop_\n'
              '_atom_site.group_PDB\n_atom_site.Cartn_x\n_atom_site.aaa\n'
              'a1 b1 c1 \na2 b2 c2 \n#\n'
          ),
      ),
      (
          'atom_site_last',  # _atom_site is placed at the end, others sorted.
          {
              '_atom_site.group_PDB': ['a1', 'a2'],
              '_zzz.zzz': ['z1', 'z2'],
              '_entry.id': ['TEST'],
              '_bbb.bbb': ['b1', 'b2'],
          },
          (
              'data_TEST\n#\n'
              '_entry.id TEST\n#\n'
              'loop_\n_bbb.bbb\nb1 \nb2 \n#\n'
              'loop_\n_zzz.zzz\nz1 \nz2 \n#\n'
              'loop_\n_atom_site.group_PDB\na1 \na2 \n#\n'
          ),
      ),
  )
  def test_mmcif_to_string(self, raw_dict, exp_output):
    raw_dict.update({'data_': ['TEST']})
    mmcif_parsed = mmcif.Mmcif(raw_dict)
    self.assertEqual(mmcif_parsed.to_string(), exp_output)

  def test_invalid_data_block_name(self):
    raw_dict = {'data_': ['TEST WITH SPACES']}
    mmcif_parsed = mmcif.Mmcif(raw_dict)
    with self.assertRaises(ValueError):
      mmcif_parsed.to_string()

  def test_missing_data(self):
    raw_dict = {'a': ['b']}
    mmcif_parsed = mmcif.Mmcif(raw_dict)
    with self.assertRaises(ValueError):
      mmcif_parsed.to_string()

  def test_mmcif_to_string_and_back(self):
    # Test we can serialize to string and then load back the same data.
    cif_string = _read_test_cif_file('2pv7')
    mmcif_original = mmcif.from_string(cif_string)
    mmcif_as_string = mmcif_original.to_string()
    mmcif_from_to_string = mmcif.from_string(mmcif_as_string)

    self.assertCountEqual(mmcif_from_to_string.keys(), mmcif_original.keys())
    for key in mmcif_original.keys():
      self.assertListEqual(mmcif_from_to_string[key], mmcif_original[key])

  def test_mmcif_to_string_deterministic(self):
    cif_string = _read_test_cif_file('2pv7')
    mmcif_parsed = mmcif.from_string(cif_string)
    self.assertEqual(mmcif_parsed.to_string(), mmcif_parsed.to_string())

  def test_mmcif_get_data_name(self):
    cif_string = _read_test_cif_file('2pv7')
    mmcif_parsed = mmcif.from_string(cif_string)
    self.assertEqual(mmcif_parsed.get_data_name(), '2PV7')

  def test_parse_multi_data_cif(self):
    cif_string = (
        'data_000\n'
        '#\n'
        '_chem_comp.name     "methyl hydrogen carbonate"\n'
        '_chem_comp.type     NON-POLYMER\n'
        '_chem_comp.formula  "C2 H4 O3"\n'
        '#\n'
        '\n'
        'data_NAD\n'
        '#\n'
        '_chem_comp.id       NAD\n'
        '_chem_comp.name     NICOTINAMIDE-ADENINE-DINUCLEOTIDE\n'
        '_chem_comp.type     NON-POLYMER\n'
        '_chem_comp.formula  "C21 H27 N7 O14 P2"\n'
        '#\n'
        '\n'
    )
    data_blocks = mmcif.parse_multi_data_cif(cif_string)
    self.assertLen(data_blocks, 2)
    # Each Mmcif has its records + a special data_ record.
    self.assertLen(data_blocks['000'], 4)
    self.assertEqual(data_blocks['000']['data_'], ['000'])
    self.assertLen(data_blocks['NAD'], 5)
    self.assertEqual(data_blocks['NAD']['data_'], ['NAD'])

  def test_get_chain_type_by_entity_id(self):
    cif_string = _read_test_cif_file('2pv7')
    mmcif_parsed = mmcif.from_string(cif_string)
    mapping = mmcif.get_chain_type_by_entity_id(mmcif_parsed)
    exp_mapping = {
        '1': 'polypeptide(L)',  # Changed from 'polymer' using _entity_poly.
        '2': 'non-polymer',
        '3': 'non-polymer',
        '4': 'water',
    }

    self.assertEqual(mapping, exp_mapping)

  @parameterized.parameters(
      dict(methods=None, expected=None),
      dict(methods=('Method',), expected='method'),
      dict(methods=('Method1', 'Method2'), expected='method1,method2'),
  )
  def test_get_mmcif_method(
      self, methods: Sequence[str] | None, expected: str | None
  ):
    if methods is not None:
      mmcif_parsed = mmcif.Mmcif({'_exptl.method': methods})
    else:
      mmcif_parsed = mmcif.Mmcif({})
    actual = mmcif.get_experimental_method(mmcif_parsed)
    self.assertEqual(actual, expected)

  def test_get_release_date_string(self):
    raw_dict = {'_pdbx_audit_revision_history.revision_date': ['2019-01-20']}
    mmcif_parsed = mmcif.Mmcif(raw_dict)
    release_date = mmcif.get_release_date(mmcif_parsed)
    self.assertEqual(release_date, '2019-01-20')

  def test_get_release_date_list(self):
    raw_dict = {
        '_pdbx_audit_revision_history.revision_date': [
            '2019-01-20',
            '2017-02-20',
            '2018-03-20',
        ]
    }
    mmcif_parsed = mmcif.Mmcif(raw_dict)
    release_date = mmcif.get_release_date(mmcif_parsed)
    self.assertEqual(release_date, '2017-02-20')

  @parameterized.named_parameters(
      ('only_one', (None, '3.6', None), 3.6),
      ('multiple', ('1.7', None, '2.1'), 1.7),
      ('non_float', ('?', '4.2', None), 4.2),
      ('non_float_2', ('.', '4.2', None), 4.2),
      ('none', (None, None, None), None),
  )
  def test_get_resolution(self, resolutions, expected):
    res_fields = [
        '_refine.ls_d_res_high',
        '_em_3d_reconstruction.resolution',
        '_reflns.d_resolution_high',
    ]
    raw_dict = {
        res_field: [res]
        for res_field, res in zip(res_fields, resolutions, strict=True)
        if res is not None
    }
    mmcif_parsed = mmcif.Mmcif(raw_dict)
    resolution = mmcif.get_resolution(mmcif_parsed)
    self.assertEqual(resolution, expected)

  @parameterized.parameters(
      ('(1)', [['1']]),
      ('1', [['1']]),
      ('P', [['P']]),
      ('(P)', [['P']]),
      ('(1,2,5)', [['1'], ['2'], ['5']]),
      ('1-4', [['1'], ['2'], ['3'], ['4']]),
      ('(1-4)', [['1'], ['2'], ['3'], ['4']]),
      ('10-12', [['10'], ['11'], ['12']]),
      ('(10-12)', [['10'], ['11'], ['12']]),
      ('(1,2)(3,4)', [['1', '3'], ['1', '4'], ['2', '3'], ['2', '4']]),
      ('1,2,4-7,9', [['1'], ['2'], ['4'], ['5'], ['6'], ['7'], ['9']]),
      ('(1,2,4-7,9)', [['1'], ['2'], ['4'], ['5'], ['6'], ['7'], ['9']]),
      (
          '(1-3)(7-8)',
          [
              ['1', '7'],
              ['1', '8'],
              ['2', '7'],
              ['2', '8'],
              ['3', '7'],
              ['3', '8'],
          ],
      ),
      ('(P)(7-9)', [['P', '7'], ['P', '8'], ['P', '9']]),
      ('(P0)(12)', [['P0', '12']]),
      ('(X0)(9-11)', [['X0', '9'], ['X0', '10'], ['X0', '11']]),
      ('(X0)(1-2,3-4)', [['X0', '1'], ['X0', '2'], ['X0', '3'], ['X0', '4']]),
      (
          '(9-11)(10-12)',
          [
              ['9', '10'],
              ['9', '11'],
              ['9', '12'],
              ['10', '10'],
              ['10', '11'],
              ['10', '12'],
              ['11', '10'],
              ['11', '11'],
              ['11', '12'],
          ],
      ),
  )
  def test_parse_oper_expr(self, oper_expression, expected_result):
    result = mmcif.parse_oper_expr(oper_expression)
    self.assertLen(result, len(expected_result))
    for actual, expected in zip(result, expected_result, strict=True):
      self.assertSequenceEqual(actual, expected)

  @parameterized.parameters(
      ('((1-3)(7-10),(4,6)(5,11))',),  # Nested parentheses not supported.
      ('(1,2,3)(4,5,6)(7,8,9)',),  # Only up to 2 parentheses supported.
  )
  def test_parse_oper_expr_unsupported(self, oper_expression):
    with self.assertRaises(ValueError):
      mmcif.parse_oper_expr(oper_expression)

  @parameterized.parameters(
      (np.array([1.07]), 0, ['1']),
      (np.array([1.07]), 1, ['1.1']),
      (np.array([1.07]), 2, ['1.07']),
      (np.array([1.07]), 3, ['1.070']),
      (np.array([1.07]), 4, ['1.0700']),
  )
  def test_format_float_array(self, values, decimals, exp_output):
    self.assertSequenceEqual(
        mmcif.format_float_array(values=values, num_decimal_places=decimals),
        exp_output,
    )

  def test_format_float_array_error(self):
    with self.assertRaises(ValueError):
      mmcif.format_float_array(values=np.zeros((2, 2)), num_decimal_places=2)


if __name__ == '__main__':
  absltest.main()
