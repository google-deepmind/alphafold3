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

import functools

from absl.testing import absltest
from alphafold3.common import resources
import msgpack


_CCD_MSGPACK_FILE = resources.filename(
    resources.ROOT / 'constants/converters/ccd.msgpack'
)


@functools.cache
def _load_ccd():
  with open(_CCD_MSGPACK_FILE, 'rb') as f:
    return msgpack.unpack(f)


class CcdGenTest(absltest.TestCase):

  def test_types_are_trivial(self):
    ccd = _load_ccd()
    self.assertEqual(type(ccd), dict)
    first = next(iter(ccd))
    self.assertEqual(type(first), str)
    table = ccd[first]
    self.assertEqual(type(table), dict)
    keys = iter(table)
    # Skip the first key, which is the cif name.
    next(keys)
    table_key = next(keys)
    self.assertEqual(type(table_key), str)
    table_value = table[table_key]
    self.assertEqual(type(table_value), list, f'{table}')
    self.assertNotEmpty(table_value)
    self.assertEqual(type(table_value[0]), str)

  def test_arg(self):
    ccd = _load_ccd()
    self.assertEqual(ccd['ARG']['data_'], ['ARG'])
    self.assertEqual(ccd['ARG']['_chem_comp.type'], ['L-PEPTIDE LINKING'])
    self.assertEqual(
        ccd['ARG']['_pdbx_chem_comp_descriptor.descriptor'],
        [
            'O=C(O)C(N)CCCN\\C(=[NH2+])N',
            'N[C@@H](CCCNC(N)=[NH2+])C(O)=O',
            'N[CH](CCCNC(N)=[NH2+])C(O)=O',
            'C(C[C@@H](C(=O)O)N)CNC(=[NH2+])N',
            'C(CC(C(=O)O)N)CNC(=[NH2+])N',
            'InChI=1S/C6H14N4O2/c7-4(5(11)12)2-1-3-10-6(8)9/h4H,1-3,7H2,(H,11,12)(H4,8,9,10)/p+1/t4-/m0/s1',
            'ODKSFYDXXFIFQN-BYPYZUCNSA-O',
        ],
    )

  def test_extended_5_char_codes(self):
    ccd = _load_ccd()
    for code in ('A1AV2', 'A1IJ5'):
      self.assertIn(code, ccd)
      self.assertEqual(ccd[code]['data_'], [code])


if __name__ == '__main__':
  absltest.main()
