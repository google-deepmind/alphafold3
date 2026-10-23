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
from alphafold3.constants import chemical_components
from alphafold3.cpp import cif_dict
from alphafold3.cpp import mmcif_atom_site

_TEST_CIF_STR_MISSING = R"""
data_test
loop_
_atom_site.id
_atom_site.label_comp_id
_atom_site.label_atom_id
1 ARG CA   # Should map to C
2 FOO BAR  # Should map to unknown: ?
"""


_TEST_CIF_STR_PRESENT = R"""
data_test
loop_
_atom_site.id
_atom_site.label_comp_id
_atom_site.label_atom_id
_atom_site.type_symbol
1 ARG CA   R
2 FOO BAR  S
"""


class MmcifToArraysTest(absltest.TestCase):

  def test_get_or_infer_type_symbol_missing(self):
    cif = cif_dict.from_string(_TEST_CIF_STR_MISSING)
    type_symbol_fn = (
        lambda res_name, atom_name: chemical_components.type_symbol(
            chemical_components.Ccd(), res_name, atom_name
        )
    )
    self.assertEqual(
        mmcif_atom_site.get_or_infer_type_symbol(cif, type_symbol_fn),
        ['C', '?'],
    )

  def test_get_or_infer_type_symbol_present(self):
    cif = cif_dict.from_string(_TEST_CIF_STR_PRESENT)
    type_symbol_fn = (
        lambda res_name, atom_name: chemical_components.type_symbol(
            chemical_components.Ccd(), res_name, atom_name
        )
    )
    self.assertEqual(
        mmcif_atom_site.get_or_infer_type_symbol(cif, type_symbol_fn),
        ['R', 'S'],
    )


if __name__ == '__main__':
  absltest.main()
