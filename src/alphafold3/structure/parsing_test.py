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

import datetime
import os

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3 import structure
from alphafold3.common import resources
from alphafold3.common.testing import data
from alphafold3.constants import chemical_components
from alphafold3.constants import mmcif_names
from alphafold3.constants import residue_names
from alphafold3.structure import chemical_components as struc_chem_comps
from alphafold3.structure import mmcif
from alphafold3.structure import parsing
from alphafold3.structure import test_utils
import numpy as np

_STRUCTURE_TEST_DATA = data.Data(
    resources.ROOT / 'test_data/miniature_databases/pdb_mmcif'
)
_PARSERS_TEST_DATA = data.Data(
    resources.ROOT / 'test_data/miniature_databases/pdb_mmcif'
)
_EXAMPLE_SMILES = (
    'c1ccc(cc1)CSC[C@@H](C(=O)N[C@H](c2ccccc2)C(=O)O)NC(=O)CC[C@@H](C(=O)O)N'
)


def _walk_test_data(path: str, suffix: str):
  for dname, _, fnames in os.walk(_STRUCTURE_TEST_DATA.path(path)):
    for fname in fnames:
      if fname.endswith(suffix):
        yield dname, fname


class ParsingTest(test_utils.StructureTestCase):

  def make_test_structure(self, name: str, **kwargs):
    data_string = _STRUCTURE_TEST_DATA.load(f'{name}.cif')
    return parsing.from_mmcif(data_string, **kwargs)

  def test_2pv7(self):
    struc = self.make_test_structure('2pv7', include_water=True)
    self.assertEqual(struc.name, '2PV7')
    self.assertEqual(struc.resolution, 2.0)
    self.assertEqual(struc.release_date, datetime.date(2007, 5, 22))
    self.assertEqual(struc.structure_method, 'X-RAY DIFFRACTION')

    self.assertEqual(struc.num_chains, 8)
    self.assertEqual(struc.filter(chain_type='polypeptide(L)').num_chains, 2)
    self.assertEqual(struc.filter(chain_type='water').num_chains, 2)
    self.assertEqual(struc.filter_to_entity_type(ligand=True).num_chains, 4)
    self.assertEqual(
        struc.filter(chain_id='A').num_residues(count_unresolved=False), 276
    )
    self.assertEqual(
        struc.filter(chain_id='B').num_residues(count_unresolved=False), 280
    )
    self.assertEqual(
        struc.filter(chain_id='C').num_residues(count_unresolved=False), 1
    )
    self.assertEqual(
        struc.filter(chain_id='D').num_residues(count_unresolved=False), 1
    )
    self.assertEqual(
        struc.filter(chain_id='E').num_residues(count_unresolved=False), 1
    )
    self.assertEqual(
        struc.filter(chain_id='F').num_residues(count_unresolved=False), 1
    )
    self.assertEqual(
        struc.filter(chain_id='G').num_residues(count_unresolved=False), 156
    )
    self.assertEqual(
        struc.filter(chain_id='H').num_residues(count_unresolved=False), 237
    )
    self.assertEqual(
        struc.filter(chain_id='A').num_residues(count_unresolved=True), 298
    )
    self.assertEqual(
        struc.filter(chain_id='B').num_residues(count_unresolved=True), 298
    )
    self.assertEqual(
        struc.filter(chain_id='C').num_residues(count_unresolved=True), 1
    )
    self.assertEqual(
        struc.filter(chain_id='D').num_residues(count_unresolved=True), 1
    )
    self.assertEqual(
        struc.filter(chain_id='E').num_residues(count_unresolved=True), 1
    )
    self.assertEqual(
        struc.filter(chain_id='F').num_residues(count_unresolved=True), 1
    )
    self.assertEqual(
        struc.filter(chain_id='G').num_residues(count_unresolved=True), 156
    )
    self.assertEqual(
        struc.filter(chain_id='H').num_residues(count_unresolved=True), 237
    )
    self.assertEqual(struc.filter(chain_id='A').num_atoms, 2205)
    self.assertEqual(struc.filter(chain_id='B').num_atoms, 2222)
    self.assertEqual(struc.author_naming_scheme.auth_asym_id['A'], 'A')
    self.assertEqual(struc.author_naming_scheme.auth_seq_id['A'][216], '295')
    self.assertIsNone(struc.author_naming_scheme.insertion_code['A'][216])
    self.assertEqual(struc.author_naming_scheme.entity_id['B'], '1')
    self.assertEqual(
        struc.author_naming_scheme.entity_desc['1'].upper(),
        'T-PROTEIN [INCLUDES: CHORISMATE MUTASE (EC 5.4.99.5) '
        '(CM) AND PREPHENATE DEHYDROGENASE (EC 1.3.1.12) (PDH)]',
    )
    self.assertEqual(struc.all_residues['A'][1], ('MSE', 2))  # Unresolved.
    self.assertEqual(struc.author_naming_scheme.auth_seq_id['A'][2], '81')

    single_atom_struc = struc.filter(chain_id='B', res_id=216, atom_name='CB')
    atom_dict = next(single_atom_struc.iter_atoms())
    self.assertEqual(atom_dict['chain_type'], 'polypeptide(L)')
    self.assertEqual(atom_dict['res_name'], 'ILE')
    self.assertEqual(atom_dict['atom_element'], 'C')
    self.assertAlmostEqual(atom_dict['atom_x'], 23.903, places=5)
    self.assertAlmostEqual(atom_dict['atom_y'], 54.350, places=5)
    self.assertAlmostEqual(atom_dict['atom_z'], 43.732, places=5)
    self.assertAlmostEqual(atom_dict['atom_b_factor'], 27.58, places=5)
    self.assertAlmostEqual(atom_dict['atom_occupancy'], 1.0, places=5)

    ligand_atom_struc = struc.filter(chain_id='C', atom_name='OXT')
    atom_dict = next(ligand_atom_struc.iter_atoms())
    self.assertEqual(atom_dict['chain_type'], 'non-polymer')
    self.assertEqual(atom_dict['res_name'], 'TYR')
    self.assertEqual(atom_dict['res_id'], 1)
    self.assertEqual(atom_dict['atom_element'], 'O')
    self.assertAlmostEqual(atom_dict['atom_x'], 14.758, places=5)
    self.assertAlmostEqual(atom_dict['atom_y'], 49.079, places=5)
    self.assertAlmostEqual(atom_dict['atom_z'], 56.052, places=5)
    self.assertAlmostEqual(atom_dict['atom_b_factor'], 41.11, places=5)
    self.assertAlmostEqual(atom_dict['atom_occupancy'], 1.0, places=5)
    self.assertEqual(
        ligand_atom_struc.author_naming_scheme.auth_asym_id['C'], 'A'
    )
    self.assertEqual(
        ligand_atom_struc.author_naming_scheme.auth_seq_id['C'][1], '401'
    )

    ligand_atom_struc = struc.filter(chain_id='D', atom_name='C7N')
    atom_dict = next(ligand_atom_struc.iter_atoms())
    self.assertEqual(atom_dict['chain_type'], 'non-polymer')
    self.assertEqual(atom_dict['res_name'], 'NAD')
    self.assertEqual(atom_dict['res_id'], 1)
    self.assertEqual(atom_dict['atom_element'], 'C')
    self.assertAlmostEqual(atom_dict['atom_x'], 13.883, places=5)
    self.assertAlmostEqual(atom_dict['atom_y'], 42.814, places=5)
    self.assertAlmostEqual(atom_dict['atom_z'], 59.653, places=5)
    self.assertAlmostEqual(atom_dict['atom_b_factor'], 45.46, places=5)
    self.assertAlmostEqual(atom_dict['atom_occupancy'], 1.0, places=5)
    self.assertEqual(
        ligand_atom_struc.author_naming_scheme.auth_asym_id['D'], 'A'
    )
    self.assertEqual(
        ligand_atom_struc.author_naming_scheme.auth_seq_id['D'][1], '400'
    )

    ligand_atom_struc = struc.filter(chain_id='E', atom_name='OXT')
    atom_dict = next(ligand_atom_struc.iter_atoms())
    self.assertEqual(atom_dict['chain_type'], 'non-polymer')
    self.assertEqual(atom_dict['res_name'], 'TYR')
    self.assertEqual(atom_dict['res_id'], 1)
    self.assertEqual(atom_dict['atom_element'], 'O')
    self.assertAlmostEqual(atom_dict['atom_x'], 20.779, places=5)
    self.assertAlmostEqual(atom_dict['atom_y'], 66.403, places=5)
    self.assertAlmostEqual(atom_dict['atom_z'], 44.005, places=5)
    self.assertAlmostEqual(atom_dict['atom_b_factor'], 34.79, places=5)
    self.assertAlmostEqual(atom_dict['atom_occupancy'], 1.0, places=5)
    self.assertEqual(
        ligand_atom_struc.author_naming_scheme.auth_asym_id['E'], 'B'
    )
    self.assertEqual(
        ligand_atom_struc.author_naming_scheme.auth_seq_id['E'][1], '401'
    )

    ligand_atom_struc = struc.filter(chain_id='F', atom_name='C7N')
    atom_dict = next(ligand_atom_struc.iter_atoms())
    self.assertEqual(atom_dict['chain_type'], 'non-polymer')
    self.assertEqual(atom_dict['res_name'], 'NAD')
    self.assertEqual(atom_dict['res_id'], 1)
    self.assertEqual(atom_dict['atom_element'], 'C')
    self.assertAlmostEqual(atom_dict['atom_x'], 23.893, places=5)
    self.assertAlmostEqual(atom_dict['atom_y'], 71.949, places=5)
    self.assertAlmostEqual(atom_dict['atom_z'], 40.170, places=5)
    self.assertAlmostEqual(atom_dict['atom_b_factor'], 43.81, places=5)
    self.assertAlmostEqual(atom_dict['atom_occupancy'], 1.0, places=5)
    self.assertEqual(
        ligand_atom_struc.author_naming_scheme.auth_asym_id['F'], 'B'
    )
    self.assertEqual(
        ligand_atom_struc.author_naming_scheme.auth_seq_id['F'][1], '400'
    )

    atom_dict = next(struc.filter(chain_id='G').iter_atoms())
    self.assertEqual(atom_dict['chain_type'], 'water')
    self.assertEqual(atom_dict['res_name'], 'HOH')
    self.assertEqual(atom_dict['res_id'], 1)
    self.assertEqual(atom_dict['atom_element'], 'O')
    self.assertAlmostEqual(atom_dict['atom_x'], 24.581, places=5)
    self.assertAlmostEqual(atom_dict['atom_y'], 65.8, places=5)
    self.assertAlmostEqual(atom_dict['atom_z'], 54.075, places=5)
    self.assertAlmostEqual(atom_dict['atom_b_factor'], 25.92, places=5)
    self.assertAlmostEqual(atom_dict['atom_occupancy'], 1.0, places=5)
    self.assertEqual(struc.author_naming_scheme.auth_asym_id['G'], 'A')
    self.assertEqual(struc.author_naming_scheme.auth_seq_id['G'][1], '402')

    atom_dict = next(struc.filter(chain_id='H').iter_atoms())
    self.assertEqual(atom_dict['chain_type'], 'water')
    self.assertEqual(atom_dict['res_name'], 'HOH')
    self.assertEqual(atom_dict['res_id'], 1)
    self.assertEqual(atom_dict['atom_element'], 'O')
    self.assertAlmostEqual(atom_dict['atom_x'], 32.234, places=5)
    self.assertAlmostEqual(atom_dict['atom_y'], 71.815, places=5)
    self.assertAlmostEqual(atom_dict['atom_z'], 47.447, places=5)
    self.assertAlmostEqual(atom_dict['atom_b_factor'], 26.38, places=5)
    self.assertAlmostEqual(atom_dict['atom_occupancy'], 1.0, places=5)
    self.assertEqual(struc.author_naming_scheme.auth_asym_id['H'], 'B')
    self.assertEqual(struc.author_naming_scheme.auth_seq_id['H'][1], '402')

  def test_3vvw(self):
    # Both chains in 3vvw have auth_seq_ids different to label_seq_ids, and
    # they do not share the same mapping, making it a good test case.
    struc = self.make_test_structure('3vvw')
    self.assertEqual(struc.name, '3VVW')
    self.assertEqual(struc.resolution, 2.5)
    self.assertEqual(struc.release_date, datetime.date(2013, 2, 27))
    self.assertEqual(struc.structure_method, 'X-RAY DIFFRACTION')

    self.assertEqual(struc.num_chains, 2)
    self.assertEqual(
        struc.filter(chain_id='A').num_residues(count_unresolved=False), 110
    )
    self.assertEqual(
        struc.filter(chain_id='B').num_residues(count_unresolved=False), 114
    )
    self.assertEqual(struc.filter(chain_id='A').num_atoms, 929)
    self.assertEqual(struc.filter(chain_id='B').num_atoms, 932)
    self.assertEqual(struc.author_naming_scheme.auth_asym_id['B'], 'B')
    self.assertEqual(struc.author_naming_scheme.auth_seq_id['A'][16], '34')
    self.assertIsNone(struc.author_naming_scheme.insertion_code['A'][16])
    self.assertEqual(struc.author_naming_scheme.entity_id['B'], '2')
    self.assertEqual(
        struc.author_naming_scheme.entity_desc['2'].upper(),
        'MICROTUBULE-ASSOCIATED PROTEINS 1A/1B LIGHT CHAIN 3C',
    )
    self.assertEqual(struc.all_residues['A'][0], ('GLY', 1))  # Unresolved.
    self.assertEqual(struc.author_naming_scheme.auth_seq_id['A'][1], '19')

    single_atom_struc = struc.filter(chain_id='B', res_id=32, atom_name='CD')
    atom_dict = next(single_atom_struc.iter_atoms())
    self.assertEqual(atom_dict['chain_type'], 'polypeptide(L)')
    self.assertEqual(atom_dict['res_name'], 'ARG')
    self.assertEqual(atom_dict['atom_element'], 'C')
    self.assertAlmostEqual(atom_dict['atom_x'], 70.368, places=5)
    self.assertAlmostEqual(atom_dict['atom_y'], 0.354, places=5)
    self.assertAlmostEqual(atom_dict['atom_z'], 2.050, places=5)
    self.assertAlmostEqual(atom_dict['atom_b_factor'], 34.88, places=5)
    self.assertAlmostEqual(atom_dict['atom_occupancy'], 1.0, places=5)

  def test_multimodel(self):
    struc = self.make_test_structure('1wib', model_id=parsing.ModelID.ALL)
    self.assertEqual(struc.coords.shape, (20, 1326, 3))
    self.assertEqual(struc.atom_b_factor.shape, (20, 1326))
    self.assertEqual(struc.atom_occupancy.shape, (20, 1326))
    self.assertEqual(struc.res_id.shape, (1326,))

  def test_bad_multimodel_first_model_ok(self):
    # 2icy is a multimodel structure that has bad ligand naming for models > 1.
    # Another examples are e.g. 2q44 or 4xq2.
    struc = self.make_test_structure(
        '2icy_trimmed', model_id=parsing.ModelID.FIRST
    )
    self.assertEqual(struc.coords.shape, (7228, 3))
    self.assertEqual(struc.atom_b_factor.shape, (7228,))
    self.assertEqual(struc.atom_occupancy.shape, (7228,))
    self.assertEqual(struc.res_id.shape, (7228,))

  def test_bad_multimodel(self):
    # 2icy is a multimodel structure that has bad ligand naming for models > 1.
    with self.assertRaisesRegex(
        ValueError, '(label_asym_id, auth_seq_id, res_name, ins_code)'
    ):
      self.make_test_structure('2icy_trimmed', model_id=2)

  def test_bad_multimodel_all_models_error(self):
    # 2icy is a multimodel structure that has bad ligand naming for models > 1.
    with self.assertRaisesRegex(
        ValueError, 'Each model must have the same number of atoms'
    ):
      self.make_test_structure('2icy_trimmed', model_id=parsing.ModelID.ALL)

  def test_mmcif_with_pdbx_branch_scheme_only(self):
    # 5nj3 has only _pdbx_branch_scheme, no _pdbx_nonpoly_scheme.
    struc = self.make_test_structure('5nj3')
    self.assertEqual(struc.chains, ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'))
    self.assertEqual(struc.num_residues(count_unresolved=True), 2_202)
    self.assertEqual(struc.num_residues(count_unresolved=False), 1_962)
    self.assertEqual(struc.num_atoms, 12_662)

  def test_3b9f(
      self,
  ):
    # 3b9f has a wide range of chain types and all of _pdbx_poly_seq_scheme,
    # _pdbx_nonpoly_scheme and _pdbx_branch_scheme are in the mmCIF.
    struc = self.make_test_structure('3b9f')
    self.assertEqual(struc.name, '3B9F')
    self.assertEqual(struc.resolution, 1.6)
    self.assertEqual(struc.release_date, datetime.date(2008, 4, 22))
    self.assertEqual(struc.structure_method, 'X-RAY DIFFRACTION')

    self.assertEqual(struc.num_chains, 13)
    self.assertEqual(
        struc.filter(chain_id='A').num_residues(count_unresolved=False), 44
    )
    self.assertEqual(
        struc.filter(chain_id='B').num_residues(count_unresolved=False), 253
    )
    self.assertEqual(
        struc.filter(chain_id='C').num_residues(count_unresolved=False), 356
    )
    self.assertEqual(
        struc.filter(chain_id='D').num_residues(count_unresolved=False), 2
    )
    self.assertEqual(
        struc.filter(chain_id='E').num_residues(count_unresolved=False), 2
    )
    self.assertEqual(
        struc.filter(chain_id='F').num_residues(count_unresolved=False), 1
    )
    self.assertEqual(
        struc.filter(chain_id='G').num_residues(count_unresolved=False), 1
    )
    self.assertEqual(
        struc.filter(chain_id='H').num_residues(count_unresolved=False), 1
    )
    self.assertEqual(
        struc.filter(chain_id='I').num_residues(count_unresolved=False), 1
    )
    self.assertEqual(
        struc.filter(chain_id='J').num_residues(count_unresolved=False), 1
    )
    self.assertEqual(
        struc.filter(chain_id='K').num_residues(count_unresolved=False), 1
    )
    self.assertEqual(
        struc.filter(chain_id='L').num_residues(count_unresolved=False), 1
    )
    self.assertEqual(
        struc.filter(chain_id='M').num_residues(count_unresolved=False), 1
    )
    self.assertEqual(struc.filter(chain_id='A').num_atoms, 362)
    self.assertEqual(struc.filter(chain_id='B').num_atoms, 2025)
    self.assertEqual(struc.filter(chain_id='C').num_atoms, 2726)
    self.assertEqual(struc.filter(chain_id='D').num_atoms, 24)
    self.assertEqual(struc.filter(chain_id='E').num_atoms, 31)
    self.assertEqual(struc.filter(chain_id='F').num_atoms, 5)
    self.assertEqual(struc.author_naming_scheme.auth_asym_id['A'], 'L')
    self.assertEqual(struc.author_naming_scheme.auth_seq_id['A'][1], '-2')
    self.assertEqual(struc.author_naming_scheme.auth_seq_id['A'][4], '1')
    self.assertEqual(struc.author_naming_scheme.insertion_code['A'][4], 'R')
    self.assertEqual(struc.author_naming_scheme.entity_id['E'], '5')
    self.assertEqual(
        struc.author_naming_scheme.entity_desc['5'],
        (
            '2-O-sulfo-alpha-L-idopyranuronic acid-(1-4)-2-deoxy-6-O-sulfo-2-'
            '(sulfoamino)-alpha-D-glucopyranose'
        ),
    )
    self.assertEqual(struc.author_naming_scheme.entity_desc['6'], 'SULFATE ION')
    self.assertEqual(struc.all_residues['A'][0], ('THR', 1))  # Unresolved.
    self.assertEqual(struc.all_residues['B'][149], ('ALA', 150))  # Unresolved.

    # A protein atom.
    single_atom_struc = struc.filter(chain_id='A', res_id=7, atom_name='OE1')
    atom_dict = next(single_atom_struc.iter_atoms())
    self.assertEqual(atom_dict['chain_type'], 'polypeptide(L)')
    self.assertEqual(atom_dict['res_name'], 'GLN')
    self.assertEqual(atom_dict['atom_element'], 'O')
    self.assertAlmostEqual(atom_dict['atom_x'], -39.497, places=5)
    self.assertAlmostEqual(atom_dict['atom_y'], 40.072, places=5)
    self.assertAlmostEqual(atom_dict['atom_z'], -54.952, places=5)
    self.assertAlmostEqual(atom_dict['atom_b_factor'], 23.86, places=5)
    self.assertAlmostEqual(atom_dict['atom_occupancy'], 0.21, places=5)

    # A branched ligand (oligosaccharide) atom.
    single_atom_struc = struc.filter(chain_id='E', res_id=2, atom_name='O6B')
    atom_dict = next(single_atom_struc.iter_atoms())
    self.assertEqual(atom_dict['chain_type'], 'branched')
    self.assertEqual(atom_dict['res_name'], 'IDS')
    self.assertEqual(atom_dict['atom_element'], 'O')
    self.assertAlmostEqual(atom_dict['atom_x'], -35.312, places=5)
    self.assertAlmostEqual(atom_dict['atom_y'], 27.972, places=5)
    self.assertAlmostEqual(atom_dict['atom_z'], -26.818, places=5)
    self.assertAlmostEqual(atom_dict['atom_b_factor'], 72.14, places=5)
    self.assertAlmostEqual(atom_dict['atom_occupancy'], 1.0, places=5)

    # A non-branched ligand (sulphate ion) atom.
    single_atom_struc = struc.filter(chain_id='J', res_id=1, atom_name='S')
    atom_dict = next(single_atom_struc.iter_atoms())
    self.assertEqual(atom_dict['chain_type'], 'non-polymer')
    self.assertEqual(atom_dict['res_name'], 'SO4')
    self.assertEqual(atom_dict['atom_element'], 'S')
    self.assertAlmostEqual(atom_dict['atom_x'], -0.070, places=5)
    self.assertAlmostEqual(atom_dict['atom_y'], 9.018, places=5)
    self.assertAlmostEqual(atom_dict['atom_z'], -31.640, places=5)
    self.assertAlmostEqual(atom_dict['atom_b_factor'], 70.55, places=5)
    self.assertAlmostEqual(atom_dict['atom_occupancy'], 0.76, places=5)

  def test_5pti(self):
    # 5pti is a crystal structure of a single chain with two small ligands, PO4
    # and UNK, in heavy water solvent (DOD). All of these chains are given the
    # same author chain ID (A). The protein chain also contains altlocs.
    struc = self.make_test_structure('5pti', include_water=True)
    self.assertEqual(struc.name, '5PTI')
    self.assertEqual(struc.resolution, 1.0)
    self.assertEqual(struc.release_date, datetime.date(1984, 10, 29))
    self.assertEqual(
        struc.structure_method, 'X-RAY DIFFRACTION,NEUTRON DIFFRACTION'
    )

    self.assertEqual(struc.num_chains, 4)
    self.assertEqual(
        struc.filter(chain_id='A').num_residues(count_unresolved=False), 58
    )
    self.assertEqual(
        struc.filter(chain_id='B').num_residues(count_unresolved=False), 1
    )
    self.assertEqual(
        struc.filter(chain_id='C').num_residues(count_unresolved=False), 1
    )
    self.assertEqual(
        struc.filter(chain_id='D').num_residues(count_unresolved=False), 63
    )
    self.assertEqual(struc.filter(chain_id='A').num_atoms, 892)
    self.assertEqual(struc.filter(chain_id='B').num_atoms, 5)
    self.assertEqual(struc.filter(chain_id='C').num_atoms, 1)
    self.assertEqual(struc.filter(chain_id='D').num_atoms, 189)
    self.assertEqual(struc.author_naming_scheme.auth_asym_id['A'], 'A')
    self.assertEqual(struc.author_naming_scheme.auth_asym_id['B'], 'A')
    self.assertEqual(struc.author_naming_scheme.auth_asym_id['C'], 'A')
    self.assertEqual(struc.author_naming_scheme.auth_asym_id['D'], 'A')
    self.assertEqual(struc.author_naming_scheme.auth_seq_id['B'][1], '70')
    self.assertIsNone(struc.author_naming_scheme.insertion_code['B'][1])
    self.assertEqual(struc.author_naming_scheme.entity_id['B'], '2')
    self.assertEqual(
        struc.author_naming_scheme.entity_desc['2'].upper(), 'PHOSPHATE ION'
    )

    single_atom_struc = struc.filter(chain_id='C', res_id=1, atom_name='UNK')
    atom_dict = next(single_atom_struc.iter_atoms())
    self.assertEqual(atom_dict['chain_type'], 'non-polymer')
    self.assertEqual(atom_dict['res_name'], 'UNX')
    self.assertEqual(atom_dict['atom_element'], 'X')
    self.assertAlmostEqual(atom_dict['atom_x'], 23.599, places=5)
    self.assertAlmostEqual(atom_dict['atom_y'], 13.195, places=5)
    self.assertAlmostEqual(atom_dict['atom_z'], -7.043, places=5)
    self.assertAlmostEqual(atom_dict['atom_b_factor'], 45.42, places=5)
    self.assertAlmostEqual(atom_dict['atom_occupancy'], 1.0, places=5)

    single_atom_struc = struc.filter(chain_id='D', res_id=1, atom_name='D2')
    atom_dict = next(single_atom_struc.iter_atoms())
    self.assertEqual(atom_dict['chain_type'], 'water')
    self.assertEqual(atom_dict['res_name'], 'DOD')
    self.assertEqual(atom_dict['atom_element'], 'D')
    self.assertAlmostEqual(atom_dict['atom_x'], 24.153, places=5)
    self.assertAlmostEqual(atom_dict['atom_y'], 15.561, places=5)
    self.assertAlmostEqual(atom_dict['atom_z'], 8.326, places=5)
    self.assertAlmostEqual(atom_dict['atom_b_factor'], 22.60, places=5)
    self.assertAlmostEqual(atom_dict['atom_occupancy'], 0.4, places=5)

  def test_from_mmcif_include_other(self):
    # 5inz has a D-peptide polymer chain with ID "D".
    mmcif_str = _PARSERS_TEST_DATA.load('5inz.cif')
    struc = parsing.from_mmcif(mmcif_str, include_other=True)
    struc_no_other = parsing.from_mmcif(mmcif_str, include_other=False)

    # 3 L-peptides, 1 D-peptide, 6 + 7 ligands. We don't include water.
    with self.subTest('with_other'):
      self.assertEqual(struc.num_chains, 3 + 1 + 6 + 7)
      self.assertEqual(struc.num_residues(count_unresolved=True), 85)
      self.assertEqual(struc.num_residues(count_unresolved=False), 85)
      self.assertEqual(struc.num_atoms, 1194)
    with self.subTest('without_other'):
      self.assertEqual(struc_no_other.num_chains, 3 + 6 + 7)
      self.assertEqual(struc_no_other.num_residues(count_unresolved=True), 67)
      self.assertEqual(struc_no_other.num_residues(count_unresolved=False), 67)
      self.assertEqual(struc_no_other.num_atoms, 935)

  @parameterized.named_parameters([
      ('include_other_include_water', True, True, 8, 32 + 140, 794),
      ('include_other_exclude_water', True, False, 4, 32, 654),
      ('exclude_other_include_water', False, True, 4, 140, 140),
      ('exclude_other_exclude_water', False, False, 0, 0, 0),
  ])
  def test_no_standard_polymers(
      self, include_other, include_water, num_chains, num_residues, num_atoms
  ):
    # 3c1p has no standard protein/RNA/DNA chains, only other and water.
    struc = self.make_test_structure(
        name='3c1p',
        include_other=include_other,
        include_water=include_water,
    )
    self.assertEqual(struc.num_chains, num_chains)
    self.assertEqual(struc.num_residues(count_unresolved=True), num_residues)
    self.assertEqual(struc.num_residues(count_unresolved=False), num_residues)
    self.assertEqual(struc.num_atoms, num_atoms)

  @parameterized.named_parameters([
      ('no_other', False),
      ('include_other', True),
  ])
  def test_170d(self, include_other):
    # 170d is a small NMR structure of a non-standard nucleic acid double
    # helix. This should be included regardless of the value of include_other.
    struc = self.make_test_structure(name='170d', include_other=include_other)
    self.assertEqual(struc.name, '170D')
    self.assertIsNone(struc.resolution)
    self.assertEqual(struc.release_date, datetime.date(1994, 7, 31))
    self.assertEqual(struc.structure_method, 'SOLUTION NMR')

    self.assertEqual(struc.num_chains, 2)
    self.assertEqual(
        struc.filter(chain_id='A').num_residues(count_unresolved=False), 12
    )
    self.assertEqual(
        struc.filter(chain_id='B').num_residues(count_unresolved=False), 12
    )
    self.assertLen(struc.all_residues['A'], 12)
    self.assertLen(struc.all_residues['B'], 12)
    self.assertTupleEqual(struc.all_residues['A'][8], ('CAR', 9))
    self.assertTupleEqual(struc.all_residues['B'][8], ('CAR', 9))
    self.assertEqual(struc.author_naming_scheme.auth_asym_id['A'], 'A')
    self.assertEqual(struc.author_naming_scheme.auth_asym_id['B'], 'B')
    self.assertEqual(struc.author_naming_scheme.auth_seq_id['B'][1], '13')
    self.assertIsNone(struc.author_naming_scheme.insertion_code['B'][1])
    self.assertEqual(struc.author_naming_scheme.entity_id['B'], '1')
    # In the following check we have to remove spaces since the PDB version has
    # an extra space corresponding to a line break in the PDB file. These line
    # breaks usually match work breaks, however in this case it doesn't so the
    # extra space is erroneous, causing a mismatch with the mmCIF version.
    self.assertEqual(
        struc.author_naming_scheme.entity_desc['1'].upper().replace(' ', ''),
        "DNA/RNA(5'-D(*CP*GP*CP*GP*AP*AP*TP*T)-R(P*CAR)-D(P*GP*CP*G)-3')",
    )

    single_atom_struc = struc.filter(chain_id='B', res_id=9, atom_name="HO2'")
    atom_dict = next(single_atom_struc.iter_atoms())
    self.assertEqual(atom_dict['chain_type'], mmcif_names.DNA_RNA_HYBRID_CHAIN)
    self.assertEqual(atom_dict['res_name'], 'CAR')
    self.assertEqual(atom_dict['atom_element'], 'H')
    self.assertAlmostEqual(atom_dict['atom_x'], -6.614, places=5)
    self.assertAlmostEqual(atom_dict['atom_y'], -0.334, places=5)
    self.assertAlmostEqual(atom_dict['atom_z'], -7.746, places=5)
    self.assertAlmostEqual(atom_dict['atom_b_factor'], 0.17, places=5)
    self.assertAlmostEqual(atom_dict['atom_occupancy'], 1.0, places=5)

  def test_6ysc_heterogeneous_residues(self):
    # 6ysc has heterogeneous residues, i.e. residues where there are two
    # possible residue types for a single residue ID. These should be
    # resolved during parsing so that only one type is selected, and this
    # should be done consistently in both the atom lists and the lists of
    # all residues (including missing ones).
    struc = self.make_test_structure('6ysc', include_water=True)
    single_res = struc.filter(res_id=134, chain_id='B')

    # We don't mind which of the two residues was chosen here, we just need
    # to ensure only one residue name has been used.
    first_res_name = single_res.res_name[0]
    self.assertTrue(np.all(single_res.res_name == first_res_name))

    matches_in_all_res = [
        res_name
        for res_name, res_id in struc.all_residues['B']
        if res_id == 134
    ]
    self.assertLen(matches_in_all_res, 1)
    self.assertEqual(matches_in_all_res[0], first_res_name)

  def test_4uz3_heterogeneous_ligand_residues(self):
    # 4uz3 has heterogeneous ligand residues in the first residues of
    # chains D and F. They both have altloc B having higher occupancy, and
    # for D that is NAG, and for F it is NDG (the order of NAG/NDG is flipped).
    struc = self.make_test_structure('4uz3', include_water=True)
    np.testing.assert_array_equal(
        struc.filter(chain_id='D').group_by_residue.res_name,
        ['NAG', 'NAG', 'NAG', 'NAG', 'NAG', 'NAG'],
    )
    np.testing.assert_array_equal(
        struc.filter(chain_id='E').group_by_residue.res_name,
        ['NAG', 'NAG', 'NAG', 'NAG', 'NAG', 'NAG'],
    )
    np.testing.assert_array_equal(
        struc.filter(chain_id='F').group_by_residue.res_name,
        ['NDG', 'NAG', 'NAG', 'NAG', 'NAG', 'NAG'],
    )
    np.testing.assert_array_equal(
        struc.filter(chain_id='G').group_by_residue.res_name,
        ['NDG'],
    )
    expected_residues = {
        'D': [
            ('NAG', 1),  # Heterogeneous
            ('NAG', 2),
            ('NAG', 3),
            ('NAG', 4),
            ('NAG', 5),
            ('NAG', 6),
        ],
        'E': [
            ('NAG', 1),
            ('NAG', 2),
            ('NAG', 3),
            ('NAG', 4),
            ('NAG', 5),
            ('NAG', 6),
        ],
        'F': [
            ('NDG', 1),  # Heterogeneous
            ('NAG', 2),
            ('NAG', 3),
            ('NAG', 4),
            ('NAG', 5),
            ('NAG', 6),
        ],
        'G': [('NDG', 1)],
    }
    self.assertEqual(struc.all_residues, struc.all_residues | expected_residues)

  @parameterized.product(fix_mse_residues=(True, False))
  def test_heterogeneous_mse_residue(self, fix_mse_residues):
    struc = self.make_test_structure('5by8', fix_mse_residues=fix_mse_residues)
    res_name = 'MET' if fix_mse_residues else 'MSE'
    self.assertEqual(struc.all_residues['A'][107], ('GLU', 108))
    self.assertEqual(struc.all_residues['A'][108], (res_name, 109))
    self.assertEqual(struc.all_residues['A'][109], ('ASN', 110))
    atom = struc.filter(chain_id='A', res_id=109, atom_name='CA')
    self.assertAlmostEqual(atom.atom_x[0], 8.826, places=3)

  def test_fix_mse(self):
    struc_unfixed = self.make_test_structure('2pv7', fix_mse_residues=False)
    struc_no_mses = self.make_test_structure('2pv7', fix_mse_residues=True)

    mse_residue_ids = (93, 109, 123, 131, 173, 215, 232)

    changed_se_atom_count = 0
    for chain_id in ('A', 'B'):
      # Check data inferred from _pdbx_poly_seq_scheme has been fixed.
      for res in struc_unfixed.all_residues[chain_id]:
        if res[1] in mse_residue_ids:
          self.assertEqual(res[0], 'MSE')

      for res in struc_no_mses.all_residues[chain_id]:
        if res[1] in mse_residue_ids:
          self.assertEqual(res[0], 'MET')

      # Check that atom site data has been fixed.
      for res_id in mse_residue_ids:
        chain_unfixed = struc_unfixed.filter(chain_id=chain_id, res_id=res_id)
        chain_no_mses = struc_no_mses.filter(chain_id=chain_id, res_id=res_id)

        for atom_unfixed, atom_no_mse in zip(
            chain_unfixed.iter_atoms(), chain_no_mses.iter_atoms(), strict=True
        ):
          self.assertEqual(atom_unfixed['res_name'], 'MSE')
          self.assertEqual(atom_no_mse['res_name'], 'MET')
          if atom_unfixed['atom_name'] == 'SE':
            changed_se_atom_count += 1
            self.assertEqual(atom_no_mse['atom_name'], 'SD')
            self.assertEqual(atom_unfixed['atom_name'], 'SE')
            self.assertEqual(atom_no_mse['atom_element'], 'S')
          else:
            atom_unfixed = dict(atom_unfixed.items())
            atom_no_mse = dict(atom_no_mse.items())
            del atom_unfixed['res_name']
            del atom_no_mse['res_name']
            self.assertDictEqual(atom_unfixed, atom_no_mse)
    self.assertEqual(changed_se_atom_count, 2 * len(mse_residue_ids))
    ccdata = struc_no_mses.chemical_components_data
    assert ccdata is not None
    self.assertIn('MSE', ccdata.chem_comp)
    self.assertIn('MET', ccdata.chem_comp)

  @parameterized.named_parameters(('chain_J', 'J'), ('chain_K', 'K'))
  def test_fix_unknown_dna(self, dna_chain_id):
    # 5tzs has 2 DNA chains: J and K.
    struc_unfixed = self.make_test_structure(
        '5tzs_trimmed', fix_unknown_dna=False
    )
    struc_fixed = self.make_test_structure('5tzs_trimmed', fix_unknown_dna=True)

    # Check all_residues.
    all_res = [res[0] for res in struc_unfixed.all_residues[dna_chain_id]]
    with self.subTest('unfixed_all_res_uses_unk_rna'):
      self.assertEqual(all_res, [residue_names.UNK_RNA] * len(all_res))

    all_res = [res[0] for res in struc_fixed.all_residues[dna_chain_id]]
    with self.subTest('fixed_all_res_uses_unk_dna'):
      self.assertEqual(all_res, [residue_names.UNK_DNA] * len(all_res))

    # Check res_name.
    with self.subTest('unfixed_res_name_uses_unk_rna'):
      self.assertSetEqual(
          set(struc_unfixed.filter(chain_id=dna_chain_id).res_name),
          {residue_names.UNK_RNA},
      )
    with self.subTest('fixed_res_name_uses_unk_dna'):
      self.assertSetEqual(
          set(struc_fixed.filter(chain_id=dna_chain_id).res_name),
          {residue_names.UNK_DNA},
      )
    with self.subTest('no_impact_on_unaffected'):
      dna_mask = struc_unfixed.chain_type == mmcif_names.DNA_CHAIN
      unk_rna_mask = struc_unfixed.res_name == residue_names.UNK_RNA
      self.assertSequenceEqual(
          list(struc_unfixed.res_name[~(dna_mask & unk_rna_mask)]),
          list(struc_fixed.res_name[~(dna_mask & unk_rna_mask)]),
      )
    with self.subTest('fixed_chemical_components'):
      ccdata = struc_fixed.chemical_components_data
      assert ccdata is not None
      self.assertIn('DN', ccdata.chem_comp)
    with self.subTest('unfixed_chemical_components'):
      ccdata = struc_unfixed.chemical_components_data
      assert ccdata is not None
      self.assertNotIn('DN', ccdata.chem_comp)
      self.assertIn('N', ccdata.chem_comp)

  @parameterized.named_parameters(
      ('with_water', True, ('A', 'B', 'C', 'D')),
      ('without_water', False, ('A', 'B', 'D')),
  )
  def test_water_in_branch_scheme(self, include_water, exp_chains):
    mmcif_str = _PARSERS_TEST_DATA.load('water_in_branch_scheme.cif')
    struc = parsing.from_mmcif(mmcif_str, include_water=include_water)

    self.assertEqual(struc.chains, exp_chains)
    self.assertEqual(include_water, mmcif_names.WATER in struc.chain_type)

  def test_chemical_components_data(self):
    struc = self.make_test_structure('2pv7', model_id=parsing.ModelID.ALL)
    # fmt: off
    self.assertSameElements(
        struc['chemical_components_data'].chem_comp.keys(),
        ('ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'HOH',
         'ILE', 'LEU', 'LYS', 'MET', 'MSE', 'NAD', 'PHE', 'PRO', 'SER', 'THR',
         'TRP', 'TYR', 'VAL'))
    # fmt: on

  def test_missing_entity_table(self):
    struc = self.make_test_structure('6vs4_missing_entity_table')
    self.assertTrue(np.all(struc.chain_type == 'polypeptide(L)'))
    self.assertEqual(struc.num_atoms, 7172)
    self.assertEqual(struc.num_residues(count_unresolved=False), 900)
    self.assertEqual(struc.num_residues(count_unresolved=True), 900)
    self.assertEqual(struc.num_chains, 4)

  def test_unknown_atoms_in_residue(self):
    struc = self.make_test_structure('6g2k_unknown_atoms_in_residue')
    assert struc.atom_element is not None
    self.assertSequenceEqual(
        struc.atom_element.tolist(),
        ['N', 'C', 'C', 'C', 'O', 'S', '?', '?', '?', '?'],
    )

  @parameterized.named_parameters(
      ('with_water', True), ('without_water', False)
  )
  def test_pymol_files_generated_from_mmcif(self, include_water):
    # This tests files generated by PyMOL after loading data into PyMOL from
    # mmCIF.
    # Source structure taken from PDB. We remove all_residues data since the
    # pymol structure lacks this information, so the fair comparison is to the
    # PDB entry without all_residues.
    src_struc = self.make_test_structure('2pv7', include_water=include_water)
    src_struc = src_struc.copy_and_update(residues=src_struc.present_residues)
    # Produced using 2pv7 mmCIF and opening + exporting it with PyMOL.
    struc = self.make_test_structure(
        'pymol_2pv7_from_mmcif', include_water=include_water
    )
    self.assertEqual(struc.num_atoms, src_struc.num_atoms)  # 4934
    self.assertEqual(
        struc.num_residues(count_unresolved=True),
        src_struc.num_residues(count_unresolved=True),
    )  # 953
    self.assertEqual(struc.num_chains, src_struc.num_chains)  # 8
    self.assertEqual(struc.chains, src_struc.chains)
    np.testing.assert_array_equal(struc.res_name, src_struc.res_name)
    # Coordinates do not match exactly, PyMOL reorders atoms within residues.
    src_res_arrays = src_struc.to_res_arrays(include_missing_residues=True)
    res_arrays = struc.to_res_arrays(include_missing_residues=True)
    np.testing.assert_array_equal(
        res_arrays.atom_positions, src_res_arrays.atom_positions
    )
    np.testing.assert_array_equal(
        res_arrays.atom_mask, src_res_arrays.atom_mask
    )

    # Check numer of atoms match for protein, ligand, water subsets separately.
    for filter_type in ({'protein': True}, {'ligand': True}, {'water': True}):
      struc_subset = struc.filter_to_entity_type(**filter_type)
      src_struc_subset = src_struc.filter_to_entity_type(**filter_type)
      self.assertEqual(struc_subset.num_atoms, src_struc_subset.num_atoms)
      self.assertEqual(
          struc_subset.num_residues(count_unresolved=True),
          src_struc_subset.num_residues(count_unresolved=True),
      )
      self.assertEqual(struc_subset.num_chains, src_struc_subset.num_chains)
      np.testing.assert_array_equal(
          struc_subset.chain_id, src_struc_subset.chain_id
      )
      np.testing.assert_array_equal(
          struc_subset.res_name, src_struc_subset.res_name
      )

  def test_pymol_mmcif_generated_from_pdb_multiple_chains(self):
    # When PyMOL exports cifs after loading from PDB it makes no attempt to
    # guess values for label_asym_id or label_entity_id, and instead writes
    # "." for all values. When there are multiple chains this leads to invalid
    # mmCIF files since there will almost certainly be two residues which share
    # the same res ID and chain ID, which is not supposed to happen.
    with self.assertRaises(ValueError):
      self.make_test_structure('pymol_2pv7_from_pdb')

  @parameterized.named_parameters(
      ('include_water', True, 8913, ('C', 'D')),
      ('exclude_water', False, 219, ()),
  )
  def test_openmm_unset_label_entity_id(
      self, include_water, exp_num_res, exp_water_chains
  ):
    mmcif_str = _PARSERS_TEST_DATA.load('openmm_unset_label_entity_id.cif')
    struc = parsing.from_mmcif(
        mmcif_str, include_water=include_water, include_bonds=False
    )
    protein_chains = struc.filter_to_entity_type(protein=True).chains
    ligand_chains = struc.filter_to_entity_type(ligand=True).chains
    water_chains = struc.filter_to_entity_type(water=True).chains

    self.assertEqual(struc.num_residues(count_unresolved=False), exp_num_res)
    self.assertEqual(struc.num_residues(count_unresolved=True), exp_num_res)
    self.assertEqual(protein_chains, ('A',))
    self.assertEqual(ligand_chains, ('B', 'E', 'F'))
    self.assertEqual(water_chains, exp_water_chains)

  @parameterized.parameters(('2grz', 8, 6), ('2g14', 4, 3))
  def test_bonds_in_multimodel_bad_altloc_chains(
      self, pdb_id, num_all_bonds, num_valid_bonds
  ):
    # Affects 2evw, 2g0v, 2g0x, 2g0z, 2g10, 2g11, 2g12, 2g14, 2grz, 2ntw.
    cif = mmcif.from_string(_PARSERS_TEST_DATA.load(f'{pdb_id}.cif'))
    struc = parsing.from_parsed_mmcif(cif, include_bonds=True)
    self.assertLen(cif['_struct_conn.id'], num_all_bonds)
    self.assertLen(struc.bonds, num_valid_bonds)

  def test_bond_parsing_with_no_atoms(self):
    cif = mmcif.from_string(_STRUCTURE_TEST_DATA.load('2kql_trimmed.cif'))
    # 2kql has 1 polypeptide(D) which is not included unless include_other=True.
    struc = parsing.from_parsed_mmcif(cif, include_bonds=True)
    self.assertIsNotNone(struc)
    self.assertEqual(struc.num_atoms, 0)
    self.assertEmpty(struc.bonds)

  def test_parsing_atomless_mmcif(self):
    mmcif_str = _PARSERS_TEST_DATA.load('2pv7.cif')
    struc = parsing.from_mmcif(mmcif_str, include_bonds=False)
    struc = struc.filter(mask=np.zeros(struc.num_atoms, dtype=bool))
    with self.assertRaises(parsing.NoAtomsError):
      structure.from_mmcif(struc.to_mmcif(), include_bonds=False)

  @parameterized.named_parameters(
      ('include_water', True), ('exclude_water', False)
  )
  def test_bad_water_entity(self, include_water):
    mmcif_str = _PARSERS_TEST_DATA.load('bad_water_entity.cif')
    with self.assertRaisesRegex(ValueError, 'non-water entity has water'):
      parsing.from_mmcif(mmcif_str, include_water=include_water)

  @parameterized.named_parameters(
      ('include_water', True), ('exclude_water', False)
  )
  def test_hydroxycysteine(self, include_water):
    # 8cmw has a DOD (heavy water) as part of a polymer chain.
    mmcif_str = _PARSERS_TEST_DATA.load('8cmw_trimmed.cif')
    struc = parsing.from_mmcif(
        mmcif_str, include_water=include_water, include_bonds=False
    )
    self.assertEqual(struc.all_residues['C'], [('CYS', 1), ('DOD', 2)])

  def test_from_mmcif_with_name(self):
    mmcif_str = _PARSERS_TEST_DATA.load('2pv7.cif')
    struc = parsing.from_mmcif(mmcif_str, name='custom', include_bonds=False)
    self.assertEqual(struc.name, 'custom')

  @parameterized.parameters('2pv7', '3b9f')
  def test_bonds_with_altlocs(self, pdb_id):
    # All of these PDBs have alt locs which get resolved when creating
    # Structure. We check that the bonds point to correct atom keys since they
    # are constructed before we resolve them.
    mmcif_str = _PARSERS_TEST_DATA.load(f'{pdb_id}.cif')
    cif = mmcif.from_string(mmcif_str)
    struc = parsing.from_parsed_mmcif(
        cif, include_bonds=True, include_other=True, include_water=True
    )
    assert struc.bonds is not None

    atom_key_to_index = {key: index for index, key in enumerate(struc.atom_key)}

    from_asym_id = []
    from_comp_id = []
    from_seq_id = []
    from_atom_id = []
    dest_asym_id = []
    dest_comp_id = []
    dest_seq_id = []
    dest_atom_id = []
    for from_atom_key, dest_atom_key in zip(
        struc.bonds.from_atom_key, struc.bonds.dest_atom_key
    ):
      from_atom_idx = atom_key_to_index[from_atom_key]
      dest_atom_idx = atom_key_to_index[dest_atom_key]

      from_asym_id.append(struc.chain_id[from_atom_idx])
      from_comp_id.append(struc.res_name[from_atom_idx])
      from_res_id = str(struc.res_id[from_atom_idx])
      if struc.chain_type[from_atom_idx] in mmcif_names.NON_POLYMER_CHAIN_TYPES:
        from_res_id = '.'
      from_seq_id.append(from_res_id)
      from_atom_id.append(struc.atom_name[from_atom_idx])

      dest_asym_id.append(struc.chain_id[dest_atom_idx])
      dest_comp_id.append(struc.res_name[dest_atom_idx])
      dest_res_id = str(struc.res_id[dest_atom_idx])
      if struc.chain_type[dest_atom_idx] in mmcif_names.NON_POLYMER_CHAIN_TYPES:
        dest_res_id = '.'
      dest_seq_id.append(dest_res_id)
      dest_atom_id.append(struc.atom_name[dest_atom_idx])

    self.assertEqual(cif['_struct_conn.ptnr1_label_asym_id'], from_asym_id)
    self.assertEqual(cif['_struct_conn.ptnr1_label_comp_id'], from_comp_id)
    self.assertEqual(cif['_struct_conn.ptnr1_label_seq_id'], from_seq_id)
    self.assertEqual(cif['_struct_conn.ptnr1_label_atom_id'], from_atom_id)
    self.assertEqual(cif['_struct_conn.ptnr2_label_asym_id'], dest_asym_id)
    self.assertEqual(cif['_struct_conn.ptnr2_label_comp_id'], dest_comp_id)
    self.assertEqual(cif['_struct_conn.ptnr2_label_seq_id'], dest_seq_id)
    self.assertEqual(cif['_struct_conn.ptnr2_label_atom_id'], dest_atom_id)

  def test_openmm_mmcif_bonds(self):
    mmcif_str = _PARSERS_TEST_DATA.load('openmm_bonds.cif')
    struc = parsing.from_mmcif(
        mmcif_str, include_bonds=True, include_water=True
    )
    with self.subTest('name_parsed'):
      self.assertEqual(struc.name, 'cell')
    with self.subTest('bonds_parsed'):
      self.assertLen(struc.bonds, 70)
    with self.subTest('atom_sites_parsed'):
      self.assertEqual(struc.num_atoms, 8899)
      self.assertEqual(struc.num_residues(count_unresolved=False), 555)
      self.assertEqual(struc.num_residues(count_unresolved=True), 555)
      self.assertEqual(struc.num_chains, 5)
    with self.subTest('round_trips_via_mmcif'):
      self.assertStructuresEqual(
          parsing.from_mmcif(
              struc.to_mmcif(), include_bonds=True, include_water=True
          ),
          struc,
          ignore_fields=[
              'all_residues',
              'chemical_components_data',
              # OpenMM uses 4 decimal places, but we write 3.
              'atom_x',
              'atom_y',
              'atom_z',
          ],
      )

  def test_missing_type_symbol(self):
    mmcif_str = _PARSERS_TEST_DATA.load('5w6l_missing_type_symbol.cif')
    struc = parsing.from_mmcif(mmcif_str, include_bonds=False)
    exp_element = []
    for res_name, atom_name in zip(struc.res_name, struc.atom_name):
      element = chemical_components.type_symbol(
          ccd=chemical_components.Ccd(),
          res_name=res_name,
          atom_name=atom_name,
      )
      exp_element.append(element or '?')
    self.assertSequenceEqual(list(struc.atom_element), exp_element)

  @parameterized.parameters(True, False)
  def test_no_label_entity_id(self, include_water):
    # Test label_entity_id is correctly filled in from the struct_asym table.
    data_string = _STRUCTURE_TEST_DATA.load('6gnx.cif')
    old_cif = mmcif.from_string(data_string)
    struc = parsing.from_parsed_mmcif(
        old_cif, include_water=include_water, include_bonds=False
    )
    new_cif = mmcif.from_string(struc.to_mmcif())
    self.assertNotIn('_atom_site.label_entity_id', old_cif)
    self.assertSameElements(
        new_cif['_atom_site.label_entity_id'], {'1', '2', '3', '4'}
    )

  def test_missing_pdbx_poly_seq_scheme(self):
    # This mmCIF has a missing _pdbx_poly_seq_scheme table with two chains which
    # have residues with the same internal ID, but mapping to a different author
    # residue ID.
    mmcif_str = _PARSERS_TEST_DATA.load('missing_pdbx_poly_seq_scheme.cif')
    struc = parsing.from_mmcif(mmcif_str)
    self.assertLen(struc.chains, 2)
    self.assertEqual(struc.num_residues(count_unresolved=True), 2)
    self.assertEqual(struc.num_residues(count_unresolved=False), 2)
    self.assertEqual(struc.num_atoms, 2)

  def test_previous_cif_outputs(self):
    """Tests that we can parse cifs produced by old versions of Structure.

    This ensures backwards compatibility when new parts of mmCIF files
    are used when writing/reading mmCIFs.
    """
    for d, name in _walk_test_data('previous_version_outputs', suffix='.cif'):
      with self.subTest(name):
        cif_path = os.path.join(d, name)
        with open(cif_path, 'rt') as f:
          cif_str = f.read()
        struc = parsing.from_mmcif(cif_str, include_bonds=False)
        self.assertGreater(struc.num_atoms, 0)

  def test_from_parsed_mmcif(self):
    cif_str = _STRUCTURE_TEST_DATA.load('2pv7.cif')
    cif = mmcif.from_string(cif_str)
    struc_from_cif_obj = parsing.from_parsed_mmcif(cif, include_bonds=False)
    struc_from_string = parsing.from_mmcif(cif_str, include_bonds=False)
    self.assertStructuresEqual(struc_from_string, struc_from_cif_obj)

  def test_from_res_arrays(self):
    atom_mask = np.array([[1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 0]], dtype=bool)
    res_name = np.array(['ALA', 'ALA', 'GLY'], dtype=object)
    chain_id = np.array(['A', 'B', 'B'], dtype=object)
    chain_entity_id = np.array(['1', '2', '2'], dtype=object)
    atom_name = np.array(['C', 'CA', 'O', 'N', 'CB'], dtype=object)[
        np.array([[0, 1, 2, 3], [4, 2, 3, 1], [0, 2, 1, 4]])
    ]
    atom_x = np.array(
        [
            [[1.0, 0.0, 0.0, 4.0], [0.0, 2.0, 0.0, 3.0], [0.0, 0.0, 7.0, 0.0]],
            [[8.0, 0.0, 0.0, 1.0], [0.0, 3.0, 0.0, 9.0], [0.0, 0.0, 5.0, 0.0]],
        ],
        dtype=np.float32,
    )

    struc = parsing.from_res_arrays(
        name='test_structure',
        res_name=res_name,
        chain_id=chain_id,
        chain_entity_id=chain_entity_id,
        atom_name=atom_name,
        atom_x=atom_x,
        atom_mask=atom_mask,
    )

    self.assertEqual(struc.name, 'test_structure')
    np.testing.assert_array_equal(
        struc.res_name,
        np.array(['ALA', 'ALA', 'ALA', 'ALA', 'GLY'], dtype=object),
        strict=True,
    )
    np.testing.assert_array_equal(
        struc.chain_id,
        np.array(['A', 'A', 'B', 'B', 'B'], dtype=object),
        strict=True,
    )
    np.testing.assert_array_equal(
        struc.chain_entity_id,
        np.array(['1', '1', '2', '2', '2'], dtype=object),
        strict=True,
    )
    np.testing.assert_array_equal(
        struc.atom_name,
        np.array(['C', 'N', 'O', 'CA', 'CA'], dtype=object),
        strict=True,
    )
    np.testing.assert_array_equal(
        struc.atom_x,
        np.array(
            [[1.0, 4.0, 2.0, 3.0, 7.0], [8.0, 1.0, 3.0, 9.0, 5.0]],
            dtype=np.float32,
        ),
        strict=True,
    )

  def test_from_res_arrays_non_alphabetical_chains(self):
    atom_mask = np.array([[1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 0]], dtype=bool)
    res_name = np.array(['ALA', 'ALA', 'GLY'], dtype=object)
    chain_id = np.array(['B', 'A', 'A'], dtype=object)
    atom_name = np.array(['C', 'CA', 'O', 'N', 'CB'], dtype=object)[
        np.array([[0, 1, 2, 3], [4, 2, 3, 1], [0, 2, 1, 4]])
    ]

    struc = parsing.from_res_arrays(
        name='test_structure',
        res_name=res_name,
        chain_id=chain_id,
        atom_name=atom_name,
        atom_mask=atom_mask,
    )
    np.testing.assert_array_equal(
        struc.chain_id,
        np.array(['B', 'B', 'A', 'A', 'A'], dtype=object),
        strict=True,
    )

  def test_from_res_arrays_non_contiguous_chains(self):
    atom_mask = np.array([[1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 0]], dtype=bool)
    res_name = np.array(['ALA', 'ALA', 'GLY'], dtype=object)
    chain_id = np.array(['A', 'B', 'A'], dtype=object)
    atom_name = np.array(['C', 'CA', 'O', 'N', 'CB'], dtype=object)[
        np.array([[0, 1, 2, 3], [4, 2, 3, 1], [0, 2, 1, 4]])
    ]

    with self.assertRaises(ValueError):
      parsing.from_res_arrays(
          name='test_structure',
          res_name=res_name,
          chain_id=chain_id,
          atom_name=atom_name,
          atom_mask=atom_mask,
      )

  def test_from_sequences(self):
    release_date = datetime.date(2007, 5, 22)
    struc = parsing.from_sequences_and_bonds(
        sequences=['AGAG', 'AGCU', 'AGCT', '(NA)(VWW)', '1TU', '(GLY)(MSE)'],
        chain_types=[
            mmcif_names.PROTEIN_CHAIN,
            mmcif_names.RNA_CHAIN,
            mmcif_names.DNA_CHAIN,
            mmcif_names.NON_POLYMER_CHAIN,
            mmcif_names.NON_POLYMER_CHAIN,
            mmcif_names.PROTEIN_CHAIN,
        ],
        sequence_formats=[
            parsing.SequenceFormat.FASTA,
            parsing.SequenceFormat.FASTA,
            parsing.SequenceFormat.FASTA,
            parsing.SequenceFormat.CCD_CODES,
            parsing.SequenceFormat.CCD_CODES,
            parsing.SequenceFormat.CCD_CODES,
        ],
        bonded_atom_pairs=None,
        ccd=chemical_components.Ccd(),
        release_date=release_date,
    )
    self.assertEqual(struc.release_date, release_date)
    np.testing.assert_array_equal(struc.chain_id, [
        'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'C', 'C', 'C', 'C', 'D',
        'D', 'E', 'F', 'F'
    ])  # pyformat: disable
    np.testing.assert_array_equal(struc.res_name, [
        'ALA', 'GLY', 'ALA', 'GLY', 'A', 'G', 'C', 'U', 'DA', 'DG', 'DC', 'DT',
        'NA', 'VWW', '1TU', 'GLY', 'MSE',
    ])  # pyformat: disable
    np.testing.assert_array_equal(
        struc.atom_name, [
            'CA', 'CA', 'CA', 'CA', "C1'", "C1'", "C1'", "C1'", "C1'", "C1'",
            "C1'", "C1'", 'NA', 'N', 'C1', 'N', 'N',
        ],
    )  # pyformat: disable
    np.testing.assert_array_equal(
        struc.res_id, [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 1, 1, 2]
    )
    self.assertSequenceEqual(
        struc.all_residues['A'],
        [('ALA', 1), ('GLY', 2), ('ALA', 3), ('GLY', 4)],
    )
    self.assertSequenceEqual(
        struc.all_residues['B'], [('A', 1), ('G', 2), ('C', 3), ('U', 4)]
    )
    self.assertSequenceEqual(
        struc.all_residues['C'], [('DA', 1), ('DG', 2), ('DC', 3), ('DT', 4)]
    )
    self.assertSequenceEqual(struc.all_residues['D'], [('NA', 1), ('VWW', 2)])
    self.assertIsNotNone(struc.chemical_components_data)
    assert struc.chemical_components_data is not None
    self.assertIn('1TU', struc.chemical_components_data.chem_comp)
    self.assertEqual(
        struc.chemical_components_data.chem_comp['1TU'].name,
        '4-(3,5-difluoro-4-hydroxybenzyl)-1,2-dimethyl-1H-imidazol-5-ol',
    )
    self.assertEqual(
        struc.chemical_components_data.chem_comp['1TU'].formula,
        'C12 H12 F2 N2 O2',
    )
    self.assertSameElements(
        set(struc.res_name), struc.chemical_components_data.chem_comp
    )

    parsed_struct = parsing.from_mmcif(struc.to_mmcif(), include_bonds=False)
    self.assertStructuresEqual(
        struc, parsed_struct, ignore_fields=['author_naming_scheme']
    )

  def test_from_sequences_with_smiles_ligands(self):
    struc = parsing.from_sequences_and_bonds(
        sequences=[
            'DM_VWW:O=C(O)C(N)CCC(=O)NC(C(=O)NC(c1ccccc1)C(=O)O)CSCc2ccccc2',
            'DM_MES:[O-]S(=O)(=O)CC[NH+]1CCOCC1',
        ],
        chain_types=[mmcif_names.NON_POLYMER_CHAIN] * 2,
        sequence_formats=[parsing.SequenceFormat.LIGAND_SMILES] * 2,
        bonded_atom_pairs=None,
        ccd=chemical_components.Ccd(),
    )
    np.testing.assert_array_equal(struc.chain_id, ['A', 'B'])
    np.testing.assert_array_equal(struc.res_name, ['DM_VWW', 'DM_MES'])
    np.testing.assert_array_equal(struc.atom_name, ['', ''])
    np.testing.assert_array_equal(struc.res_id, [1, 1])
    self.assertSequenceEqual(struc.all_residues['A'], [('DM_VWW', 1)])
    self.assertSequenceEqual(struc.all_residues['B'], [('DM_MES', 1)])
    assert struc.chemical_components_data is not None
    self.assertDictEqual(
        struc.chemical_components_data.chem_comp,
        {
            'DM_VWW': struc_chem_comps.ChemCompEntry(
                type='non-polymer',
                pdbx_smiles=(
                    'O=C(O)C(N)CCC(=O)NC(C(=O)NC(c1ccccc1)C(=O)O)CSCc2ccccc2'
                ),
            ),
            'DM_MES': struc_chem_comps.ChemCompEntry(
                type='non-polymer',
                pdbx_smiles='[O-]S(=O)(=O)CC[NH+]1CCOCC1',
            ),
        },
    )

    parsed_struct = parsing.from_mmcif(struc.to_mmcif(), include_bonds=False)
    self.assertStructuresEqual(
        parsed_struct, struc, ignore_fields=['author_naming_scheme']
    )

  def test_from_sequences_with_smiles_ligands_repetitions(self):
    struc = parsing.from_sequences_and_bonds(
        sequences=[
            'DM_MES:[O-]S(=O)(=O)CC[NH+]1CCOCC1',
            'DM_MES:[O-]S(=O)(=O)CC[NH+]1CCOCC1',
            'DM_MES:[O-]S(=O)(=O)CC[NH+]1CCOCC1',
        ],
        chain_types=[mmcif_names.NON_POLYMER_CHAIN] * 3,
        sequence_formats=[parsing.SequenceFormat.LIGAND_SMILES] * 3,
        bonded_atom_pairs=None,
        ccd=chemical_components.Ccd(),
    )
    assert struc.chemical_components_data is not None
    self.assertDictEqual(
        struc.chemical_components_data.chem_comp,
        {
            'DM_MES': struc_chem_comps.ChemCompEntry(
                type='non-polymer',
                pdbx_smiles='[O-]S(=O)(=O)CC[NH+]1CCOCC1',
            ),
        },
    )
    np.testing.assert_array_equal(struc.chain_id, ['A', 'B', 'C'])

  def test_from_sequences_with_o1_first_atoms(self):
    struc = parsing.from_sequences_and_bonds(
        sequences=['NAG', 'XYP', 'O'],  # first atoms: C1, O, O
        chain_types=[mmcif_names.NON_POLYMER_CHAIN] * 3,
        sequence_formats=[parsing.SequenceFormat.CCD_CODES] * 3,
        bonded_atom_pairs=None,
        ccd=chemical_components.Ccd(),
    )
    np.testing.assert_array_equal(struc.chain_id, ['A', 'B', 'C'])
    np.testing.assert_array_equal(struc.res_name, ['NAG', 'XYP', 'O'])
    np.testing.assert_array_equal(struc.atom_name, ['C1', 'C1', 'O'])

  def test_from_sequences_with_smiles_ligands_mismatch(self):
    with self.assertRaisesRegex(
        ValueError, 'Mismatching data for ligand DM_MES'
    ):
      _ = parsing.from_sequences_and_bonds(
          sequences=[
              'DM_MES:[O-]S(=O)(=O)CC[NH+]1CCOCC1',
              'DM_MES:C1COCC[NH+]1CCS(=O)(=O)[O-]',
          ],
          chain_types=[mmcif_names.NON_POLYMER_CHAIN] * 2,
          sequence_formats=[parsing.SequenceFormat.LIGAND_SMILES] * 2,
          bonded_atom_pairs=None,
          ccd=chemical_components.Ccd(),
      )

  def test_from_sequences_with_smiles_ligands_ccdname(self):
    with self.assertRaisesRegex(ValueError, 'Ligand name BEN is in CCD'):
      parsing.from_sequences_and_bonds(
          sequences=['BEN:c1ccccc1'],
          chain_types=[mmcif_names.NON_POLYMER_CHAIN],
          sequence_formats=[parsing.SequenceFormat.LIGAND_SMILES],
          bonded_atom_pairs=None,
          ccd=chemical_components.Ccd(),
      )

  def test_from_sequences_and_bonds(self):
    release_date = datetime.date(2007, 5, 22)
    struc = parsing.from_sequences_and_bonds(
        sequences=[
            'AGAG',
            'AGCU',
            'AGCT',
            '(NA)(VWW)',
            '1TU',
            '(GLY)(MSE)',
            f'DM_S:{_EXAMPLE_SMILES}',
        ],
        chain_types=[
            mmcif_names.PROTEIN_CHAIN,
            mmcif_names.RNA_CHAIN,
            mmcif_names.DNA_CHAIN,
            mmcif_names.NON_POLYMER_CHAIN,
            mmcif_names.NON_POLYMER_CHAIN,
            mmcif_names.PROTEIN_CHAIN,
            mmcif_names.NON_POLYMER_CHAIN,
        ],
        sequence_formats=[
            parsing.SequenceFormat.FASTA,
            parsing.SequenceFormat.FASTA,
            parsing.SequenceFormat.FASTA,
            parsing.SequenceFormat.CCD_CODES,
            parsing.SequenceFormat.CCD_CODES,
            parsing.SequenceFormat.CCD_CODES,
            parsing.SequenceFormat.LIGAND_SMILES,
        ],
        bonded_atom_pairs=[
            ((0, 1, 'N'), (3, 0, 'NA')),  # inter-chain.
            ((3, 0, 'NA'), (3, 1, 'SG')),  # intra-chain.
            ((4, 0, 'F8'), (4, 0, 'F7')),  # intra-residue.
        ],
        ccd=chemical_components.Ccd(),
        release_date=release_date,
        name='test_from_sequences_and_bonds',
        bond_type='covale',
    )
    self.assertEqual(struc.release_date, release_date)
    self.assertEqual(struc.name, 'test_from_sequences_and_bonds')
    self.assertEqual(struc.num_chains, 7)
    self.assertEqual(
        struc.num_residues(count_unresolved=True), 4 + 4 + 4 + 2 + 1 + 2 + 1
    )
    self.assertDictEqual(
        struc.chain_single_letter_sequence(),
        {
            'A': 'AGAG',
            'B': 'AGCU',
            'C': 'AGCT',
            'D': 'XX',
            'E': 'X',
            'F': 'GM',
            'G': 'X',
        },
    )
    self.assertDictEqual(
        struc.chain_res_name_sequence(),
        {
            'A': ('ALA', 'GLY', 'ALA', 'GLY'),
            'B': ('A', 'G', 'C', 'U'),
            'C': ('DA', 'DG', 'DC', 'DT'),
            'D': ('NA', 'VWW'),
            'E': ('1TU',),
            'F': ('GLY', 'MSE'),
            'G': ('DM_S',),
        },
    )
    self.assertDictEqual(
        struc.all_residues,
        {
            'A': [('ALA', 1), ('GLY', 2), ('ALA', 3), ('GLY', 4)],
            'B': [('A', 1), ('G', 2), ('C', 3), ('U', 4)],
            'C': [('DA', 1), ('DG', 2), ('DC', 3), ('DT', 4)],
            'D': [('NA', 1), ('VWW', 2)],
            'E': [('1TU', 1)],
            'F': [('GLY', 1), ('MSE', 2)],
            'G': [('DM_S', 1)],
        },
    )
    self.assertSameElements(
        struc.filter(chain_id='A', res_id=1).atom_name, ['CA']
    )
    self.assertSameElements(
        struc.filter(chain_id='A', res_id=2).atom_name, ['N']
    )
    self.assertSameElements(
        struc.filter(chain_id='A', res_id=3).atom_name, ['CA']
    )
    self.assertSameElements(
        struc.filter(chain_id='A', res_id=4).atom_name, ['CA']
    )
    self.assertSameElements(
        struc.filter(chain_id='B', res_id=1).atom_name, ["C1'"]
    )
    self.assertSameElements(
        struc.filter(chain_id='B', res_id=2).atom_name, ["C1'"]
    )
    self.assertSameElements(
        struc.filter(chain_id='B', res_id=3).atom_name, ["C1'"]
    )
    self.assertSameElements(
        struc.filter(chain_id='B', res_id=4).atom_name, ["C1'"]
    )
    self.assertSameElements(
        struc.filter(chain_id='C', res_id=1).atom_name, ["C1'"]
    )
    self.assertSameElements(
        struc.filter(chain_id='C', res_id=2).atom_name, ["C1'"]
    )
    self.assertSameElements(
        struc.filter(chain_id='C', res_id=3).atom_name, ["C1'"]
    )
    self.assertSameElements(
        struc.filter(chain_id='C', res_id=4).atom_name, ["C1'"]
    )
    self.assertSameElements(
        struc.filter(chain_id='D', res_id=1).atom_name, ['NA']
    )
    self.assertSameElements(
        struc.filter(chain_id='D', res_id=2).atom_name, ['SG']
    )
    self.assertSameElements(
        struc.filter(chain_id='E', res_id=1).atom_name, ['F8', 'F7']
    )
    self.assertSameElements(
        struc.filter(chain_id='F', res_id=1).atom_name, ['N']
    )
    self.assertSameElements(
        struc.filter(chain_id='F', res_id=2).atom_name, ['N']
    )
    self.assertSameElements(
        struc.filter(chain_id='G', res_id=1).atom_name, ['']
    )
    self.assertNotEmpty(struc.to_mmcif())

    exp_from_atoms = [
        {'chain_id': 'A', 'res_id': 2, 'atom_name': 'N', 'res_name': 'GLY'},
        {'chain_id': 'D', 'res_id': 1, 'atom_name': 'NA', 'res_name': 'NA'},
        {'chain_id': 'E', 'res_id': 1, 'atom_name': 'F8', 'res_name': '1TU'},
    ]
    exp_dest_atoms = [
        {'chain_id': 'D', 'res_id': 1, 'atom_name': 'NA', 'res_name': 'NA'},
        {'chain_id': 'D', 'res_id': 2, 'atom_name': 'SG', 'res_name': 'VWW'},
        {'chain_id': 'E', 'res_id': 1, 'atom_name': 'F7', 'res_name': '1TU'},
    ]
    for bond_i, (from_atom, dest_atom, bond_info) in enumerate(
        struc.iter_bonds()
    ):
      self.assertEqual(from_atom, from_atom | exp_from_atoms[bond_i])  # pyrefly: ignore[unsupported-operation]
      self.assertEqual(dest_atom, dest_atom | exp_dest_atoms[bond_i])  # pyrefly: ignore[unsupported-operation]
      self.assertEqual(bond_info['type'], 'covale')

  def test_from_sequences_and_bonds_smiles_with_bonded_atoms(self):
    release_date = datetime.date(2007, 5, 22)
    struc = parsing.from_sequences_and_bonds(
        sequences=[
            'AGAG',
            f'DM_B:{_EXAMPLE_SMILES}',
        ],
        chain_types=[
            mmcif_names.PROTEIN_CHAIN,
            mmcif_names.NON_POLYMER_CHAIN,
        ],
        sequence_formats=[
            parsing.SequenceFormat.FASTA,
            parsing.SequenceFormat.LIGAND_SMILES,
        ],
        bonded_atom_pairs=[
            ((0, 0, 'N'), (1, 0, 'O1')),  # inter-residue
        ],
        ccd=chemical_components.Ccd(),
        release_date=release_date,
        name='test_from_sequences_and_bonds',
        bond_type='covale',
    )
    self.assertEqual(struc.release_date, release_date)
    self.assertEqual(struc.name, 'test_from_sequences_and_bonds')
    self.assertEqual(struc.num_chains, 2)
    self.assertEqual(struc.num_residues(count_unresolved=True), 4 + 1)
    self.assertDictEqual(
        struc.chain_single_letter_sequence(), {'A': 'AGAG', 'B': 'X'}
    )
    self.assertDictEqual(
        struc.chain_res_name_sequence(),
        {'A': ('ALA', 'GLY', 'ALA', 'GLY'), 'B': ('DM_B',)},
    )
    self.assertDictEqual(
        struc.all_residues,
        {
            'A': [('ALA', 1), ('GLY', 2), ('ALA', 3), ('GLY', 4)],
            'B': [('DM_B', 1)],
        },
    )
    assert struc.atom_name is not None
    self.assertListEqual(
        struc.atom_name.tolist(), ['N', 'CA', 'CA', 'CA', 'O1']
    )
    self.assertNotEmpty(struc.to_mmcif())

    exp_from_atoms = [
        {'chain_id': 'A', 'res_id': 1, 'atom_name': 'N', 'res_name': 'ALA'},
    ]
    exp_dest_atoms = [
        {'chain_id': 'B', 'res_id': 1, 'atom_name': 'O1', 'res_name': 'DM_B'},
    ]
    for bond_i, (from_atom, dest_atom, bond_info) in enumerate(
        struc.iter_bonds()
    ):
      self.assertEqual(from_atom, from_atom | exp_from_atoms[bond_i])  # pyrefly: ignore[unsupported-operation]
      self.assertEqual(dest_atom, dest_atom | exp_dest_atoms[bond_i])  # pyrefly: ignore[unsupported-operation]
      self.assertEqual(bond_info['type'], 'covale')
    with self.subTest('chemical_components_data'):
      assert struc.chemical_components_data is not None
      expected = {
          'DM_B': struc_chem_comps.ChemCompEntry(
              type='non-polymer',
              name='?',
              pdbx_synonyms='?',
              formula='?',
              formula_weight='?',
              mon_nstd_flag='?',
              pdbx_smiles=_EXAMPLE_SMILES,
          )
      }
      self.assertEqual(
          struc.chemical_components_data.chem_comp,
          {**struc.chemical_components_data.chem_comp, **expected},
      )
    with self.subTest('smiles_survives_mmcif_roundtrip'):
      parsed_struc = parsing.from_mmcif(struc.to_mmcif(), include_bonds=True)
      self.assertStructuresEqual(
          parsed_struc,
          struc,
          ignore_fields=['author_naming_scheme', 'chemical_components_data'],
      )
      self.assertEqual(
          parsed_struc.chemical_components_data.chem_comp['DM_B'],  # pyrefly: ignore[missing-attribute]
          struc.chemical_components_data.chem_comp['DM_B'],
      )

  def test_from_sequences_and_bonds_with_smiles_ligands(self):
    struc = parsing.from_sequences_and_bonds(
        sequences=[
            'DM_V:O=C(O)C',
            'DM_M:[O-]S(=O)(=O)CC[NH+]',
        ],
        chain_types=[mmcif_names.NON_POLYMER_CHAIN] * 2,
        sequence_formats=[parsing.SequenceFormat.LIGAND_SMILES] * 2,
        bonded_atom_pairs=None,
        ccd=chemical_components.Ccd(),
    )
    with self.subTest('chain_id'):
      np.testing.assert_array_equal(struc.chain_id, ['A', 'B'])
    with self.subTest('res_name'):
      np.testing.assert_array_equal(struc.res_name, ['DM_V', 'DM_M'])
    with self.subTest('atom_name'):
      np.testing.assert_array_equal(struc.atom_name, ['', ''])
    with self.subTest('res_id'):
      np.testing.assert_array_equal(struc.res_id, [1, 1])
    with self.subTest('all_residues'):
      self.assertSequenceEqual(struc.all_residues['A'], [('DM_V', 1)])
      self.assertSequenceEqual(struc.all_residues['B'], [('DM_M', 1)])
    with self.subTest('chemical_components_data'):
      assert struc.chemical_components_data is not None
      self.assertDictEqual(
          struc.chemical_components_data.chem_comp,
          {
              'DM_V': struc_chem_comps.ChemCompEntry(
                  type='non-polymer',
                  pdbx_smiles='O=C(O)C',
              ),
              'DM_M': struc_chem_comps.ChemCompEntry(
                  type='non-polymer',
                  pdbx_smiles='[O-]S(=O)(=O)CC[NH+]',
              ),
          },
      )

    parsed_struct = parsing.from_mmcif(struc.to_mmcif(), include_bonds=False)
    self.assertStructuresEqual(
        parsed_struct, struc, ignore_fields=['author_naming_scheme']
    )

  def test_parse_multimodel_cif(self):
    mmcif_str = _PARSERS_TEST_DATA.load('TS_RR_C1.cif')
    struc_first = parsing.from_mmcif(
        mmcif_str, model_id=parsing.ModelID.FIRST, include_bonds=False
    )
    struc_all = parsing.from_mmcif(
        mmcif_str, model_id=parsing.ModelID.ALL, include_bonds=False
    )

    with self.subTest('first_model'):
      self.assertEqual(struc_first.num_chains, 3)
      self.assertEqual(struc_first.num_residues(count_unresolved=True), 3)
      self.assertEqual(struc_first.num_residues(count_unresolved=False), 3)
      self.assertEqual(struc_first.num_atoms, 67)
      self.assertEqual(struc_first.coords.shape, (67, 3))

    with self.subTest('all_model'):
      self.assertEqual(struc_all.num_chains, 3)
      self.assertEqual(struc_all.num_residues(count_unresolved=True), 3)
      self.assertEqual(struc_all.num_residues(count_unresolved=False), 3)
      self.assertEqual(struc_all.num_atoms, 67)
      self.assertEqual(struc_all.coords.shape, (239, 67, 3))

  def test_chain_with_only_unresolved_residues(self):
    # Chain U has only unresolved residues.
    mmcif_str = _PARSERS_TEST_DATA.load('7a4p_trimmed.cif')
    struc = parsing.from_mmcif(
        mmcif_str, include_other=True, include_bonds=False
    )
    unresolved_chain = struc.filter(
        chain_id='U', cascade_delete=structure.CascadeDelete.NONE
    )
    # Chain U has residue data.
    self.assertEqual(
        unresolved_chain.all_residues['U'], [('LEU', 1), ('PRO', 2), ('THR', 3)]
    )
    # Chain U has no atom data.
    self.assertEmpty(unresolved_chain.atom_name)

  @parameterized.named_parameters(
      ('shuffled', ('C', 'B', 'A')),
      ('same_order_as_chain_id', ('A', 'B', 'C')),
  )
  def test_from_atom_arrays_all_residues(self, chain_order):
    expected_all_residues = {
        'A': [('ALA', 1), ('ILE', 2), ('GLY', 3)],  # ILE is unresolved
        'B': [('DA', 1)],
        'C': [('UNK', 1), ('UNK', 2)],  # Chain C is resolved but unk.
    }
    expected_all_residues = {k: expected_all_residues[k] for k in chain_order}

    struc_from_atoms = parsing.from_atom_arrays(
        name='unset',
        res_id=np.array([1, 1, 1, 3, 1, 1, 2], dtype=np.int32),
        chain_id=np.array(['A', 'A', 'A', 'A', 'B', 'C', 'C'], dtype=object),
        atom_name=np.array(
            ['N', 'CA', 'C', 'CA', "C1'", 'CA', 'CA'], dtype=object
        ),
        atom_element=np.array(
            ['N', 'C', 'C', 'C', 'C', 'C', 'C'], dtype=object
        ),
        all_residues=expected_all_residues,
        chain_type=np.array(
            [
                mmcif_names.PROTEIN_CHAIN,
                mmcif_names.PROTEIN_CHAIN,
                mmcif_names.PROTEIN_CHAIN,
                mmcif_names.PROTEIN_CHAIN,
                mmcif_names.DNA_CHAIN,
                mmcif_names.PROTEIN_CHAIN,
                mmcif_names.PROTEIN_CHAIN,
            ],
            dtype=object,
        ),
        res_name=np.array(
            ['ALA', 'ALA', 'ALA', 'GLY', 'DA', 'UNK', 'UNK'], dtype=object
        ),
        atom_x=np.zeros(7, dtype=np.float32),
        atom_y=np.zeros(7, dtype=np.float32),
        atom_z=np.zeros(7, dtype=np.float32),
        atom_b_factor=np.zeros(7, dtype=np.float32),
        atom_occupancy=np.ones(7, dtype=np.float32),
    )

    self.assertDictEqual(expected_all_residues, struc_from_atoms.all_residues)

  @parameterized.named_parameters(
      ('res_id', [('ALA', 1), ('ILE', 2), ('GLY', 4)]),  # Bad GLY id.
      ('res_name', [('ALA', 1), ('ILE', 2), ('TYR', 3)]),  # Bad res 3 name.
      ('shuffle', [('GLY', 3), ('ILE', 2), ('ALA', 1)]),  # Order doesn't match.
  )
  def test_from_atom_arrays_all_residue_order_inconsistent(self, chain_a_res):
    all_residues = {
        'A': chain_a_res,
        'B': [('DA', 1)],
        'C': [('UNK', 1), ('UNK', 2)],  # Chain C is resolved but unk.
    }

    with self.assertRaises(ValueError):
      parsing.from_atom_arrays(
          name='unset',
          res_id=np.array([1, 1, 1, 3, 1, 1, 2], dtype=np.int32),
          chain_id=np.array(['A', 'A', 'A', 'A', 'B', 'C', 'C'], dtype=object),
          atom_name=np.array(
              ['N', 'CA', 'C', 'CA', "C1'", 'CA', 'CA'], dtype=object
          ),
          atom_element=np.array(
              ['N', 'C', 'C', 'C', 'C', 'C', 'C'], dtype=object
          ),
          all_residues=all_residues,
          chain_type=np.array(
              [
                  mmcif_names.PROTEIN_CHAIN,
                  mmcif_names.PROTEIN_CHAIN,
                  mmcif_names.PROTEIN_CHAIN,
                  mmcif_names.PROTEIN_CHAIN,
                  mmcif_names.DNA_CHAIN,
                  mmcif_names.PROTEIN_CHAIN,
                  mmcif_names.PROTEIN_CHAIN,
              ],
              dtype=object,
          ),
          res_name=np.array(
              ['ALA', 'ALA', 'ALA', 'GLY', 'DA', 'UNK', 'UNK'], dtype=object
          ),
          atom_x=np.zeros(7, dtype=np.float32),
          atom_y=np.zeros(7, dtype=np.float32),
          atom_z=np.zeros(7, dtype=np.float32),
          atom_b_factor=np.zeros(7, dtype=np.float32),
          atom_occupancy=np.ones(7, dtype=np.float32),
      )

  @parameterized.parameters('2pv7', '3b9f', '1wib', '5inz')
  def test_from_atom_arrays_real(self, pdb_id):
    mmcif_str = _PARSERS_TEST_DATA.load(f'{pdb_id}.cif')
    original_struc = parsing.from_mmcif(
        mmcif_str, include_bonds=True, include_other=True, include_water=True
    )

    struc_from_atoms = parsing.from_atom_arrays(
        res_id=original_struc.res_id,
        name=original_struc.name,
        release_date=original_struc.release_date,
        resolution=original_struc.resolution,
        structure_method=original_struc.structure_method,
        all_residues=original_struc.all_residues,
        bioassembly_data=original_struc.bioassembly_data,
        chemical_components_data=original_struc.chemical_components_data,
        bond_table=original_struc.bonds_table,
        chain_id=original_struc.chain_id,
        chain_type=original_struc.chain_type,
        res_name=original_struc.res_name,
        atom_key=original_struc.atom_key,
        atom_name=original_struc.atom_name,
        atom_element=original_struc.atom_element,
        atom_x=original_struc.atom_x,
        atom_y=original_struc.atom_y,
        atom_z=original_struc.atom_z,
        atom_b_factor=original_struc.atom_b_factor,
        atom_occupancy=original_struc.atom_occupancy,
    )
    without_author_naming_scheme = original_struc.reset_author_naming_scheme()
    self.assertStructuresEqual(without_author_naming_scheme, struc_from_atoms)


if __name__ == '__main__':
  absltest.main()
