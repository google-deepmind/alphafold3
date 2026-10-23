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
from alphafold3 import structure
from alphafold3.constants import chemical_components
from alphafold3.constants import mmcif_names
from alphafold3.structure import chemical_components as struc_chem_comps
from alphafold3.structure import mmcif
from alphafold3.structure import test_utils
import rdkit.Chem as rd_chem

ChemCompAtom = struc_chem_comps.ChemCompAtom
ChemCompBond = struc_chem_comps.ChemCompBond
ChemicalComponentsData = struc_chem_comps.ChemicalComponentsData
ChemCompEntry = struc_chem_comps.ChemCompEntry


class ChemCompEntryTest(parameterized.TestCase):

  def test_structure_to_cif_dict(self):
    smiles = 'Cc1c(c2cc(ccc2o1)NS(=O)(=O)c3ccc(cc3)Cl)C(=O)O'
    expected_atom_symbols = [
        'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'O', 'N', 'S', 'O', 'O',
        'C', 'C', 'C', 'C', 'C', 'C', 'Cl', 'C', 'O', 'O', 'H', 'H', 'H', 'H',
        'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
    ]  # pyformat: disable
    expected_atom_charges = ['0'] * 36
    struc = structure.from_sequences_and_bonds(
        sequences=[f'foo:{smiles}'],
        chain_types=[mmcif_names.NON_POLYMER_CHAIN],
        sequence_formats=[structure.SequenceFormat.LIGAND_SMILES],
        bonded_atom_pairs=None,
        ccd=chemical_components.Ccd(),
        name='struc_with_ligand',
    )
    chemical_components_data = struc.chemical_components_data
    self.assertIsNotNone(chemical_components_data)
    chem_comp_entry = chemical_components_data.chem_comp['foo']
    mol = chem_comp_entry.rdkit_mol
    mol = rd_chem.AddHs(mol)
    charges = [str(a.GetFormalCharge()) for a in mol.GetAtoms()]
    self.assertEqual(charges, expected_atom_charges)
    atom_idxs = [a.GetIdx() for a in mol.GetAtoms()]
    self.assertEqual(atom_idxs, list(range(36)))
    atom_symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    self.assertEqual(atom_symbols, expected_atom_symbols)
    self.assertLen(mol.GetBonds(), 38)


class ChemicalComponentsDataTest(parameterized.TestCase):

  def test_from_mmcif(self):
    mmcif_str = test_utils.read_test_cif_file('2pv7')
    parsed_mmcif = mmcif.from_string(mmcif_str)
    chemical_components_data = (
        struc_chem_comps.ChemicalComponentsData.from_mmcif(
            parsed_mmcif, fix_mse=False, fix_unknown_dna=False
        )
    )

    self.assertLen(chemical_components_data.chem_comp, 23)
    self.assertEqual(
        chemical_components_data.chem_comp['NAD'],
        ChemCompEntry(
            type='non-polymer',
            name='NICOTINAMIDE-ADENINE-DINUCLEOTIDE',
            pdbx_synonyms='?',
            formula='C21 H27 N7 O14 P2',
            formula_weight='663.425',
            mon_nstd_flag='.',
        ),
    )

  def test_to_mmcif(self):
    chemical_components_data = ChemicalComponentsData(
        chem_comp={
            'ALN': ChemCompEntry(
                type='L-PEPTIDE LINKING',
                name='Alanine',
                formula='C13 H13 N O2',
            ),
            'DM-X': ChemCompEntry(
                type='non-polymer',
                name='Secret',
                formula='TBD',
            ),
            'CHEMBL58951': ChemCompEntry(
                type='peptide-like',
                name='too long',
                formula='C23 H27 N3 O6 S',
            ),
        }
    )
    mmcif_dict = chemical_components_data.to_mmcif_dict()

    self.assertCountEqual(
        mmcif_dict.keys(),
        (
            '_chem_comp.id',
            '_chem_comp.type',
            '_chem_comp.name',
            '_chem_comp.pdbx_synonyms',
            '_chem_comp.formula',
            '_chem_comp.formula_weight',
            '_chem_comp.mon_nstd_flag',
        ),
    )
    self.assertSameStructure(
        mmcif_dict['_chem_comp.id'], ('ALN', 'CHEMBL58951', 'DM-X')
    )
    self.assertSameStructure(
        mmcif_dict['_chem_comp.name'], ('Alanine', 'too long', 'Secret')
    )
    self.assertSameStructure(
        mmcif_dict['_chem_comp.formula'],
        ('C13 H13 N O2', 'C23 H27 N3 O6 S', 'TBD'),
    )

  def test_to_mmcif_pdbx_smiles(self):
    chemical_components_data = ChemicalComponentsData(
        chem_comp={
            'ALN': ChemCompEntry(
                type='L-PEPTIDE LINKING',
                name='Alanine',
            ),
            'DM-X': ChemCompEntry(
                type='non-polymer',
                name='Secret',
                pdbx_smiles='under construction',
            ),
        }
    )
    mmcif_dict = chemical_components_data.to_mmcif_dict()

    self.assertCountEqual(
        mmcif_dict.keys(),
        (
            '_chem_comp.id',
            '_chem_comp.type',
            '_chem_comp.name',
            '_chem_comp.pdbx_synonyms',
            '_chem_comp.formula',
            '_chem_comp.formula_weight',
            '_chem_comp.mon_nstd_flag',
            '_chem_comp.pdbx_smiles',
        ),
    )
    self.assertSameStructure(mmcif_dict['_chem_comp.id'], ('ALN', 'DM-X'))
    self.assertSameStructure(
        mmcif_dict['_chem_comp.name'], ('Alanine', 'Secret')
    )
    self.assertSameStructure(
        mmcif_dict['_chem_comp.pdbx_smiles'], ('?', 'under construction')
    )

  def test_to_mmcif_empty(self):
    chemical_components_data = ChemicalComponentsData(chem_comp={})
    mmcif_dict = chemical_components_data.to_mmcif_dict()
    self.assertEmpty(mmcif_dict)

  def test_unset_smiles(self):
    vww_entry = ChemCompEntry(
        type='peptide-like',
        name='TEST',
        pdbx_synonyms='?',
        formula='C23 H27 N3 O6 S',
        formula_weight='473.542',
        mon_nstd_flag='y',
        pdbx_smiles=None,
    )
    exp_chem_comp_data = ChemicalComponentsData(chem_comp={'VWW': vww_entry})

    chem_comp_cif_dict = mmcif.Mmcif({
        '_chem_comp.id': ['VWW'],
        '_chem_comp.name': ['TEST'],
        '_chem_comp.type': ['peptide-like'],
        '_chem_comp.formula_weight': ['473.542'],
        '_chem_comp.formula': ['C23 H27 N3 O6 S'],
        '_chem_comp.mon_nstd_flag': ['y'],
        '_chem_comp.pdbx_synonyms': ['?'],
        '_chem_comp.pdbx_smiles': ['?'],
    })

    chem_comp_data = ChemicalComponentsData.from_mmcif(
        chem_comp_cif_dict, fix_mse=False, fix_unknown_dna=False
    )
    self.assertEqual(chem_comp_data, exp_chem_comp_data)

  def test_chem_comp_atoms(self):
    vww_entry = ChemCompEntry(
        type='peptide-like',
        name='TEST',
        pdbx_synonyms='?',
        formula='C23 H27 N3 O6 S',
        formula_weight='473.542',
        mon_nstd_flag='y',
        pdbx_smiles=None,
        chem_comp_atoms={
            'C1': ChemCompAtom(
                type_symbol='C',
                charge=0,
                leaving_atom_flag='N',
                model_ideal_x=1.0,
                model_ideal_y=2.0,
                model_ideal_z=3.0,
            ),
            'C2': ChemCompAtom(
                type_symbol='C',
                charge=1,
                leaving_atom_flag='N',
                model_ideal_x=4.0,
                model_ideal_y=5.0,
                model_ideal_z=6.0,
            ),
            'N3': ChemCompAtom(
                type_symbol='N',
                charge=-1,
                leaving_atom_flag='Y',
                model_ideal_x=7.0,
                model_ideal_y=8.0,
                model_ideal_z=9.0,
            ),
            'O4': ChemCompAtom(
                type_symbol='O',
                leaving_atom_flag='?',
                model_ideal_x=10.0,
                model_ideal_y=11.201,
                model_ideal_z=12.1,
            ),
            'S5': ChemCompAtom(
                type_symbol='S',
                charge=0,
                leaving_atom_flag='?',
                model_ideal_x=2.0,
                model_ideal_y=2.0,
                model_ideal_z=8.2,
            ),
        },
        chem_comp_bonds=[
            ChemCompBond(
                atom_id_1='C1',
                atom_id_2='C2',
                value_order='SING',
                aromatic_flag='N',
                stereo_config='N',
            ),
            ChemCompBond(
                atom_id_1='C1',
                atom_id_2='N3',
                value_order='DOUB',
                aromatic_flag='Y',
                stereo_config='Y',
            ),
            ChemCompBond(
                atom_id_1='C2',
                atom_id_2='O4',
                value_order='SING',
                aromatic_flag='?',
                stereo_config='?',
            ),
            ChemCompBond(
                atom_id_1='N3',
                atom_id_2='S5',
                value_order='ARO',
                aromatic_flag='?',
                stereo_config='?',
            ),
        ],
    )
    exp_chem_comp_data = ChemicalComponentsData(chem_comp={'VWW': vww_entry})

    chem_comp_cif_dict = mmcif.Mmcif({
        '_chem_comp.id': ['VWW'],
        '_chem_comp.name': ['TEST'],
        '_chem_comp.type': ['peptide-like'],
        '_chem_comp.formula_weight': ['473.542'],
        '_chem_comp.formula': ['C23 H27 N3 O6 S'],
        '_chem_comp.mon_nstd_flag': ['y'],
        '_chem_comp.pdbx_synonyms': ['?'],
        # Atom fields
        '_chem_comp_atom.comp_id': ['VWW', 'VWW', 'VWW', 'VWW', 'VWW'],
        '_chem_comp_atom.atom_id': ['C1', 'C2', 'N3', 'O4', 'S5'],
        '_chem_comp_atom.type_symbol': ['C', 'C', 'N', 'O', 'S'],
        '_chem_comp_atom.charge': ['0', '1', '-1', '?', '0'],
        '_chem_comp_atom.pdbx_ordinal': ['?', '?', '?', '?', '?'],
        '_chem_comp_atom.pdbx_leaving_atom_flag': ['N', 'N', 'Y', '?', '?'],
        '_chem_comp_atom.pdbx_model_Cartn_x_ideal': [
            '1.000',
            '4.000',
            '7.000',
            '10.000',
            '2.000',
        ],
        '_chem_comp_atom.pdbx_model_Cartn_y_ideal': [
            '2.000',
            '5.000',
            '8.000',
            '11.201',
            '2.000',
        ],
        '_chem_comp_atom.pdbx_model_Cartn_z_ideal': [
            '3.000',
            '6.000',
            '9.000',
            '12.100',
            '8.200',
        ],
        # Bond fields
        '_chem_comp_bond.comp_id': ['VWW', 'VWW', 'VWW', 'VWW'],
        '_chem_comp_bond.atom_id_1': ['C1', 'C1', 'C2', 'N3'],
        '_chem_comp_bond.atom_id_2': ['C2', 'N3', 'O4', 'S5'],
        '_chem_comp_bond.value_order': ['SING', 'DOUB', 'SING', 'ARO'],
        '_chem_comp_bond.pdbx_aromatic_flag': ['N', 'Y', '?', '?'],
        '_chem_comp_bond.pdbx_stereo_config': ['N', 'Y', '?', '?'],
    })

    self.assertEqual(
        exp_chem_comp_data.to_mmcif_dict(),
        {k: v for k, v in chem_comp_cif_dict.items()},
    )

    # with 2 decimal places we lose precision, so we can't compare directly.
    chem_comp_data = ChemicalComponentsData.from_mmcif(
        chem_comp_cif_dict, fix_mse=False, fix_unknown_dna=False
    )
    self.assertSameStructure(chem_comp_data, exp_chem_comp_data)

    # Round trip
    self.assertSameStructure(
        exp_chem_comp_data,
        ChemicalComponentsData.from_mmcif(
            mmcif.Mmcif({**exp_chem_comp_data.to_mmcif_dict()}),
            fix_mse=False,
            fix_unknown_dna=False,
        ),
    )

  def test_chem_comp_atoms_with_missing(self):
    vww_entry = ChemCompEntry(
        type='peptide-like',
        name='TEST',
        pdbx_synonyms='?',
        formula='C23 H27 N3 O6 S',
        formula_weight='473.542',
        mon_nstd_flag='y',
        pdbx_smiles=None,
        chem_comp_atoms={
            'C1': ChemCompAtom(
                type_symbol='C',
                charge=0,
                leaving_atom_flag='N',
                model_ideal_x=1.0,
                model_ideal_y=2.0,
                model_ideal_z=3.0,
                ordinal=1,
            ),
            'C2': ChemCompAtom(
                type_symbol='C',
                charge=1,
                leaving_atom_flag='N',
                model_ideal_x=None,
                model_ideal_y=None,
                model_ideal_z=None,
                ordinal=2,
            ),
        },
        chem_comp_bonds=[
            ChemCompBond(
                atom_id_1='C1',
                atom_id_2='C2',
                value_order='SING',
                aromatic_flag='N',
                stereo_config='N',
            ),
        ],
    )
    exp_chem_comp_data = ChemicalComponentsData(chem_comp={'VWW': vww_entry})

    chem_comp_cif_dict = mmcif.Mmcif({
        '_chem_comp.id': ['VWW'],
        '_chem_comp.name': ['TEST'],
        '_chem_comp.type': ['peptide-like'],
        '_chem_comp.formula_weight': ['473.542'],
        '_chem_comp.formula': ['C23 H27 N3 O6 S'],
        '_chem_comp.mon_nstd_flag': ['y'],
        '_chem_comp.pdbx_synonyms': ['?'],
        # Atom fields
        '_chem_comp_atom.comp_id': ['VWW', 'VWW'],
        '_chem_comp_atom.atom_id': ['C1', 'C2'],
        '_chem_comp_atom.pdbx_ordinal': ['1', '2'],
        '_chem_comp_atom.type_symbol': ['C', 'C'],
        '_chem_comp_atom.charge': ['0', '1'],
        '_chem_comp_atom.pdbx_leaving_atom_flag': ['N', 'N'],
        '_chem_comp_atom.pdbx_model_Cartn_x_ideal': [
            '1.000',
            '?',
        ],
        '_chem_comp_atom.pdbx_model_Cartn_y_ideal': [
            '2.000',
            '?',
        ],
        '_chem_comp_atom.pdbx_model_Cartn_z_ideal': [
            '3.000',
            '?',
        ],
        # Bond fields
        '_chem_comp_bond.comp_id': ['VWW'],
        '_chem_comp_bond.atom_id_1': ['C1'],
        '_chem_comp_bond.atom_id_2': ['C2'],
        '_chem_comp_bond.value_order': ['SING'],
        '_chem_comp_bond.pdbx_aromatic_flag': ['N'],
        '_chem_comp_bond.pdbx_stereo_config': ['N'],
    })

    self.assertEqual(
        exp_chem_comp_data.to_mmcif_dict(),
        {k: v for k, v in chem_comp_cif_dict.items()},
    )

    # with 2 decimal places we lose precision, so we can't compare directly.
    chem_comp_data = ChemicalComponentsData.from_mmcif(
        chem_comp_cif_dict, fix_mse=False, fix_unknown_dna=False
    )
    self.assertSameStructure(chem_comp_data, exp_chem_comp_data)

    # Round trip
    self.assertSameStructure(
        exp_chem_comp_data,
        ChemicalComponentsData.from_mmcif(
            mmcif.Mmcif({**exp_chem_comp_data.to_mmcif_dict()}),
            fix_mse=False,
            fix_unknown_dna=False,
        ),
    )

  def test_chem_comp_atoms_from_mmcif(self):
    mmcif_2pv7 = mmcif.from_string(test_utils.read_test_cif_file('2pv7'))
    chem_comp_data = ChemicalComponentsData.from_mmcif(
        mmcif_2pv7, fix_mse=False, fix_unknown_dna=False
    )
    chem_comp_data_dict = chem_comp_data.to_mmcif_dict()
    self.assertContainsSubset(chem_comp_data_dict.keys(), mmcif_2pv7.keys())
    for key, value in chem_comp_data_dict.items():
      self.assertEqual(value, mmcif_2pv7[key], msg=key)

  def test_chem_comp_atoms_decimal(self):
    vww_entry = ChemCompEntry(
        type='peptide-like',
        name='TEST',
        pdbx_synonyms='?',
        formula='C23 H27 N3 O6 S',
        formula_weight='473.542',
        mon_nstd_flag='y',
        pdbx_smiles=None,
        chem_comp_atoms={
            'C1': ChemCompAtom(
                type_symbol='C',
                ordinal=1,
                charge=0,
                leaving_atom_flag='N',
                model_ideal_x=1 / 3,
                model_ideal_y=2.0,
                model_ideal_z=3.0,
            )
        },
    )
    chem_comp_data = ChemicalComponentsData(chem_comp={'VWW': vww_entry})

    exp_chem_comp = mmcif.Mmcif({
        '_chem_comp.id': ['VWW'],
        '_chem_comp.name': ['TEST'],
        '_chem_comp.type': ['peptide-like'],
        '_chem_comp.formula_weight': ['473.542'],
        '_chem_comp.formula': ['C23 H27 N3 O6 S'],
        '_chem_comp.mon_nstd_flag': ['y'],
        '_chem_comp.pdbx_synonyms': ['?'],
        '_chem_comp_atom.comp_id': ['VWW'],
        '_chem_comp_atom.atom_id': ['C1'],
        '_chem_comp_atom.type_symbol': ['C'],
        '_chem_comp_atom.charge': ['0'],
        '_chem_comp_atom.pdbx_ordinal': ['1'],
        '_chem_comp_atom.pdbx_leaving_atom_flag': ['N'],
        '_chem_comp_atom.pdbx_model_Cartn_x_ideal': [
            '0.333',
        ],
        '_chem_comp_atom.pdbx_model_Cartn_y_ideal': [
            '2.000',
        ],
        '_chem_comp_atom.pdbx_model_Cartn_z_ideal': [
            '3.000',
        ],
    })
    self.assertEqual(
        chem_comp_data.to_mmcif_dict(),
        {k: v for k, v in exp_chem_comp.items()},
    )

  def test_chem_comp_atoms_missing_fields(self):
    vww_entry = ChemCompEntry(
        type='peptide-like',
        name='TEST',
        pdbx_synonyms='?',
        formula='C23 H27 N3 O6 S',
        formula_weight='473.542',
        mon_nstd_flag='y',
        pdbx_smiles=None,
        chem_comp_atoms={
            'C1': ChemCompAtom(
                type_symbol='C',
                leaving_atom_flag='?',
                ordinal=1,
            ),
        },
        chem_comp_bonds=[
            ChemCompBond(
                atom_id_1='C1',
                atom_id_2='C2',
                value_order='?',
                aromatic_flag='?',
                stereo_config='?',
            ),
        ],
    )
    exp_chem_comp_data = ChemicalComponentsData(chem_comp={'VWW': vww_entry})

    chem_comp_cif_dict = mmcif.Mmcif({
        '_chem_comp.id': ['VWW'],
        '_chem_comp.name': ['TEST'],
        '_chem_comp.type': ['peptide-like'],
        '_chem_comp.formula_weight': ['473.542'],
        '_chem_comp.formula': ['C23 H27 N3 O6 S'],
        '_chem_comp.mon_nstd_flag': ['y'],
        '_chem_comp.pdbx_synonyms': ['?'],
        # Atom fields
        '_chem_comp_atom.comp_id': ['VWW'],
        '_chem_comp_atom.atom_id': ['C1'],
        '_chem_comp_atom.type_symbol': ['C'],
        '_chem_comp_atom.pdbx_ordinal': ['1'],
        # Bond fields
        '_chem_comp_bond.comp_id': ['VWW'],
        '_chem_comp_bond.atom_id_1': ['C1'],
        '_chem_comp_bond.atom_id_2': ['C2'],
    })

    chem_comp_data = ChemicalComponentsData.from_mmcif(
        chem_comp_cif_dict, fix_mse=False, fix_unknown_dna=False
    )
    self.assertEqual(chem_comp_data, exp_chem_comp_data)


class ChemicalComponentsTest(parameterized.TestCase):

  def test_chem_comp_entry_empty_field(self):
    with self.assertRaisesRegex(ValueError, 'empty string'):
      _ = ChemCompEntry(type='non-polymer', pdbx_smiles='')

  @parameterized.parameters([
      (
          ChemCompEntry(type='non-polymer', name='specified'),
          ChemCompEntry(type='non-polymer', name='specified'),
          True,
      ),
      (
          ChemCompEntry(type='non-polymer', name='specified'),
          ChemCompEntry(type='non-polymer', name='.'),
          True,
      ),
      (
          ChemCompEntry(type='non-polymer', name='specified'),
          ChemCompEntry(type='non-polymer'),
          True,
      ),
      (
          ChemCompEntry(type='non-polymer'),
          ChemCompEntry(type='non-polymer', name='specified'),
          False,
      ),
      (
          ChemCompEntry(type='non-polymer', name='value'),
          ChemCompEntry(type='non-polymer', name='another value'),
          False,
      ),
  ])
  def test_chem_comp_entry_extends(self, entry, other_entry, expected_result):
    self.assertEqual(entry.extends(other_entry), expected_result)

  def test_get_data_for_ccd_components(self):
    chem_comp = struc_chem_comps.get_data_for_ccd_components(
        ccd=chemical_components.Ccd(),
        chemical_component_ids=['VWW', 'ALN', 'VWY'],
    ).chem_comp

    self.assertEqual(
        chem_comp,
        {
            'ALN': ChemCompEntry(
                type='L-PEPTIDE LINKING',
                name='NAPHTHALEN-2-YL-3-ALANINE',
                pdbx_synonyms='?',
                formula='C13 H13 N O2',
                formula_weight='215.248',
                mon_nstd_flag='n',
            ),
            'VWW': ChemCompEntry(
                type='peptide-like',
                name='L-gamma-glutamyl-S-benzyl-N-[(S)-carboxy(phenyl)methyl]-L-cysteinamide',
                pdbx_synonyms='?',
                formula='C23 H27 N3 O6 S',
                formula_weight='473.542',
                mon_nstd_flag='y',
            ),
            'VWY': ChemCompEntry(
                type='non-polymer',
                name='N-ethyl-4-[(methylsulfonyl)amino]benzamide',
                pdbx_synonyms='?',
                formula='C10 H14 N2 O3 S',
                formula_weight='242.295',
                mon_nstd_flag='.',
            ),
        },
    )

  def test_get_data_for_ccd_components_with_smiles(self):
    chem_comp = struc_chem_comps.get_data_for_ccd_components(
        ccd=chemical_components.Ccd(),
        chemical_component_ids=['VWY'],
        populate_pdbx_smiles=True,
    ).chem_comp

    self.assertEqual(
        chem_comp,
        {
            'VWY': ChemCompEntry(
                type='non-polymer',
                name='N-ethyl-4-[(methylsulfonyl)amino]benzamide',
                pdbx_synonyms='?',
                formula='C10 H14 N2 O3 S',
                formula_weight='242.295',
                mon_nstd_flag='.',
                pdbx_smiles='CCNC(=O)c1ccc(cc1)NS(=O)(=O)C',
            ),
        },
    )

  def test_get_data_for_ccd_components_ignore_unknowns(self):
    self.assertEmpty(
        struc_chem_comps.get_data_for_ccd_components(
            ccd=chemical_components.Ccd(),
            chemical_component_ids=['UNKNOWN_LIGAND'],
        ).chem_comp
    )

  @parameterized.parameters([True, False])
  def test_populate_missing_ccd_data(self, populate_pdbx_smiles):
    chemical_components_data = ChemicalComponentsData(
        chem_comp={
            'ALN': ChemCompEntry(
                type='L-PEPTIDE LINKING',
                name='Alanine',  # value doesn't match CCD content
                pdbx_smiles='custom value',
            ),
            'DM-X': ChemCompEntry(type='non-polymer', name='Secret'),
        }
    )

    updated_chem_comp = struc_chem_comps.populate_missing_ccd_data(
        ccd=chemical_components.Ccd(),
        chemical_components_data=chemical_components_data,
        chemical_component_ids=['ALN'],
        populate_pdbx_smiles=populate_pdbx_smiles,
    ).chem_comp
    self.assertEqual(
        updated_chem_comp,
        {
            'ALN': ChemCompEntry(
                type='L-PEPTIDE LINKING',
                name='Alanine',  # value doesn't match CCD content
                pdbx_synonyms='?',
                formula='C13 H13 N O2',
                formula_weight='215.248',
                mon_nstd_flag='n',
                pdbx_smiles='custom value',
            ),
            'DM-X': ChemCompEntry(type='non-polymer', name='Secret'),
        },
    )

    # Doesn't fail for non-CCD component id, no changes to initial content.
    self.assertEqual(
        struc_chem_comps.populate_missing_ccd_data(
            ccd=chemical_components.Ccd(),
            chemical_components_data=chemical_components_data,
            chemical_component_ids=['DM-X'],
            populate_pdbx_smiles=populate_pdbx_smiles,
        ),
        chemical_components_data,
    )

  @parameterized.parameters(
      [(True, 'c1ccc2c(c1)cccc2C[C@@H](C(=O)O)N'), (False, None)]
  )
  def test_populate_missing_ccd_data_no_pdbx_smiles_set(
      self, populate_pdbx_smiles, expected_smiles
  ):
    chemical_components_data = ChemicalComponentsData(
        chem_comp={
            'ALN': ChemCompEntry(
                type='L-PEPTIDE LINKING',
                name='Alanine',  # Value doesn't match CCD content.
            ),
            'DM-X': ChemCompEntry(type='non-polymer', name='Secret'),
        }
    )

    updated_chem_comp = struc_chem_comps.populate_missing_ccd_data(
        ccd=chemical_components.Ccd(),
        chemical_components_data=chemical_components_data,
        chemical_component_ids=['ALN'],
        populate_pdbx_smiles=populate_pdbx_smiles,
    ).chem_comp

    # It's a long list of atoms/bonds to compare, so compare just the length.
    self.assertEqual(
        updated_chem_comp,
        {
            'ALN': ChemCompEntry(
                type='L-PEPTIDE LINKING',
                name='Alanine',  # Value doesn't match CCD content.
                pdbx_synonyms='?',
                formula='C13 H13 N O2',
                formula_weight='215.248',
                mon_nstd_flag='n',
                pdbx_smiles=expected_smiles,
            ),
            'DM-X': ChemCompEntry(type='non-polymer', name='Secret'),
        },
    )

    # Doesn't fail for non-CCD component id, no changes to initial content.
    self.assertEqual(
        struc_chem_comps.populate_missing_ccd_data(
            ccd=chemical_components.Ccd(),
            chemical_components_data=chemical_components_data,
            chemical_component_ids=['DM-X'],
            populate_pdbx_smiles=populate_pdbx_smiles,
        ),
        chemical_components_data,
    )

  @parameterized.parameters([
      ('GLY', 11, 10),
      ('PRO', 19, 19),
      ('ARG', 28, 27),
      ('TRP', 28, 29),
      ('UNK', 11, 10),
      ('A', 37, 39),
      ('G', 38, 40),
      ('C', 35, 36),
      ('U', 34, 35),
      ('DA', 36, 38),
      ('DG', 37, 39),
      ('DC', 34, 35),
      ('DT', 36, 37),
  ])
  def test_get_all_atoms_in_entry(self, res_name, exp_atoms, exp_bonds):
    res = struc_chem_comps.get_all_atoms_in_entry(
        ccd=chemical_components.Ccd(), res_name=res_name
    )
    self.assertLen(res['_chem_comp_atom.atom_id'], exp_atoms)
    self.assertLen(res['_chem_comp_atom.type_symbol'], exp_atoms)
    self.assertLen(res['_chem_comp_bond.atom_id_1'], exp_bonds)
    self.assertLen(res['_chem_comp_bond.atom_id_2'], exp_bonds)

  @parameterized.parameters([
      (
          'GLY',
          {'N', 'CA', 'C', 'O', 'OXT', 'H', 'H2', 'HA2', 'HA3', 'HXT', 'H3'},
      ),
      ('ZN', {'ZN'}),
      ('UNL', set()),
  ])
  def test_get_res_atom_names(self, res_name, expected_atoms):
    atoms = struc_chem_comps.get_res_atom_names(
        ccd=chemical_components.Ccd(), res_name=res_name
    )
    self.assertEqual(expected_atoms, atoms)

  def test_get_res_atom_names_error(self):
    with self.assertRaisesRegex(ValueError, 'Unknown residue type GGGG'):
      struc_chem_comps.get_res_atom_names(
          ccd=chemical_components.Ccd(), res_name='GGGG'
      )


if __name__ == '__main__':
  absltest.main()
