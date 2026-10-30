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

import os

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3 import structure
from alphafold3.common import resources
from alphafold3.constants import chemical_components
from alphafold3.constants import mmcif_names
from alphafold3.model.pipeline import inter_chain_bonds
from alphafold3.model.pipeline import structure_cleaning
from alphafold3.structure import test_utils
import numpy as np

_MINIATURE_PDB_MMCIF = (
    resources.ROOT / 'test_data/miniature_databases/pdb_mmcif'
)


def _read_test_cif_file(
    filename: str, data_dir: str | os.PathLike[str] = _MINIATURE_PDB_MMCIF
) -> structure.Structure:
  data_path = os.path.join(data_dir, f'{filename}.cif')
  with resources.open_resource(data_path, 'rt') as f:
    return structure.from_mmcif(f.read(), include_bonds=True)


class StructureCleaningTest(test_utils.StructureTestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.input_struc = structure.from_atom_arrays(
        chain_id=np.array(['A', 'A', 'B', 'B', 'B', 'B'], dtype=object),
        chain_type=np.array(
            [
                mmcif_names.PROTEIN_CHAIN,
                mmcif_names.PROTEIN_CHAIN,
                mmcif_names.BRANCHED_CHAIN,
                mmcif_names.BRANCHED_CHAIN,
                mmcif_names.BRANCHED_CHAIN,
                mmcif_names.BRANCHED_CHAIN,
            ],
            dtype=object,
        ),
        res_id=np.array([1, 1, 1, 1, 2, 2], dtype=np.int32),
        res_name=np.array(
            ['ASN', 'ASN', 'NAG', 'NAG', 'NAG', 'NAG'], dtype=object
        ),
        atom_name=np.array(
            ['ND2', 'OXT', 'C1', 'O1', 'C1', 'O1'], dtype=object
        ),
        atom_element=np.array(['N', 'O', 'C', 'O', 'C', 'O'], dtype=object),
        atom_x=np.zeros(6, dtype=np.float32),
        atom_y=np.zeros(6, dtype=np.float32),
        atom_z=np.zeros(6, dtype=np.float32),
    )

  @parameterized.parameters(
      ('6ww6', _MINIATURE_PDB_MMCIF),
  )
  def test_clean_structure(self, target_name, data_dir):
    input_struc = _read_test_cif_file(target_name, data_dir=data_dir)
    clean_struc = structure_cleaning.clean_structure(
        input_struc,
        ccd=chemical_components.Ccd(),
        drop_missing_sequence=True,
        drop_non_standard_atoms=True,
        filter_waters=True,
        filter_hydrogens=True,
        filter_leaving_atoms=True,
        only_glycan_ligands_for_leaving_atoms=True,
        covalent_bonds_only=False,
        remove_polymer_polymer_bonds=False,
        remove_bad_bonds=False,
        fix_standalone_glycans=False,
    )
    self.assertCountEqual(
        clean_struc.author_naming_scheme.auth_asym_id.keys(), clean_struc.chains
    )
    self.assertCountEqual(
        clean_struc.author_naming_scheme.auth_seq_id.keys(), clean_struc.chains
    )
    self.assertCountEqual(
        clean_struc.author_naming_scheme.insertion_code.keys(),
        clean_struc.chains,
    )
    self.assertCountEqual(clean_struc.all_residues.keys(), clean_struc.chains)

    with self.subTest(name='CleanedAndRecleanedEqual'):
      # No cleaning missing or clashing chains.
      cleaned_clean_struc = structure_cleaning.clean_structure(
          clean_struc,
          ccd=chemical_components.Ccd(),
          drop_missing_sequence=False,
          drop_non_standard_atoms=False,
          filter_waters=True,
          filter_hydrogens=True,
          filter_leaving_atoms=True,
          only_glycan_ligands_for_leaving_atoms=True,
          covalent_bonds_only=False,
          remove_polymer_polymer_bonds=False,
          remove_bad_bonds=False,
          fix_standalone_glycans=False,
      )
      self.assertStructuresEqual(clean_struc, cleaned_clean_struc)

  def test_create_empty_output_struc_and_layout(self):
    input_struc = _read_test_cif_file('6ww6')
    clean_struc = structure_cleaning.clean_structure(
        input_struc,
        ccd=chemical_components.Ccd(),
        drop_missing_sequence=True,
        drop_non_standard_atoms=True,
        filter_waters=True,
        filter_hydrogens=True,
        filter_leaving_atoms=True,
        only_glycan_ligands_for_leaving_atoms=True,
        covalent_bonds_only=False,
        remove_polymer_polymer_bonds=False,
        remove_bad_bonds=False,
        fix_standalone_glycans=False,
    )
    empty_struc, _ = structure_cleaning.create_empty_output_struc_and_layout(
        clean_struc, ccd=chemical_components.Ccd()
    )
    self.assertEqual(empty_struc.coords.sum(), 0)
    self.assertEqual(empty_struc.atom_x.sum(), 0)
    self.assertEqual(empty_struc.atom_y.sum(), 0)
    self.assertEqual(empty_struc.atom_z.sum(), 0)

  def test_filter_leaving_atoms_slow(self):
    input_struc = self.input_struc.add_bonds(
        bonded_atom_pairs=[
            (('A', 1, 'ND2'), ('B', 1, 'C1')),
            (('B', 1, 'O1'), ('B', 2, 'C1')),
        ],
        bond_type=mmcif_names.COVALENT_BOND,
    )

    clean_struc = structure_cleaning.clean_structure(
        input_struc,
        ccd=chemical_components.Ccd(),
        drop_missing_sequence=False,
        drop_non_standard_atoms=False,
        filter_waters=False,
        filter_hydrogens=False,
        filter_leaving_atoms=True,
        only_glycan_ligands_for_leaving_atoms=False,
        covalent_bonds_only=False,
        remove_polymer_polymer_bonds=False,
        remove_bad_bonds=False,
        fix_standalone_glycans=False,
    )

    np.testing.assert_array_equal(
        clean_struc.chain_id, ['A', 'A', 'B', 'B', 'B']
    )
    np.testing.assert_array_equal(clean_struc.res_id, [1, 1, 1, 1, 2])
    np.testing.assert_array_equal(
        clean_struc.res_name, ['ASN', 'ASN', 'NAG', 'NAG', 'NAG']
    )
    np.testing.assert_array_equal(
        clean_struc.atom_name, ['ND2', 'OXT', 'C1', 'O1', 'C1']
    )

  def test_filter_leaving_atoms_fast(self):
    input_struc = self.input_struc.add_bonds(
        bonded_atom_pairs=[
            (('A', 1, 'ND2'), ('B', 1, 'C1')),
            (('B', 1, 'C1'), ('B', 2, 'C1')),
        ],
        bond_type=mmcif_names.COVALENT_BOND,
    )
    clean_struc = structure_cleaning.clean_structure(
        input_struc,
        ccd=chemical_components.Ccd(),
        drop_missing_sequence=False,
        drop_non_standard_atoms=False,
        filter_waters=False,
        filter_hydrogens=False,
        filter_leaving_atoms=True,
        only_glycan_ligands_for_leaving_atoms=True,
        covalent_bonds_only=False,
        remove_polymer_polymer_bonds=False,
        remove_bad_bonds=False,
        fix_standalone_glycans=False,
    )
    np.testing.assert_array_equal(clean_struc.chain_id, ['A', 'A', 'B', 'B'])
    np.testing.assert_array_equal(clean_struc.res_id, [1, 1, 1, 2])
    np.testing.assert_array_equal(
        clean_struc.res_name, ['ASN', 'ASN', 'NAG', 'NAG']
    )
    np.testing.assert_array_equal(
        clean_struc.atom_name, ['ND2', 'OXT', 'C1', 'C1']
    )

  def test_filter_leaving_atoms_branched(self):
    # Create a structure with 5 NAG residues, each with a C1 and O1 atom with
    # a branched structure:
    # 1-2-3
    #   |
    #   4-5
    input_struc = structure.from_atom_arrays(
        chain_id=np.array(['A'] * 10, dtype=object),
        chain_type=np.array([mmcif_names.BRANCHED_CHAIN] * 10, dtype=object),
        res_id=np.array([1, 1, 2, 2, 3, 3, 4, 4, 5, 5], dtype=np.int32),
        res_name=np.array(['NAG'] * 10, dtype=object),
        atom_name=np.array(['C1', 'O1'] * 5, dtype=object),
        atom_element=np.array(['C', 'O'] * 5, dtype=object),
        atom_x=np.zeros(10, dtype=np.float32),
        atom_y=np.zeros(10, dtype=np.float32),
        atom_z=np.zeros(10, dtype=np.float32),
    ).add_bonds(
        bonded_atom_pairs=[
            (('A', 1, 'O1'), ('A', 2, 'C1')),  # Keeps O1 for residue 1.
            (('A', 2, 'C1'), ('A', 3, 'C1')),  # No effect - C1 always kept.
            (('A', 2, 'O1'), ('A', 4, 'C1')),  # Keeps O1 for residue 2.
            # This final bond is modified below to have type "?" which means
            # that it won't cause the cleaner to keep hold of this O1 atom.
            (('A', 4, 'O1'), ('A', 5, 'C1')),  # Not kept because of "?" type.
        ],
        bond_type=mmcif_names.COVALENT_BOND,
    )

    # Set the final bond to have type "?".
    new_bond_type = input_struc.bonds.type.copy()
    new_bond_type[-1] = '?'
    input_struc = input_struc.copy_and_update(
        bonds=input_struc.bonds.copy_and_update(type=new_bond_type)
    )

    clean_struc = structure_cleaning.clean_structure(
        input_struc,
        ccd=chemical_components.Ccd(),
        drop_missing_sequence=False,
        drop_non_standard_atoms=False,
        filter_waters=False,
        filter_hydrogens=False,
        filter_leaving_atoms=True,
        only_glycan_ligands_for_leaving_atoms=True,
        covalent_bonds_only=False,
        remove_polymer_polymer_bonds=False,
        remove_bad_bonds=False,
        fix_standalone_glycans=False,
    )
    np.testing.assert_array_equal(clean_struc.res_id, [1, 1, 2, 2, 3, 4, 5])
    np.testing.assert_array_equal(
        clean_struc.atom_name, ['C1', 'O1', 'C1', 'O1', 'C1', 'C1', 'C1']
    )

  @parameterized.parameters([
      ('5s9m', _MINIATURE_PDB_MMCIF),  # One chain with bonds.
  ])
  def test_remove_nonsymmetric_bonds_on_symmetrical_structure(
      self, mmcif_id, data_dir
  ):
    struc = _read_test_cif_file(mmcif_id, data_dir=data_dir)
    struc = struc.generate_bioassembly()
    cleaned_struc = structure_cleaning.clean_structure(
        struc,
        ccd=chemical_components.Ccd(),
        drop_non_standard_atoms=False,
        drop_missing_sequence=False,
        filter_waters=False,
        filter_hydrogens=False,
        filter_leaving_atoms=False,
        only_glycan_ligands_for_leaving_atoms=False,
        covalent_bonds_only=False,
        remove_polymer_polymer_bonds=False,
        remove_bad_bonds=False,
        fix_standalone_glycans=False,
    )
    # Symmetric bond patterns should result in no filtered bonds.
    self.assertEqual(len(cleaned_struc.bonds), len(struc.bonds))

  @parameterized.parameters([
      ('5s9m', _MINIATURE_PDB_MMCIF, 9),
  ])
  def test_covalent_bond_filter(
      self, mmcif_id, data_dir, expected_num_covalent_bonds
  ):
    struc = _read_test_cif_file(mmcif_id, data_dir=data_dir)
    struc = struc.generate_bioassembly()
    cleaned_struc = structure_cleaning.clean_structure(
        struc,
        ccd=chemical_components.Ccd(),
        drop_non_standard_atoms=False,
        drop_missing_sequence=False,
        filter_waters=False,
        filter_hydrogens=False,
        filter_leaving_atoms=False,
        only_glycan_ligands_for_leaving_atoms=False,
        covalent_bonds_only=True,
        remove_polymer_polymer_bonds=False,
        remove_bad_bonds=False,
        fix_standalone_glycans=False,
    )
    self.assertLess(len(cleaned_struc.bonds), len(struc.bonds))
    self.assertLen(cleaned_struc.bonds, expected_num_covalent_bonds)

  @parameterized.parameters([
      ('170d', _MINIATURE_PDB_MMCIF, 0),  # Nucleic covalent bonds.
      ('1al4', _MINIATURE_PDB_MMCIF, 0),  # Protein covalent bonds.
  ])
  def test_polymer_polymer_bond_filter(
      self, mmcif_id, data_dir, expected_num_bonds
  ):
    struc = _read_test_cif_file(mmcif_id, data_dir=data_dir)
    struc = struc.generate_bioassembly()

    # Remove non-covalent bonds, not testing covalent bond filtering here.
    is_covalent = np.isin(struc.bonds.type, ['covale'])
    struc_covale_only = struc.copy_and_update(bonds=struc.bonds[is_covalent])  # pyrefly: ignore[bad-argument-type]

    cleaned_struc = structure_cleaning.clean_structure(
        struc_covale_only,
        ccd=chemical_components.Ccd(),
        drop_non_standard_atoms=False,
        drop_missing_sequence=False,
        filter_waters=False,
        filter_hydrogens=False,
        filter_leaving_atoms=False,
        only_glycan_ligands_for_leaving_atoms=False,
        covalent_bonds_only=False,
        remove_polymer_polymer_bonds=True,
        remove_bad_bonds=False,
        fix_standalone_glycans=False,
    )
    num_bonds = len(cleaned_struc.bonds)
    self.assertEqual(num_bonds, expected_num_bonds)

  def test_bad_bond_filter(self):
    struc = _read_test_cif_file('1h9h', data_dir=_MINIATURE_PDB_MMCIF)
    struc = struc.generate_bioassembly()
    struc = struc.copy_and_update(bonds=structure.Bonds.make_empty())
    self.assertEmpty(struc.bonds)

    # Add in fake bad bond and check it is filtered out later.
    struc = struc.add_bonds(
        ((('A', 1, 'C'), ('A', 5, 'N')),), bond_type=mmcif_names.COVALENT_BOND
    )
    self.assertLen(struc.bonds, 1)

    cleaned_struc = structure_cleaning.clean_structure(
        struc,
        ccd=chemical_components.Ccd(),
        drop_non_standard_atoms=False,
        drop_missing_sequence=False,
        filter_waters=False,
        filter_hydrogens=False,
        filter_leaving_atoms=False,
        only_glycan_ligands_for_leaving_atoms=False,
        covalent_bonds_only=False,
        remove_polymer_polymer_bonds=False,
        remove_bad_bonds=True,
        fix_standalone_glycans=False,
    )
    self.assertEmpty(cleaned_struc.bonds)

  @parameterized.parameters([
      # 7pxi - has one water bond that is always removed.
      ('7pxi', _MINIATURE_PDB_MMCIF, True, True, 1, 1, 0),
  ])
  def test_combined_bond_filters(
      self,
      mmcif_id,
      data_dir,
      covalent_bonds_only,
      remove_polymer_polymer_bonds,
      expected_num_bonds,
      expected_num_covale_bonds,
      expected_num_polymer_polymer_bonds,
  ):
    struc = _read_test_cif_file(mmcif_id, data_dir=data_dir)
    struc = structure_cleaning.clean_structure(
        struc,
        ccd=chemical_components.Ccd(),
        drop_non_standard_atoms=False,
        drop_missing_sequence=False,
        filter_waters=False,
        filter_hydrogens=False,
        filter_leaving_atoms=False,
        only_glycan_ligands_for_leaving_atoms=False,
        covalent_bonds_only=covalent_bonds_only,
        remove_polymer_polymer_bonds=remove_polymer_polymer_bonds,
        remove_bad_bonds=True,
        fix_standalone_glycans=False,
    )
    num_bonds = len(struc.bonds)
    self.assertEqual(num_bonds, expected_num_bonds)

    num_covale_bonds = sum(np.isin(struc.bonds.type, ['covale']))
    self.assertEqual(num_covale_bonds, expected_num_covale_bonds)

    num_pp_bonds = 0
    for bond in struc.iter_bonds():
      if (
          bond.from_atom['chain_type'] in mmcif_names.POLYMER_CHAIN_TYPES
          and bond.dest_atom['chain_type'] in mmcif_names.POLYMER_CHAIN_TYPES
      ):
        num_pp_bonds += 1
    self.assertEqual(num_pp_bonds, expected_num_polymer_polymer_bonds)

  def test_fix_standalone_glycans(self):
    struc = structure.from_sequences_and_bonds(
        sequences=[
            '(ASN)(ALA)(ALA)(ALA)(ASN)',  # A = Protein chain, 5 residues.
            '(NAG)',  # B = Standalone glycan.
            '(NAG)',  # C = Linked glycan.
            '(NAG)(FUC)',  # D = Standalone branched glycan.
            '(NAG)(FUC)',  # E = Linked branched glycan.
        ],
        chain_types=(
            [mmcif_names.PROTEIN_CHAIN] + [mmcif_names.NON_POLYMER_CHAIN] * 4
        ),
        sequence_formats=[structure.SequenceFormat.CCD_CODES] * 5,
        chain_ids=['A', 'B', 'C', 'D', 'E'],
        bonded_atom_pairs=[
            ((0, 0, 'ND2'), (2, 0, 'C1')),  # A--C: ASN--NAG.
            ((3, 0, 'O6'), (3, 1, 'C1')),  # D--D: NAG--FUC.
            ((0, 0, 'ND2'), (4, 0, 'C1')),  # A--E: ASN--NAG.
            ((4, 0, 'O6'), (4, 1, 'C1')),  # E--E: NAG--FUC.
        ],
        ccd=chemical_components.Ccd(),
        name='with_glycans',
        bond_type=mmcif_names.COVALENT_BOND,
        release_date=None,
    )
    polymer_ligand_bonds, ligand_ligand_bonds = (
        inter_chain_bonds.get_polymer_ligand_and_ligand_ligand_bonds(
            struc, only_glycan_ligands=False, allow_multiple_bonds_per_atom=True
        )
    )
    all_atom_struc, _ = structure_cleaning.create_empty_output_struc_and_layout(
        struc=struc,
        ccd=chemical_components.Ccd(),
        polymer_ligand_bonds=polymer_ligand_bonds,
        ligand_ligand_bonds=ligand_ligand_bonds,
        drop_ligand_leaving_atoms=False,  # Do not do any corrections yet.
        fix_standalone_glycans=False,  # Do not do any corrections yet.
    )

    clean_struc_unfixed = structure_cleaning.clean_structure(
        all_atom_struc,
        ccd=chemical_components.Ccd(),
        drop_non_standard_atoms=False,
        drop_missing_sequence=False,
        filter_waters=False,
        filter_hydrogens=False,
        filter_leaving_atoms=True,
        only_glycan_ligands_for_leaving_atoms=True,
        covalent_bonds_only=False,
        remove_polymer_polymer_bonds=False,
        remove_bad_bonds=False,
        fix_standalone_glycans=False,
    )
    clean_struc_fixed = structure_cleaning.clean_structure(
        all_atom_struc,
        ccd=chemical_components.Ccd(),
        drop_non_standard_atoms=False,
        drop_missing_sequence=False,
        filter_waters=False,
        filter_hydrogens=False,
        filter_leaving_atoms=True,
        only_glycan_ligands_for_leaving_atoms=True,
        covalent_bonds_only=False,
        remove_polymer_polymer_bonds=False,
        remove_bad_bonds=False,
        fix_standalone_glycans=True,
    )

    with self.subTest('all_atom_struc_atoms'):
      # 4 O1s in NAGs, 2 O1s in FUCs.
      self.assertEqual(all_atom_struc.num_atoms, 114)
      self.assertEqual(all_atom_struc.filter(atom_name='O1').num_atoms, 6)

    with self.subTest('clean_struc_unfixed_atoms'):
      # All 6 O1s are (incorrectly) removed.
      self.assertEqual(clean_struc_unfixed.num_atoms, 108)
      self.assertEqual(clean_struc_unfixed.filter(atom_name='O1').num_atoms, 0)

    with self.subTest('clean_struc_fixed_atoms'):
      # Only 5 O1s are removed (kept for the standalone NAG).
      self.assertEqual(clean_struc_fixed.num_atoms, 109)
      self.assertEqual(clean_struc_fixed.filter(atom_name='O1').num_atoms, 1)


if __name__ == '__main__':
  absltest.main()
