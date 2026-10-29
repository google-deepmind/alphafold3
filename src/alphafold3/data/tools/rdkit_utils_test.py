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

import string

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.constants import chemical_components
from alphafold3.cpp import cif_dict
from alphafold3.data.tools import rdkit_utils
import rdkit.Chem as rd_chem
from rdkit.Chem import AllChem


# pylint: disable=trailing-whitespace
# fmt: off
_CCD_PDB = """HETATM    1  C   VWW     1       4.213   5.060   0.877  1.00  0.00           C  
HETATM    2  C1  VWW     1      -1.624   0.494  -0.557  1.00  0.00           C  
HETATM    3  C2  VWW     1      -4.564   1.545  -1.122  1.00  0.00           C  
HETATM    4  CA  VWW     1       3.852   4.147  -0.266  1.00  0.00           C  
HETATM    5  CA1 VWW     1      -0.460  -0.314  -1.071  1.00  0.00           C  
HETATM    6  CA2 VWW     1      -3.984   0.728   0.003  1.00  0.00           C  
HETATM    7  CB  VWW     1       2.729   3.203   0.170  1.00  0.00           C  
HETATM    8  CB1 VWW     1       0.050  -1.236   0.038  1.00  0.00           C  
HETATM    9  CB2 VWW     1      -5.039  -0.201   0.545  1.00  0.00           C  
HETATM   10  CD  VWW     1       1.172   1.432  -0.596  1.00  0.00           C  
HETATM   11  CD1 VWW     1       1.909  -3.221   0.852  1.00  0.00           C  
HETATM   12 CD11 VWW     1      -6.687  -0.726   2.202  1.00  0.00           C  
HETATM   13  CD2 VWW     1      -6.297  -2.232   0.381  1.00  0.00           C  
HETATM   14  CE  VWW     1       3.079  -4.100   0.493  1.00  0.00           C  
HETATM   15  CE1 VWW     1      -6.976  -1.905   1.541  1.00  0.00           C  
HETATM   16  CG  VWW     1       2.278   2.361  -1.025  1.00  0.00           C  
HETATM   17  CG1 VWW     1      -5.718   0.126   1.704  1.00  0.00           C  
HETATM   18  CG2 VWW     1      -5.332  -1.378  -0.119  1.00  0.00           C  
HETATM   19  CH  VWW     1       5.226  -5.713  -0.167  1.00  0.00           C  
HETATM   20  CT1 VWW     1       5.444  -4.454   0.362  1.00  0.00           C  
HETATM   21  CT2 VWW     1       3.935  -6.167  -0.361  1.00  0.00           C  
HETATM   22  CZ1 VWW     1       4.370  -3.648   0.692  1.00  0.00           C  
HETATM   23  CZ2 VWW     1       2.862  -5.361  -0.031  1.00  0.00           C  
HETATM   24  H   VWW     1       5.360   2.799   0.117  1.00  0.00           H  
HETATM   25  H1  VWW     1       0.926   0.582  -2.408  1.00  0.00           H  
HETATM   26  H2  VWW     1       4.832   2.786  -1.460  1.00  0.00           H  
HETATM   27  H3  VWW     1      -2.987  -0.974  -0.784  1.00  0.00           H  
HETATM   28  HA  VWW     1       3.517   4.742  -1.115  1.00  0.00           H  
HETATM   29  HA1 VWW     1      -0.782  -0.913  -1.923  1.00  0.00           H  
HETATM   30  HA2 VWW     1      -3.642   1.393   0.796  1.00  0.00           H  
HETATM   31  HB2 VWW     1       3.093   2.546   0.960  1.00  0.00           H  
HETATM   32 HB21 VWW     1       0.372  -0.636   0.889  1.00  0.00           H  
HETATM   33  HB3 VWW     1       1.888   3.787   0.542  1.00  0.00           H  
HETATM   34 HB31 VWW     1      -0.750  -1.908   0.348  1.00  0.00           H  
HETATM   35  HD1 VWW     1      -7.217  -0.471   3.108  1.00  0.00           H  
HETATM   36  HD2 VWW     1       2.184  -2.572   1.684  1.00  0.00           H  
HETATM   37 HD21 VWW     1      -6.523  -3.153  -0.135  1.00  0.00           H  
HETATM   38  HD3 VWW     1       1.063  -3.844   1.142  1.00  0.00           H  
HETATM   39  HE  VWW     1      -7.732  -2.571   1.930  1.00  0.00           H  
HETATM   40  HG1 VWW     1      -5.492   1.047   2.221  1.00  0.00           H  
HETATM   41  HG2 VWW     1       1.914   3.018  -1.816  1.00  0.00           H  
HETATM   42 HG21 VWW     1      -4.805  -1.631  -1.026  1.00  0.00           H  
HETATM   43  HG3 VWW     1       3.120   1.777  -1.397  1.00  0.00           H  
HETATM   44  HH  VWW     1       6.065  -6.342  -0.425  1.00  0.00           H  
HETATM   45  HT1 VWW     1       6.452  -4.100   0.516  1.00  0.00           H  
HETATM   46  HT2 VWW     1       3.765  -7.152  -0.769  1.00  0.00           H  
HETATM   47  HXT VWW     1       3.727   6.693   1.884  1.00  0.00           H  
HETATM   48 HXT1 VWW     1      -5.954   2.851  -1.649  1.00  0.00           H  
HETATM   49  HZ1 VWW     1       4.540  -2.664   1.105  1.00  0.00           H  
HETATM   50  HZ2 VWW     1       1.853  -5.716  -0.182  1.00  0.00           H  
HETATM   51  N   VWW     1       5.030   3.360  -0.654  1.00  0.00           N  
HETATM   52  N1  VWW     1       0.615   0.590  -1.489  1.00  0.00           N  
HETATM   53  N2  VWW     1      -2.852  -0.057  -0.497  1.00  0.00           N  
HETATM   54  O   VWW     1       5.183   4.823   1.557  1.00  0.00           O  
HETATM   55  O1  VWW     1      -1.455   1.640  -0.198  1.00  0.00           O  
HETATM   56  O2  VWW     1      -4.082   1.478  -2.228  1.00  0.00           O  
HETATM   57  OE1 VWW     1       0.783   1.442   0.553  1.00  0.00           O  
HETATM   58  OXT VWW     1       3.456   6.137   1.141  1.00  0.00           O  
HETATM   59 OXT1 VWW     1      -5.616   2.347  -0.896  1.00  0.00           O  
HETATM   60  SG  VWW     1       1.451  -2.208  -0.581  1.00  0.00           S  
CONECT    1    4   54   54   58
CONECT    2    5   53   55   55
CONECT    3    6   56   56   59
CONECT    4    7   28   51
CONECT    5    8   29   52
CONECT    6    9   30   53
CONECT    7   16   31   33
CONECT    8   32   34   60
CONECT    9   17   17   18
CONECT   10   16   52   57   57
CONECT   11   14   36   38   60
CONECT   12   15   15   17   35
CONECT   13   15   18   18   37
CONECT   14   22   22   23
CONECT   15   39
CONECT   16   41   43
CONECT   17   40
CONECT   18   42
CONECT   19   20   20   21   44
CONECT   20   22   45
CONECT   21   23   23   46
CONECT   22   49
CONECT   23   50
CONECT   24   51
CONECT   25   52
CONECT   26   51
CONECT   27   53
CONECT   47   58
CONECT   48   59
END
"""
# pylint: enable=trailing-whitespace
# fmt: on


def _create_placeholder_mol():
  # Create a molecule with unordered names and explicit hydrogens.
  mol = rd_chem.MolFromSmiles('CN1CCC[C@H]1c2cccnc2')
  mol = rd_chem.AddHs(mol)
  atom_names = [string.ascii_uppercase[-i] for i in range(mol.GetNumAtoms())]

  for atom, name in zip(mol.GetAtoms(), atom_names, strict=True):
    atom.SetProp('atom_name', name)

  return mol, atom_names


class RdkitUtilTest(parameterized.TestCase):

  def test_reorder_atoms_from_smiles(self):
    mol, atom_names = _create_placeholder_mol()

    # Reorder.
    mol_reordered = rdkit_utils.sort_atoms_by_name(mol)

    # Check the names + sorting.
    new_names = [atom.GetProp('atom_name') for atom in mol_reordered.GetAtoms()]
    self.assertSequenceEqual(new_names, sorted(atom_names))

    self.assert_molecules_consistent(mol, mol_reordered)

  @parameterized.parameters([('KSC',), ('GND',), ('GP6',)])
  def test_reorder_atoms_from_ccd(self, ligand_code):
    ligand_mmcif = chemical_components.Ccd().get(ligand_code)
    mol = rdkit_utils.mol_from_ccd_cif(ligand_mmcif, force_parse=True)
    atom_names = [atom.GetProp('atom_name') for atom in mol.GetAtoms()]

    # Reorder.
    mol_reordered = rdkit_utils.sort_atoms_by_name(mol)

    # Check the names + sorting.
    new_names = [atom.GetProp('atom_name') for atom in mol_reordered.GetAtoms()]
    self.assertSequenceEqual(new_names, sorted(atom_names))
    self.assert_molecules_consistent(mol, mol_reordered)

  def assert_molecules_consistent(self, mol1, mol2):
    """Checks the atoms in two molecules are the same, discounting order."""
    atom_names1 = {
        atom.GetProp('atom_name'): i for i, atom in enumerate(mol1.GetAtoms())
    }
    atom_names2 = {
        atom.GetProp('atom_name'): i for i, atom in enumerate(mol2.GetAtoms())
    }

    # Check that the names match to the same atomic numbers.
    for name, index1 in atom_names1.items():
      atom1 = mol1.GetAtomWithIdx(index1)
      atom2 = mol2.GetAtomWithIdx(atom_names2[name])

      self.assertEqual(atom1.GetAtomicNum(), atom2.GetAtomicNum())
      self.assertEqual(atom1.GetFormalCharge(), atom2.GetFormalCharge())
      self.assertLen(atom1.GetBonds(), len(atom2.GetBonds()))

    # Name to index.
    atom_names1_reversed = {v: k for k, v in atom_names1.items()}
    atom_names2_reversed = {v: k for k, v in atom_names2.items()}

    # Check the bonds start and end on atoms with the same names.
    for bond1, bond2 in zip(mol1.GetBonds(), mol2.GetBonds(), strict=True):
      begin_name1 = atom_names1_reversed[bond1.GetBeginAtom().GetIdx()]
      begin_name2 = atom_names2_reversed[bond2.GetBeginAtom().GetIdx()]

      end_name1 = atom_names1_reversed[bond1.GetBeginAtom().GetIdx()]
      end_name2 = atom_names2_reversed[bond2.GetBeginAtom().GetIdx()]

      self.assertEqual(begin_name1, begin_name2)
      self.assertEqual(end_name1, end_name2)

  @parameterized.parameters(
      ('KSC', 27, 27), ('GND', 13, 12), ('GP6', 21, 20), ('NI', 1, 1)
  )
  def test_mol_from_ccd_cif(self, ligand_code, num_atoms, num_atoms_no_h):
    ligand_mmcif = chemical_components.Ccd().get(ligand_code)
    mol = rdkit_utils.mol_from_ccd_cif(
        ligand_mmcif,
        force_parse=False,
        sort_alphabetically=True,
        remove_hydrogens=True,
    )

    self.assertEqual(mol.GetNumAtoms(), num_atoms)
    # Shouldn't remove definitional hydrogens.
    self.assertEqual(rd_chem.RemoveHs(mol).GetNumAtoms(), num_atoms)

    # But should remove the others.
    self.assertEqual(rd_chem.RemoveAllHs(mol).GetNumAtoms(), num_atoms_no_h)

    # Should have atom names, they should be sorted by atom name.
    atom_names = [atom.GetProp('atom_name') for atom in mol.GetAtoms()]
    self.assertSequenceEqual(atom_names, sorted(atom_names))

  @parameterized.parameters(('KSC', 49), ('GND', 27), ('GP6', 33), ('NI', 1))
  def test_mol_from_ccd_cif_no_h_removal(self, ligand_code, num_atoms):
    ligand_mmcif = chemical_components.Ccd().get(ligand_code)
    mol = rdkit_utils.mol_from_ccd_cif(
        ligand_mmcif,
        force_parse=False,
        sort_alphabetically=False,
        remove_hydrogens=False,
    )

    self.assertEqual(mol.GetNumAtoms(), num_atoms)

  @parameterized.parameters('AOH', 'B1Z', 'B8B', 'BAZ', 'BBQ', 'A1H23')
  def test_mol_from_ccd_cif_fails(self, ligand_code):
    ligand_mmcif = chemical_components.Ccd().get(ligand_code)
    with self.assertRaises(
        rdkit_utils.MolFromMmcifError,
        msg=f'Expecting `MolFromMmcifError` for ligand {ligand_code}.',
    ):
      mol = rdkit_utils.mol_from_ccd_cif(ligand_mmcif, force_parse=True)
      if mol is None:
        self.fail(
            f'Ligand {ligand_code} was converted to RDKit Mol as None, '
            'but conversion is expected to fail with an exception.'
        )
      self.fail(
          f'Ligand {ligand_code} is not expected to be converted to RDKit '
          'Mol. Please check if that is an expected behavior.'
      )

  def assert_mmcif_consistent(
      self,
      ccd_cif,
      mol_mmcif,
      table_name,
      fields,
      allow_subset=False,
      include_only=lambda _: True,
  ):
    fields = set(fields)
    if not isinstance(ccd_cif, cif_dict.CifDict):
      ccd_cif = cif_dict.CifDict(ccd_cif)

    def get_records_to_compare(cif):
      table = cif.extract_loop_as_list(prefix=f'{table_name}.')
      records = []
      for record in table:
        if include_only(record):
          records.append(
              tuple((f, record[f'{table_name}.{f}']) for f in fields)
          )
      return records

    ccd_records = get_records_to_compare(ccd_cif)
    mol_records = get_records_to_compare(mol_mmcif)
    if allow_subset:
      ccd_records = set(ccd_records)
      mol_records = set(mol_records)
      self.assertCountEqual(ccd_records, ccd_records | mol_records)
    else:
      self.assertCountEqual(ccd_records, mol_records)

  @parameterized.product(
      ligand_code=('KSC', 'GND', 'GP6', 'NI'), include_hydrogens=(True, False)
  )
  def test_mol_to_ccd_cif(self, ligand_code, include_hydrogens):
    self.maxDiff = None  # pylint:disable=invalid-name

    ccd_cif = chemical_components.Ccd().get(ligand_code)
    mol = rdkit_utils.mol_from_ccd_cif(
        ccd_cif, force_parse=False, remove_hydrogens=not include_hydrogens
    )
    mol_mmcif = rdkit_utils.mol_to_ccd_cif(
        mol, ligand_code, include_hydrogens=include_hydrogens
    )
    self.assert_mmcif_consistent(
        ccd_cif,
        mol_mmcif,
        table_name='_chem_comp_atom',
        fields=(
            'atom_id',
            'charge',
            'type_symbol',
            'pdbx_model_Cartn_x_ideal',
            'pdbx_model_Cartn_y_ideal',
            'pdbx_model_Cartn_z_ideal',
        ),
        allow_subset=not include_hydrogens,
    )
    self.assert_mmcif_consistent(
        ccd_cif,
        mol_mmcif,
        table_name='_chem_comp_bond',
        fields=('atom_id_1', 'atom_id_2', 'pdbx_aromatic_flag'),
        allow_subset=not include_hydrogens,
    )
    self.assert_mmcif_consistent(
        ccd_cif,
        mol_mmcif,
        table_name='_chem_comp_bond',
        fields=('atom_id_1', 'atom_id_2', 'value_order'),
        allow_subset=not include_hydrogens,
        # Since Kekulization might differ for CCD representation and RDKit
        # molecule, we don't check bond types within aromatic rings.
        include_only=lambda x: x['_chem_comp_bond.pdbx_aromatic_flag'] == 'N',
    )

  @parameterized.parameters((False, 12), (True, 26))
  def test_mol_to_ccd_cif_generate_atom_data(
      self, include_hydrogens, expected_num_atoms
  ):
    mol = rd_chem.MolFromSmiles('CN1CCC[C@H]1c2cccnc2')
    AllChem.EmbedMultipleConfs(mol, 1, AllChem.ETKDGv3())

    cif = rdkit_utils.mol_to_ccd_cif(
        mol, 'XXX', include_hydrogens=include_hydrogens
    )
    atoms = cif.extract_loop_as_dict(
        '_chem_comp_atom.', '_chem_comp_atom.atom_id'
    )
    self.assertLen(atoms, expected_num_atoms)
    for atom_name, atom_fields in atoms.items():
      self.assertContainsSubset(
          (
              '_chem_comp_atom.pdbx_model_Cartn_x_ideal',
              '_chem_comp_atom.pdbx_model_Cartn_y_ideal',
              '_chem_comp_atom.pdbx_model_Cartn_z_ideal',
          ),
          atom_fields.keys(),
          msg=(
              'pdbx_model_Cartn_{x,y,z}_ideal expected to be populated, '
              f'but missing for {atom_name}.'
          ),
      )

  def test_mol_to_ccd_cif_embeds_smiles(self):
    mol = rd_chem.MolFromSmiles('CN1CCC[C@H]1c2cccnc2')
    AllChem.EmbedMultipleConfs(mol, 1, AllChem.ETKDGv3())

    cif = rdkit_utils.mol_to_ccd_cif(
        mol,
        'XXX',
        # A different SMILES descriptor for the same molecule.
        pdbx_smiles='n1cc(ccc1)C2N(C)CCC2',
    )
    self.assertEqual(cif['_chem_comp.pdbx_smiles'], ['n1cc(ccc1)C2N(C)CCC2'])

  def test_mol_to_ccd_cif_skip_ideal_conformer(self):
    mol = rd_chem.MolFromSmiles('CN1CCC[C@H]1c2cccnc2')

    cif = rdkit_utils.mol_to_ccd_cif(mol, 'XXX')
    atoms = cif.extract_loop_as_dict(
        '_chem_comp_atom.', '_chem_comp_atom.atom_id'
    )
    self.assertNotEmpty(atoms)
    for atom_fields in atoms.values():
      self.assertNotIn(
          '_chem_comp_atom.pdbx_model_Cartn_x_ideal', atom_fields.keys()
      )
      self.assertNotIn(
          '_chem_comp_atom.pdbx_model_Cartn_x_ideal', atom_fields.keys()
      )
      self.assertNotIn(
          '_chem_comp_atom.pdbx_model_Cartn_y_ideal', atom_fields.keys()
      )

  @parameterized.parameters([True, False])
  def test_assign_atom_names(self, keep_existing_names):
    ccd_cif = chemical_components.Ccd().get('VWW')
    # This method assigns the CCD atom names.
    mol = rdkit_utils.mol_from_ccd_cif(ccd_cif)

    # Remove a few of the atom names.
    for atom_idx in {
        8,  # CB2.
        19,  # CT1.
        29,  # OE1
    }:
      mol.GetAtomWithIdx(atom_idx).ClearProp('atom_name')

    mol = rdkit_utils.assign_atom_names_from_graph(mol, keep_existing_names)
    atom_names = [atom.GetProp('atom_name') for atom in mol.GetAtoms()]
    atomic_nums = [atom.GetAtomicNum() for atom in mol.GetAtoms()]

    # pyformat: disable
    self.assertListEqual(atomic_nums, [
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7,
        7, 7, 8, 8, 8, 8, 8, 8, 16
    ])
    # pyformat: enable

    self.assertListEqual(
        [atom.GetIdx() for atom in mol.GetAtoms()],
        list(range(mol.GetNumAtoms())),
    )

    if keep_existing_names:
      # CCD names, but:
      # CB2 -> C3 (C1, C2 are in the CCD names).
      # CT1 -> C4 (C1, C2 and C3 already specified).
      # OE1 -> O3 (O1 and O2 already specified).
      # pyformat: disable
      expected_atom_names = [
          'C', 'C1', 'C2', 'CA', 'CA1', 'CA2', 'CB', 'CB1', 'C3', 'CD', 'CD1',
          'CD11', 'CD2', 'CE', 'CE1', 'CG', 'CG1', 'CG2', 'CH', 'C4', 'CT2',
          'CZ1', 'CZ2', 'N', 'N1', 'N2', 'O', 'O1', 'O2', 'O3', 'OXT', 'OXT1',
          'SG'
      ]
      # pyformat: enable
    else:
      expected_atom_names = (
          *(f'C{i}' for i in range(1, 24)),
          *(f'N{i}' for i in range(1, 4)),
          *(f'O{i}' for i in range(1, 7)),
          'S1',
      )

    self.assertListEqual(atom_names, list(expected_atom_names))

  @parameterized.parameters([True, False])
  def test_assign_atom_names_two_letter_element(self, keep_existing_names):
    ccd_cif = chemical_components.Ccd().get('LH0')
    # This method assigns the CCD atom names.
    mol = rdkit_utils.mol_from_ccd_cif(ccd_cif)

    # Remove two of the atom names.
    for atom_idx in {
        7,  # C16
        17,  # CL
    }:
      mol.GetAtomWithIdx(atom_idx).ClearProp('atom_name')

    mol = rdkit_utils.assign_atom_names_from_graph(mol, keep_existing_names)
    atom_names = [atom.GetProp('atom_name') for atom in mol.GetAtoms()]
    atomic_nums = [atom.GetAtomicNum() for atom in mol.GetAtoms()]

    # pyformat: disable
    self.assertListEqual(atomic_nums, [
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 17, 17, 9, 9, 7, 7,
        7, 7, 8, 16
    ])
    # pyformat: enable

    self.assertListEqual(
        [atom.GetIdx() for atom in mol.GetAtoms()],
        list(range(mol.GetNumAtoms())),
    )

    if keep_existing_names:
      # CCD names, but:
      # C17 -> C16.
      # CL -> CL2.
      # pyformat: disable
      expected_atom_names = [
          'C', 'C1', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C2',
          'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'CL2', 'CL1', 'F', 'F1',
          'N', 'N1', 'N2', 'N3', 'O', 'S'
      ]
      # pyformat: enable
    else:
      expected_atom_names = (
          *(f'C{i}' for i in range(1, 18)),
          *(f'CL{i}' for i in range(1, 3)),
          *(f'F{i}' for i in range(1, 3)),
          *(f'N{i}' for i in range(1, 5)),
          'O1',
          'S1',
      )

    self.assertListEqual(atom_names, list(expected_atom_names))

  def test_mol_from_ccd_mmcif(self):
    ccd_mmcif = chemical_components.Ccd().get('VWW')
    mol = rdkit_utils.mol_from_ccd_cif(ccd_mmcif, remove_hydrogens=False)
    self.assertEqual(rd_chem.MolToPDBBlock(mol), _CCD_PDB)


if __name__ == '__main__':
  absltest.main()
