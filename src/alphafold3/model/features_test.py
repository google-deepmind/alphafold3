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

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3 import structure
from alphafold3.common import folding_input
from alphafold3.constants import chemical_components
from alphafold3.constants import mmcif_names
from alphafold3.constants import periodic_table
from alphafold3.model import features
from alphafold3.model.atom_layout import atom_layout
from alphafold3.model.pipeline import inter_chain_bonds
from alphafold3.structure import chemical_components as struc_chem_comps
import numpy as np


class FeaturesTest(parameterized.TestCase):

  @parameterized.parameters([
      dict(res_name='HEM'),
      dict(res_name='G'),
      dict(res_name='SRO'),
  ])
  def test_get_reference_ccd(self, res_name):
    ccd_cif = chemical_components.Ccd().get(res_name.split('-')[0])
    ccd_atom_names = ccd_cif['_chem_comp_atom.atom_id']  # pyrefly: ignore[unsupported-operation]
    ccd_atom_x = [
        float(x) for x in ccd_cif['_chem_comp_atom.pdbx_model_Cartn_x_ideal']  # pyrefly: ignore[unsupported-operation]
    ]
    ccd_type_symbols = ccd_cif['_chem_comp_atom.type_symbol']  # pyrefly: ignore[unsupported-operation]
    ccd_charges = [
        0 if x == '?' else int(x) for x in ccd_cif['_chem_comp_atom.charge']  # pyrefly: ignore[unsupported-operation]
    ]
    ccd_elements = [
        periodic_table.ATOMIC_NUMBER.get(elem_type.capitalize(), 0)
        for elem_type in ccd_type_symbols
    ]

    chem_comp = struc_chem_comps.ChemicalComponentsData({})
    ref, _, _ = features.get_reference(
        res_name=res_name,
        chemical_components_data=chem_comp,
        ccd=chemical_components.Ccd(),
        random_state=np.random.RandomState(),
        ref_max_modified_date=datetime.date.fromisoformat('3000-01-01'),
        conformer_max_iterations=None,
    )
    atom_x_ref = [ref[n]['positions'][0] for n in ccd_atom_names]
    charges_ref = [ref[n]['charge'] for n in ccd_atom_names]
    elements_ref = [ref[n]['element'] for n in ccd_atom_names]
    self.assertNotEqual(atom_x_ref, ccd_atom_x)  # Check augmentation applied.
    self.assertEqual(charges_ref, ccd_charges)
    self.assertEqual(elements_ref, ccd_elements)

  def test_get_reference_smiles(self):
    ligand_name = 'SRO-test'
    smiles = 'c1cc2c(cc1O)c(c[nH]2)CCN'
    chem_comp = struc_chem_comps.ChemicalComponentsData({
        ligand_name: struc_chem_comps.ChemCompEntry(
            type='non-polymer', pdbx_smiles=smiles
        )
    })
    ref, _, _ = features.get_reference(
        res_name=ligand_name,
        chemical_components_data=chem_comp,
        ccd=chemical_components.Ccd(),
        random_state=np.random.RandomState(),
        ref_max_modified_date=datetime.date.fromisoformat('3000-01-01'),
        conformer_max_iterations=None,
    )
    self.assertTrue(ref)

  def test_get_reference_smiles_with_missing_smiles(self):
    ligand_name = 'SRO-test'
    chem_comp = struc_chem_comps.ChemicalComponentsData({
        ligand_name: struc_chem_comps.ChemCompEntry(
            type='non-polymer', pdbx_smiles=None
        )
    })
    with self.assertRaises(ValueError):
      features.get_reference(
          res_name=ligand_name,
          chemical_components_data=chem_comp,
          ccd=chemical_components.Ccd(),
          random_state=np.random.RandomState(),
          ref_max_modified_date=datetime.date.fromisoformat('3000-01-01'),
          conformer_max_iterations=None,
      )

  def test_ligand_ligand_bonds_info(self):
    zeros = np.zeros(4, dtype=np.float32)
    bonds_struct = structure.from_atom_arrays(
        chain_id=np.array(['A', 'A', 'B', 'B'], dtype=object),
        chain_type=np.array(
            [
                mmcif_names.BRANCHED_CHAIN,
                mmcif_names.BRANCHED_CHAIN,
                mmcif_names.PROTEIN_CHAIN,
                mmcif_names.PROTEIN_CHAIN,
            ],
            dtype=object,
        ),
        res_id=np.array([1, 2, 1, 1], dtype=np.int32),
        res_name=np.array(['NAG', 'NAG', 'ALA', 'ALA'], dtype=object),
        atom_name=np.array(['C1', 'C1', 'ND2', 'CA'], dtype=object),
        atom_element=np.array(['C', 'C', 'N', 'C'], dtype=object),
        atom_x=zeros,
        atom_y=zeros,
        atom_z=zeros,
    ).add_bonds(
        bonded_atom_pairs=[
            (('A', 1, 'C1'), ('A', 2, 'C1')),  #  NAG1-C1 <-> NAG2-C1
            (('B', 1, 'ND2'), ('A', 2, 'C1')),  # ASN-ND2 <-> NAG2-C1
        ],
        bond_type=mmcif_names.COVALENT_BOND,
    )

    all_tokens = atom_layout.atom_layout_from_structure(bonds_struct)
    bond_layout = inter_chain_bonds.get_ligand_ligand_bonds(
        bonds_struct, only_glycan_ligands=False
    )
    padding_shapes = features.PaddingShapes(
        num_tokens=4,
        msa_size=0,
        num_chains=2,
        num_templates=0,
        num_atoms=4,
    )
    ligand_ligand_info = features.LigandLigandBondInfo.compute_features(
        all_tokens, bond_layout, padding_shapes
    )
    np.testing.assert_equal(
        np.array([[0, 1]] + [[0, 0]] * 39),
        ligand_ligand_info.tokens_to_ligand_ligand_bonds.gather_idxs,
    )
    np.testing.assert_equal(
        np.array([[True, True]] + [[False, False]] * 39),
        ligand_ligand_info.tokens_to_ligand_ligand_bonds.gather_mask,
    )
    np.testing.assert_equal(
        np.array([4]),
        ligand_ligand_info.tokens_to_ligand_ligand_bonds.input_shape,
    )

  def test_colinear_ligand_frames(self):
    all_tokens = atom_layout.AtomLayout(
        atom_name=np.array(['A', 'B', 'C', 'D'], dtype=object),
        atom_element=np.array(['A', 'B', 'C', 'D'], dtype=object),
        res_id=np.array([1, 1, 1, 1], dtype=int),
        res_name=np.array(['NAG', 'NAG', 'NAG', 'NAG'], dtype=object),
        chain_type=np.array([mmcif_names.NON_POLYMER_CHAIN] * 4, dtype=object),
        chain_id=np.array(['A', 'A', 'A', 'A'], dtype=object),
    )
    all_token_atoms_layout = atom_layout.AtomLayout(
        atom_name=all_tokens.atom_name[:, None],
        atom_element=all_tokens.atom_element[:, None],  # pyrefly: ignore[unsupported-operation]
        res_id=all_tokens.res_id[:, None],
        res_name=all_tokens.res_name[:, None],  # pyrefly: ignore[unsupported-operation]
        chain_type=all_tokens.chain_type[:, None],  # pyrefly: ignore[unsupported-operation]
        chain_id=all_tokens.chain_id[:, None],
    )
    ref_positions = np.stack(
        [
            np.array([0, 5, 1, 12]),
            np.zeros(4, dtype=np.float32),
            np.zeros(4, dtype=np.float32),
        ],
        axis=1,
    )
    padding_shapes = features.PaddingShapes(
        num_tokens=4,
        msa_size=0,
        num_chains=1,
        num_templates=0,
        num_atoms=4,
    )
    ref_structure = features.RefStructure(
        positions=ref_positions[:, None, :],
        mask=np.ones([4, 1], dtype=bool),
        element=np.zeros([4, 1]),
        charge=np.zeros([4, 1]),
        atom_name_chars=np.zeros([4, 1]),
        ref_space_uid=np.zeros([4, 1]),
    )

    frames = features.Frames.compute_features(
        all_tokens,
        all_token_atoms_layout,
        ref_structure,
        padding_shapes,
    )
    np.testing.assert_array_equal(frames.mask, False)

  def test_msa_features(self):
    num_res = 10
    msa_size = 5
    all_tokens = atom_layout.AtomLayout(
        atom_name=np.array(['CA'] * num_res, dtype=object),
        atom_element=np.array(['C'] * num_res, dtype=object),
        res_id=np.array([1, 2, 3, 4, 5, 1, 2, 3, 4, 5], dtype=int),
        res_name=np.array(['ALA'] * 5 + ['GLY'] * 5, dtype=object),
        chain_type=np.array(
            [mmcif_names.PROTEIN_CHAIN] * num_res, dtype=object
        ),
        chain_id=np.array(['A'] * 5 + ['B'] * 5, dtype=object),
    )

    padding_shapes = features.PaddingShapes(
        num_tokens=num_res,
        msa_size=msa_size,
        num_chains=2,
        num_templates=0,
        num_atoms=num_res,
    )

    msa_features = features.MSA.compute_features(
        all_tokens=all_tokens,
        standard_token_idxs=np.arange(10, dtype=np.int32),
        padding_shapes=padding_shapes,
        unpaired_msa_by_chain_id={
            'A': '>q\nAAAAA\n>h1\nAAAAC',
            'B': '>q\nGGGGG\n>g1\n-----\n>g2\n-----\n>h1\nGGGGA\n>h2\nGGGGC',
        },
        paired_msa_by_chain_id={
            'A': '>q\nAAAAA',
            'B': '>q\nGGGGG',
        },
        logging_name='test',
        max_paired_sequence_per_species=1,
    )

    with self.subTest('shape_correct'):
      self.assertEqual(msa_features.rows.shape, (msa_size, num_res))
      self.assertEqual(msa_features.mask.shape, (msa_size, num_res))
      self.assertEqual(msa_features.deletion_matrix.shape, (msa_size, num_res))
      self.assertEqual(msa_features.profile.shape, (num_res, 21 + 5 + 5))
      self.assertEqual(msa_features.deletion_mean.shape, (num_res,))
      self.assertEqual(msa_features.num_alignments.shape, ())

    with self.subTest('mask_all_true'):
      self.assertTrue(np.all(msa_features.mask))

    with self.subTest('deletion_matrix_empty'):
      self.assertEqual(np.count_nonzero(msa_features.deletion_matrix), 0)

    with self.subTest('num_alignments_correct'):
      self.assertEqual(msa_features.num_alignments, msa_size)

    with self.subTest('msa_correct'):
      # fmt: off
      # pylint: disable=bad-whitespace
      exp_msa = np.array(
          [
              [ 0,  0,  0,  0,  0,  7,  7,  7,  7,  7],  # AAAAAGGGGG
              [ 0,  0,  0,  0,  4, 21, 21, 21, 21, 21],  # AAAAC-----
              [21, 21, 21, 21, 21, 21, 21, 21, 21, 21],  # ----------
              [21, 21, 21, 21, 21,  7,  7,  7,  7,  0],  # -----GGGGA
              [21, 21, 21, 21, 21,  7,  7,  7,  7,  4],  # -----GGGGD
          ],
          dtype=np.int8,
      )
      # pylint: enable=bad-whitespace
      # fmt: on
      np.testing.assert_array_equal(msa_features.rows, exp_msa, strict=True)

  def test_template_features(self):
    num_res = 5
    max_templates = 4

    all_tokens = atom_layout.AtomLayout(
        atom_name=np.array(['CA'] * num_res, dtype=object),
        atom_element=np.array(['C'] * num_res, dtype=object),
        res_id=np.array([1, 2, 3, 4, 5], dtype=int),
        res_name=np.array(['ALA', 'GLY', 'ALA', 'GLY', 'ALA'], dtype=object),
        chain_type=np.array(
            [mmcif_names.PROTEIN_CHAIN] * num_res, dtype=object
        ),
        chain_id=np.array(['A'] * num_res, dtype=object),
    )

    padding_shapes = features.PaddingShapes(
        num_tokens=num_res,
        msa_size=0,
        num_chains=1,
        num_templates=max_templates,
        num_atoms=num_res,
    )

    # Build a template structure from a real sequence, convert to mmCIF.
    template_struc = structure.from_sequences_and_bonds(
        sequences=['AGRGMMM'],
        chain_types=[mmcif_names.PROTEIN_CHAIN],
        sequence_formats=[structure.SequenceFormat.FASTA],
        bonded_atom_pairs=None,
        ccd=chemical_components.Ccd(),
    ).copy_and_update_globals(release_date=datetime.date(1989, 11, 17))

    templates_by_chain_id = {
        'A': [
            folding_input.Template(
                mmcif=template_struc.to_mmcif(),
                query_to_template_map={0: 0, 1: 1, 2: 2, 3: 3},
            ),
        ],
    }

    template_features = features.Templates.compute_features(
        all_tokens=all_tokens,
        standard_token_idxs=np.arange(num_res, dtype=np.int32),
        padding_shapes=padding_shapes,
        templates_by_chain_id=templates_by_chain_id,
        max_templates=max_templates,
        logging_name='test',
    )

    with self.subTest('shapes_correct'):
      self.assertEqual(template_features.aatype.shape, (max_templates, num_res))
      self.assertEqual(
          template_features.atom_positions.shape,
          (max_templates, num_res, 24, 3),
      )
      self.assertEqual(
          template_features.atom_mask.shape, (max_templates, num_res, 24)
      )

    with self.subTest('aatype_aligned_residues'):
      # First 4 res are aligned (ALA=0, GLY=7, ARG=1, GLY=7), last is gap (21).
      np.testing.assert_array_equal(
          template_features.aatype[0], [0, 7, 1, 7, 21]
      )

    with self.subTest('aatype_padding_templates_zero'):
      # Templates 1-3 are padding, should be all zeros.
      np.testing.assert_array_equal(template_features.aatype[1:], 0)

    with self.subTest('atom_mask_aligned_residues'):
      # Aligned residues have 1 atom present, unaligned are all False.
      self.assertEqual(np.sum(template_features.atom_mask[0, :4]), 4)
      np.testing.assert_array_equal(template_features.atom_mask[0, 4:], False)

    with self.subTest('atom_mask_padding_templates_false'):
      np.testing.assert_array_equal(template_features.atom_mask[1:], False)


if __name__ == '__main__':
  absltest.main()
