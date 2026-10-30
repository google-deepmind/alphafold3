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
from alphafold3.common import folding_input
from alphafold3.common import resources
from alphafold3.common.testing import data
from alphafold3.constants import chemical_components
from alphafold3.data import featurisation
from alphafold3.model.atom_layout import atom_layout  # pylint: disable=unused-import
from alphafold3.structure import test_utils
import numpy as np

_MINIMAL_USER_CCD = """data_MY-X7F
#
_chem_comp.id MY-X7F
_chem_comp.name '5,8-bis(oxidanyl)naphthalene-1,4-dione'
_chem_comp.type non-polymer
_chem_comp.formula 'C10 H6 O4'
_chem_comp.mon_nstd_parent_comp_id ?
_chem_comp.pdbx_synonyms ?
_chem_comp.formula_weight 190.152
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
MY-X7F C02 C 0 -1.418 -1.260 0.018
MY-X7F C03 C 0 -0.665 -2.503 -0.247
MY-X7F C04 C 0 0.677 -2.501 -0.235
MY-X7F C05 C 0 1.421 -1.257 0.043
MY-X7F C06 C 0 0.706 0.032 0.008
MY-X7F C07 C 0 -0.706 0.030 -0.004
MY-X7F C08 C 0 -1.397 1.240 -0.037
MY-X7F C10 C 0 -0.685 2.443 -0.057
MY-X7F C11 C 0 0.679 2.445 -0.045
MY-X7F C12 C 0 1.394 1.243 -0.013
MY-X7F O01 O 0 -2.611 -1.301 0.247
MY-X7F O09 O 0 -2.752 1.249 -0.049
MY-X7F O13 O 0 2.750 1.257 -0.001
MY-X7F O14 O 0 2.609 -1.294 0.298
MY-X7F H1 H 0 -1.199 -3.419 -0.452
MY-X7F H2 H 0 1.216 -3.416 -0.429
MY-X7F H3 H 0 -1.221 3.381 -0.082
MY-X7F H4 H 0 1.212 3.384 -0.062
MY-X7F H5 H 0 -3.154 1.271 0.830
MY-X7F H6 H 0 3.151 1.241 -0.880
#
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
O01 C02 DOUB N
O09 C08 SING N
C02 C03 SING N
C02 C07 SING N
C03 C04 DOUB N
C08 C07 DOUB Y
C08 C10 SING Y
C07 C06 SING Y
C10 C11 DOUB Y
C04 C05 SING N
C06 C05 SING N
C06 C12 DOUB Y
C11 C12 SING Y
C05 O14 DOUB N
C12 O13 SING N
C03 H1 SING N
C04 H2 SING N
C10 H3 SING N
C11 H4 SING N
O09 H5 SING N
O13 H6 SING N
#
"""


class FeaturisationTest(test_utils.StructureTestCase):

  def setUp(self):
    super().setUp()
    test_json = data.Data(resources.ROOT / 'common/test_data/').load(
        'alphafold_input.json'
    )
    self._folding_input = folding_input.Input.from_json(test_json)
    self._folding_input = self._folding_input.fill_missing_fields()

  def test_validate_input_msas(self):
    fold_input = folding_input.Input(
        name='test',
        chains=[
            folding_input.ProteinChain(
                id='A',
                sequence='YYDPETGTWY',
                ptms=[],
                unpaired_msa=None,
                paired_msa='',
                templates=[],
            )
        ],
        rng_seeds=[10],
    )
    with self.assertRaisesWithLiteralMatch(
        ValueError,
        'Protein chain 1 is missing unpaired MSA.',
    ):
      featurisation.validate_fold_input(fold_input)

  def test_validate_input_templates(self):
    fold_input = folding_input.Input(
        name='test',
        chains=[
            folding_input.ProteinChain(
                id='A',
                sequence='YYDPETGTWY',
                ptms=[],
                unpaired_msa='',
                paired_msa='',
                templates=None,
            )
        ],
        rng_seeds=[10],
    )
    with self.assertRaisesWithLiteralMatch(
        ValueError,
        'Protein chain 1 is missing Templates.',
    ):
      featurisation.validate_fold_input(fold_input)

  def test_featurise_input(self):
    batch = featurisation.featurise_input(
        self._folding_input, ccd=chemical_components.Ccd(), buckets=None
    )[0]
    self.assertContainsSubset(
        ('aatype', 'msa', 'template_aatype', 'tokens_to_queries:gather_idxs'),
        batch.keys(),
    )

  def test_featurise_user_ccd_ligand(self):
    af_input = folding_input.Input(
        name='test_folding_input',
        rng_seeds=[42],
        chains=[folding_input.Ligand(id='A', ccd_ids=['MY-X7F'])],
        user_ccd=_MINIMAL_USER_CCD,
    )
    ccd = chemical_components.Ccd(user_ccd=af_input.user_ccd)
    batch = featurisation.featurise_input(af_input, ccd=ccd, buckets=None)[0]
    # fmt: off
    self.assertEqual(
        batch['flat_output_layout'][()].atom_name.tolist(),  # pyrefly: ignore[missing-attribute]
        ['C02', 'C03', 'C04', 'C05', 'C06', 'C07', 'C08', 'C10', 'C11', 'C12',
         'O01', 'O09', 'O13', 'O14']
    )
    # fmt: on
    np.testing.assert_array_almost_equal(
        batch['ref_pos'][:, 0, :],
        [
            [-1.3, -0.4, 2.1],
            [-2.2, 0.5, 2.7],
            [-3.2, 1.1, 2.1],
            [-3.4, 0.8, 0.7],
            [-2.4, -0.0, 0.0],
            [-1.4, -0.7, 0.7],
            [-0.5, -1.6, 0.1],
            [-0.7, -1.8, -1.1],
            [-1.7, -1.2, -1.8],
            [-2.6, -0.3, -1.2],
            [-0.4, -0.9, 2.8],
            [0.4, -2.2, 0.8],
            [-3.6, 0.2, -2.0],
            [-4.3, 1.4, 0.0],
        ],
        # We are not interested in the exact values of the atom positions,
        # just that they are present and roughly the same.
        decimal=0,
    )

  def test_inhomogeneous_msa(self):
    fold_input = folding_input.Input(
        name='test',
        chains=[
            folding_input.ProteinChain(
                id='A',
                sequence='YYDPETGTWY',
                ptms=[],
                unpaired_msa='>query\nYYDPETGTWY\n>pattern2\nGGGGGGGGGG',
                paired_msa='',
                templates=[],
            ),
            folding_input.ProteinChain(
                id='B',
                sequence='YYDPETGTWY',
                ptms=[],
                unpaired_msa='>query\nYYDPETGTWY\n>pattern2\nPPPPPPPPPP',
                paired_msa='',
                templates=[],
            ),
        ],
        rng_seeds=[0],
    )
    features = featurisation.featurise_input(
        fold_input, ccd=chemical_components.Ccd(), buckets=None
    )[0]
    self.assertEqual(features['msa'].shape[1], 20)
    np.testing.assert_equal(features['msa'][2, :10], [7] * 10)
    np.testing.assert_equal(features['msa'][2, 10:], [14] * 10)

  def test_rna_msa_encoded_correctly(self):
    fold_input = folding_input.Input(
        name='test',
        chains=[
            folding_input.RnaChain(
                id='A',
                sequence='GGGGGG',
                unpaired_msa='>query\nGGGGGG\n>alignment1\nGAGGGG\n',
                modifications=[],
            ),
        ],
        rng_seeds=[0],
    )
    features = featurisation.featurise_input(
        fold_input, ccd=chemical_components.Ccd(), buckets=None
    )[0]
    self.assertEqual(features['num_alignments'], 3)
    alignments = features['msa'][: features['num_alignments']]
    # Check that the alignments are encoded as RNA residues, not amino acids.
    self.assertGreater(np.min(alignments), 20)
    self.assertLess(np.min(alignments), 26)

  def test_rdkit_conformer_generation_fallback(self):
    fold_input = folding_input.Input(
        name='test',
        chains=[folding_input.Ligand(id='B', ccd_ids=['BCR'])],
        rng_seeds=[1],
    )
    ccd = chemical_components.Ccd()

    with self.subTest('fallback'):
      features_no_fallback = featurisation.featurise_input(
          fold_input,
          ccd=ccd,
          buckets=None,
          ref_max_modified_date=datetime.date(1000, 1, 1),
      )[0]
      ref_pos_no_fallback = features_no_fallback['ref_pos']

      features_with_fallback = featurisation.featurise_input(
          fold_input,
          ccd=ccd,
          buckets=None,
          ref_max_modified_date=datetime.date(3000, 1, 1),
      )[0]
      ref_pos_with_fallback = features_with_fallback['ref_pos']

      self.assertEqual(ref_pos_no_fallback.shape, ref_pos_with_fallback.shape)
      self.assertNotEqual(
          ref_pos_no_fallback[0][0].tolist(),
          ref_pos_with_fallback[0][0].tolist(),
      )

    with self.subTest('missing_ref_max_modified_date'):
      with self.assertRaises(TypeError):
        featurisation.featurise_input(fold_input, ccd=ccd, buckets=None)

  @parameterized.named_parameters(
      dict(
          testcase_name='linked_glycan',
          ligand=folding_input.Ligand(id='B', ccd_ids=['NAG']),
          bonded_atom_pairs=[(('A', 1, 'ND2'), ('B', 1, 'C1'))],
      ),
      dict(
          testcase_name='standalone_branched_oligosaccharide',
          ligand=folding_input.Ligand(id='B', ccd_ids=['NAG', 'FUC']),
          bonded_atom_pairs=[(('B', 1, 'O6'), ('B', 2, 'C1'))],
      ),
      dict(
          testcase_name='linked_branched_oligosaccharide',
          ligand=folding_input.Ligand(id='B', ccd_ids=['NAG', 'FUC']),
          bonded_atom_pairs=[
              (('A', 1, 'ND2'), ('B', 1, 'C1')),  # Polymer-ligand bond.
              (('B', 1, 'O6'), ('B', 2, 'C1')),  # Ligand-ligand bond.
          ],
      ),
  )
  def test_featurisation_with_fix_standalone_glycans_unchanged(
      self, ligand, bonded_atom_pairs
  ):
    fold_input = folding_input.Input(
        name='9SI4',
        chains=[
            folding_input.ProteinChain(
                id='A',
                sequence='NAAAN',
                ptms=[],
                unpaired_msa='',
                paired_msa='',
                templates=[],
            ),
            ligand,
        ],
        bonded_atom_pairs=bonded_atom_pairs,
        rng_seeds=[1],
    )
    ccd = chemical_components.Ccd()
    features_unfixed = featurisation.featurise_input(
        fold_input, ccd=ccd, buckets=None, fix_standalone_glycans=False
    )[0]
    features_fixed = featurisation.featurise_input(
        fold_input, ccd=ccd, buckets=None, fix_standalone_glycans=True
    )[0]
    self.assertCountEqual(features_unfixed.keys(), features_fixed.keys())

    for k in features_fixed.keys():
      fixed = features_fixed[k]
      unfixed = features_unfixed[k]
      if k in ['cleaned_struc', 'empty_output_struc']:
        self.assertStructuresEqual(fixed.item(), unfixed.item())
      elif k in ['flat_output_layout', 'token_atoms_layout']:
        self.assertEqual(fixed.item(), unfixed.item(), msg=k)
      else:
        np.testing.assert_equal(fixed, unfixed, strict=True, err_msg=k)

  def test_featurisation_with_fix_standalone_glycans_changed(self):
    fold_input = folding_input.Input(
        name='9SI4',
        chains=[
            folding_input.ProteinChain(
                id='A',
                sequence='NAAAN',
                ptms=[],
                unpaired_msa='',
                paired_msa='',
                templates=[],
            ),
            folding_input.Ligand(id='B', ccd_ids=['NAG']),
        ],
        bonded_atom_pairs=[],
        rng_seeds=[1],
    )
    ccd = chemical_components.Ccd()
    features_unfixed = featurisation.featurise_input(
        fold_input, ccd=ccd, buckets=None, fix_standalone_glycans=False
    )[0]
    features_fixed = featurisation.featurise_input(
        fold_input, ccd=ccd, buckets=None, fix_standalone_glycans=True
    )[0]
    self.assertCountEqual(features_unfixed.keys(), features_fixed.keys())

    for k in features_fixed.keys():
      fixed = features_fixed[k]
      unfixed = features_unfixed[k]
      if k in ['empty_output_struc', 'cleaned_struc']:
        self.assertEqual(fixed.item().num_atoms, unfixed.item().num_atoms + 1)
      elif k in ['token_atoms_layout']:
        self.assertEqual(fixed.item().atom_name.shape, (5 + 15, 24), msg=k)
        self.assertEqual(unfixed.item().atom_name.shape, (5 + 14, 24), msg=k)
      elif k in ['flat_output_layout']:
        self.assertEqual(fixed.item().atom_name.shape, (47,), msg=k)
        self.assertEqual(unfixed.item().atom_name.shape, (46,), msg=k)
      elif k in [
          'tokens_to_ligand_ligand_bonds:gather_idxs',
          'tokens_to_ligand_ligand_bonds:gather_mask',
      ]:
        # Check that the shape difference is exactly 10 tokens in total.
        shape_diff = np.array(fixed.shape) - np.array(unfixed.shape)
        self.assertEqual(shape_diff.sum(), 10, msg=f'{k}: {shape_diff=}')
      elif fixed.shape != unfixed.shape:
        # Check that the shape difference is exactly 1 atom in total.
        shape_diff = np.array(fixed.shape) - np.array(unfixed.shape)
        self.assertEqual(shape_diff.sum(), 1, msg=f'{k}: {shape_diff=}')
      elif not np.equal(fixed, unfixed).all():
        if fixed.size == 1:
          self.assertEqual(fixed.item(), unfixed.item() + 1, msg=k)
        elif k == 'tokens_to_keys:gather_idxs':
          self.assertEqual(np.sum(fixed != unfixed), 249, msg=k)
        elif k == 'tokens_to_keys:gather_mask':
          self.assertEqual(np.sum(fixed != unfixed), 15, msg=k)
        elif k == 'queries_to_keys:gather_idxs':
          self.assertEqual(np.sum(fixed != unfixed), 613, msg=k)
        else:
          # Only one atom is different between fixed and unfixed.
          self.assertEqual(np.sum(fixed != unfixed), 1, msg=k)
      else:
        self.assertIn(
            k,
            [
                'num_alignments',
                'queries_to_keys:gather_mask',
                'queries_to_keys:input_shape',
                'queries_to_token_atoms:input_shape',
                'polymer_ligand_bonds',
                'ligand_ligand_bonds',
            ],
        )
        np.testing.assert_equal(fixed, unfixed, strict=True, err_msg=k)


if __name__ == '__main__':
  absltest.main()
