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
from alphafold3.model import confidence_types
from alphafold3.model import model
from alphafold3.model import post_processing
from alphafold3.structure import mmcif
from etils import epath
import numpy as np
import zstandard

_DATA_ROOT = resources.ROOT / 'test_data/miniature_databases/pdb_mmcif'


def _create_structure(protein_name: str, ligand: bool) -> structure.Structure:
  cif_path = os.path.join(_DATA_ROOT, f'{protein_name}.cif')
  cif_str = resources.open_resource(cif_path, mode='rt').read()
  struc = structure.from_mmcif(cif_str, fix_arginines=True, include_bonds=False)
  return struc.filter_to_entity_type(protein=True, ligand=ligand)


class PostProcessingTest(parameterized.TestCase):

  def setUp(self):
    super().setUp()

    self.inference_result = model.InferenceResult(
        predicted_structure=_create_structure(
            protein_name='6d0i', ligand=False
        ),
        numerical_data={
            'full_pde': np.zeros((6, 6), dtype=np.float32),
            'full_pae': np.zeros((6, 6), dtype=np.float32),
            'contact_probs': np.zeros((6, 6), dtype=np.float32),
        },
        metadata={
            'predicted_distance_error': 0.0,
            'predicted_dist': np.zeros((6, 6), dtype=np.float32),
            'ranking_score': 0.0,
            'fraction_disordered': 0.0,
            'has_clash': False,
            'chain_pair_pde_mean': np.zeros((4, 4), dtype=np.float64),
            'chain_pair_pde_min': np.zeros((4, 4), dtype=np.float64),
            'intra_chain_single_pde': np.array([0.0, 0.0, 0.0, 0.0]),
            'cross_chain_single_pde': np.array([0.0, 0.0, 0.0, 0.0]),
            'chain_pair_pae_mean': np.zeros((4, 4), dtype=np.float64),
            'chain_pair_pae_min': np.zeros((4, 4), dtype=np.float64),
            'chain_pair_iptm': np.zeros((4, 4), dtype=np.float64),
            'ptm': 0.0,
            'iptm': 0.0,
            'iptm_ichain': np.array([0.0, 0.0, 0.0, 0.0]),
            'iptm_xchain': np.array([0.0, 0.0, 0.0, 0.0]),
            'token_chain_ids': ['A', 'A', 'B', 'B', 'C', 'D'],  # pyrefly: ignore[bad-assignment]
            'token_res_ids': [1, 2, 1, 2, 1, 1],  # pyrefly: ignore[bad-assignment]
        },
        debug_outputs={},
    )
    self._output_dir = epath.Path(self.create_tempdir())

  @parameterized.named_parameters(
      dict(
          testcase_name='compress',
          compress=True,
          cif_path='model.cif.zst',
          confidences_path='confidences.json.zst',
          opener=zstandard.open,
          keep_license=True,
      ),
      dict(
          testcase_name='no_compress',
          compress=False,
          cif_path='model.cif',
          confidences_path='confidences.json',
          opener=open,
          keep_license=True,
      ),
      dict(
          testcase_name='no_keep_license',
          compress=False,
          cif_path='model.cif',
          confidences_path='confidences.json',
          opener=open,
          keep_license=False,
      ),
  )
  def test_write_output(
      self,
      compress: bool,
      cif_path: str,
      confidences_path: str,
      opener,
      keep_license: bool,
  ):
    post_processing.write_output(
        self.inference_result,
        output_dir=self._output_dir,
        terms_of_use='ALPHAFOLD 3 OUTPUT TERMS OF USE',
        compress=compress,
        keep_license=keep_license,
    )

    with self.subTest('expected_files'):
      self.assertSameElements(
          os.listdir(self._output_dir),
          [
              cif_path,
              'summary_confidences.json',
              confidences_path,
              'TERMS_OF_USE.md',
          ],
      )

    with self.subTest('cif_file'):
      with opener(self._output_dir / cif_path, 'rt') as f:
        pred_cif = f.read()
      pred_cif_dict = mmcif.from_string(pred_cif)
      self.assertEqual(pred_cif_dict['_entry.id'], ['6D0I'])
      license_text = 'Non-commercial use only'
      comment_text = 'By using this file you agree to the legally binding terms'
      if keep_license:
        self.assertIn(license_text, pred_cif)
        self.assertIn(comment_text, pred_cif)
      else:
        self.assertNotIn(license_text, pred_cif)
        self.assertNotIn(comment_text, pred_cif)

    with self.subTest('summary_file'):
      with open(self._output_dir / 'summary_confidences.json', 'rt') as f:
        pred_summary = f.read()
      structure_confidence_summary = (
          confidence_types.StructureConfidenceSummary.from_json(pred_summary)
      )
      self.assertEqual(structure_confidence_summary.ptm, 0.0)
      self.assertEqual(structure_confidence_summary.iptm, 0.0)
      self.assertEqual(structure_confidence_summary.ranking_score, 0.0)

    with self.subTest('full_data_file'):
      with opener(self._output_dir / confidences_path, 'rt') as f:
        pred_full_data = f.read()
      full_data = confidence_types.StructureConfidenceFull.from_json(
          pred_full_data
      )
      self.assertEqual(set(full_data.atom_chain_ids), {'A', 'B', 'C', 'D'})

    with self.subTest('terms_of_use_file'):
      with open(self._output_dir / 'TERMS_OF_USE.md', 'rt') as f:
        terms_of_use = f.read()
      self.assertEqual(terms_of_use, 'ALPHAFOLD 3 OUTPUT TERMS OF USE')

  def test_write_embeddings(self):
    mock_embeddings = {'a': np.array([1, 2, 3]), 'b': np.array([4, 5, 6])}
    post_processing.write_embeddings(
        mock_embeddings,
        output_dir=self._output_dir,
        name='test',
    )

    npz_path = self._output_dir / 'test_embeddings.npz'
    self.assertTrue(os.path.exists(npz_path))
    with np.load(npz_path) as data:
      self.assertCountEqual(data.files, ['a', 'b'])
      np.testing.assert_array_equal(data['a'], mock_embeddings['a'])
      np.testing.assert_array_equal(data['b'], mock_embeddings['b'])

  @parameterized.named_parameters(
      ('keep_license', True), ('discard_license', True)
  )
  def test_post_process_inference_result(self, keep_license: bool):
    result = post_processing.post_process_inference_result(
        self.inference_result, keep_license=keep_license
    )

    with self.subTest('cif file'):
      pred_cif_dict = mmcif.from_string(result.cif)
      self.assertEqual(pred_cif_dict['_entry.id'], ['6D0I'])
      cif_str = result.cif.decode('utf-8')
      comment_str = 'By using this file you agree to the legally binding terms'
      if keep_license:
        self.assertIn(comment_str, cif_str)
      else:
        self.assertNotIn(comment_str, cif_str)
    with self.subTest('summary file'):
      structure_confidence_summary = (
          confidence_types.StructureConfidenceSummary.from_json(
              result.structure_confidence_summary_json.decode('utf-8')
          )
      )
      self.assertEqual(structure_confidence_summary.ptm, 0.0)
      self.assertEqual(structure_confidence_summary.iptm, 0.0)
      self.assertEqual(structure_confidence_summary.ranking_score, 0.0)
    with self.subTest('full data file'):
      full_data = confidence_types.StructureConfidenceFull.from_json(
          result.structure_full_data_json.decode('utf-8')
      )
      self.assertEqual(set(full_data.atom_chain_ids), {'A', 'B', 'C', 'D'})
    with self.subTest('mean conf1d'):
      self.assertAlmostEqual(result.mean_confidence_1d, 18.052, places=3)


if __name__ == '__main__':
  absltest.main()
