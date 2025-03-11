# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Tests end-to-end running of AlphaFold 3."""

import contextlib
import csv
import dataclasses
import datetime
import difflib
import json
import os
import pathlib
import pickle

from absl import logging
from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.common import folding_input
from alphafold3.common import resources
from alphafold3.common.testing import data as testing_data
from alphafold3.data import pipeline
from alphafold3.model.scoring import alignment
from alphafold3.structure import test_utils
import jax
import numpy as np

import run_alphafold
import shutil


_JACKHMMER_BINARY_PATH = shutil.which('jackhmmer')
_NHMMER_BINARY_PATH = shutil.which('nhmmer')
_HMMALIGN_BINARY_PATH = shutil.which('hmmalign')
_HMMSEARCH_BINARY_PATH = shutil.which('hmmsearch')
_HMMBUILD_BINARY_PATH = shutil.which('hmmbuild')


@contextlib.contextmanager
def _output(name: str):
  with open(result_path := f'{absltest.TEST_TMPDIR.value}/{name}', "wb") as f:
    yield result_path, f


jax.config.update('jax_enable_compilation_cache', False)


def _generate_diff(actual: str, expected: str) -> str:
  return '\n'.join(
      difflib.unified_diff(
          expected.split('\n'),
          actual.split('\n'),
          fromfile='expected',
          tofile='actual',
          lineterm='',
      )
  )


class InferenceTest(test_utils.StructureTestCase):
  """Test AlphaFold 3 inference."""

  def setUp(self):
    super().setUp()
    small_bfd_database_path = testing_data.Data(
        resources.ROOT
        / 'test_data/miniature_databases/bfd-first_non_consensus_sequences__subsampled_1000.fasta'
    ).path()
    mgnify_database_path = testing_data.Data(
        resources.ROOT
        / 'test_data/miniature_databases/mgy_clusters__subsampled_1000.fa'
    ).path()
    uniprot_cluster_annot_database_path = testing_data.Data(
        resources.ROOT
        / 'test_data/miniature_databases/uniprot_all__subsampled_1000.fasta'
    ).path()
    uniref90_database_path = testing_data.Data(
        resources.ROOT
        / 'test_data/miniature_databases/uniref90__subsampled_1000.fasta'
    ).path()
    ntrna_database_path = testing_data.Data(
        resources.ROOT
        / 'test_data/miniature_databases/nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq__subsampled_1000.fasta'
    ).path()
    rfam_database_path = testing_data.Data(
        resources.ROOT
        / 'test_data/miniature_databases/rfam_14_4_clustered_rep_seq__subsampled_1000.fasta'
    ).path()
    rna_central_database_path = testing_data.Data(
        resources.ROOT
        / 'test_data/miniature_databases/rnacentral_active_seq_id_90_cov_80_linclust__subsampled_1000.fasta'
    ).path()
    pdb_database_path = testing_data.Data(
        resources.ROOT / 'test_data/miniature_databases/pdb_mmcif'
    ).path()
    seqres_database_path = testing_data.Data(
        resources.ROOT
        / 'test_data/miniature_databases/pdb_seqres_2022_09_28__subsampled_1000.fasta'
    ).path()

    self._data_pipeline_config = pipeline.DataPipelineConfig(
        jackhmmer_binary_path=_JACKHMMER_BINARY_PATH,
        nhmmer_binary_path=_NHMMER_BINARY_PATH,
        hmmalign_binary_path=_HMMALIGN_BINARY_PATH,
        hmmsearch_binary_path=_HMMSEARCH_BINARY_PATH,
        hmmbuild_binary_path=_HMMBUILD_BINARY_PATH,
        small_bfd_database_path=small_bfd_database_path,
        mgnify_database_path=mgnify_database_path,
        uniprot_cluster_annot_database_path=uniprot_cluster_annot_database_path,
        uniref90_database_path=uniref90_database_path,
        ntrna_database_path=ntrna_database_path,
        rfam_database_path=rfam_database_path,
        rna_central_database_path=rna_central_database_path,
        pdb_database_path=pdb_database_path,
        seqres_database_path=seqres_database_path,
        max_template_date=datetime.date(2021, 9, 30),
    )
    test_input = {
        'name': '5tgy',
        'modelSeeds': [1234],
        'sequences': [
            {
                'protein': {
                    'id': 'P',
                    'sequence': 'SEFEKLRQTGDELVQAFQRLREIFDKGDDDSLEQVLEEIEELIQKHRQLFDNRQEAADTEAAKQGDQWVQLFQRFREAIDKGDKDSLEQLLEELEQALQKIRELAEKKN',
                    'modifications': [],
                    'unpairedMsa': None,
                    'pairedMsa': None,
                }
            },
            {'ligand': {'id': 'LL', 'ccdCodes': ['7BU']}},
        ],
        'dialect': folding_input.JSON_DIALECT,
        'version': folding_input.JSON_VERSION,
    }
    self._test_input_json = json.dumps(test_input)
    self._model_config = run_alphafold.make_model_config(
        return_embeddings=True, flash_attention_implementation='triton'
    )
    self._runner = run_alphafold.ModelRunner(
        config=self._model_config,
        device=jax.local_devices()[0],
        model_dir=pathlib.Path(run_alphafold.MODEL_DIR.value),
    )

  def test_model_inference(self):
    """Run model inference and assert that output exists."""
    featurised_examples = pickle.loads(
        (resources.ROOT / 'test_data' / 'featurised_example.pkl').read_bytes()
    )

    self.assertLen(featurised_examples, 1)
    featurised_example = featurised_examples[0]
    result = self._runner.run_inference(
        featurised_example, jax.random.PRNGKey(0)
    )
    self.assertIsNotNone(result)
    _, embeddings = self._runner.extract_inference_results_and_maybe_embeddings(
        batch=featurised_example, result=result, target_name='target'
    )
    self.assertLen(embeddings, 2)

  def test_process_fold_input_runs_only_inference(self):
    with self.assertRaisesRegex(ValueError, 'missing unpaired MSA.'):
      run_alphafold.process_fold_input(
          fold_input=folding_input.Input.from_json(self._test_input_json),
          # No data pipeline config, so featurisation will run first, and fail
          # since the input is missing MSAs.
          data_pipeline_config=None,
          model_runner=self._runner,
          output_dir=self.create_tempdir().full_path,
      )

  @parameterized.named_parameters(
      {
          'testcase_name': 'default_bucket',
          'bucket': None,
          'seed': 1,
      },
      {
          'testcase_name': 'bucket_1024',
          'bucket': 1024,
          'seed': 42,
      },
  )
  def test_inference(self, bucket, seed):
    """Run AlphaFold 3 inference."""

    ### Prepare inputs with modified seed.
    fold_input = folding_input.Input.from_json(self._test_input_json)
    fold_input = dataclasses.replace(fold_input, rng_seeds=[seed])

    output_dir = self.create_tempdir().full_path
    actual = run_alphafold.process_fold_input(
        fold_input,
        self._data_pipeline_config,
        run_alphafold.ModelRunner(
            config=self._model_config,
            device=jax.local_devices(backend='gpu')[0],
            model_dir=pathlib.Path(run_alphafold.MODEL_DIR.value),
        ),
        output_dir=output_dir,
        buckets=None if bucket is None else [bucket],
    )
    logging.info('finished get_inference_result')
    expected_model_cif_filename = f'{fold_input.sanitised_name()}_model.cif'
    expected_summary_confidences_filename = (
        f'{fold_input.sanitised_name()}_summary_confidences.json'
    )
    expected_confidences_filename = (
        f'{fold_input.sanitised_name()}_confidences.json'
    )
    expected_data_json_filename = f'{fold_input.sanitised_name()}_data.json'

    prefix = f'seed-{seed}'
    self.assertSameElements(
        os.listdir(output_dir),
        [
            # Subdirectories, one for each sample and one for embeddings.
            f'{prefix}_sample-0',
            f'{prefix}_sample-1',
            f'{prefix}_sample-2',
            f'{prefix}_sample-3',
            f'{prefix}_sample-4',
            f'{prefix}_embeddings',
            # Top ranking result.
            expected_confidences_filename,
            expected_model_cif_filename,
            expected_summary_confidences_filename,
            # Ranking scores for all samples.
            f'{fold_input.sanitised_name()}_ranking_scores.csv',
            # The input JSON defining the job.
            expected_data_json_filename,
            # The output terms of use.
            'TERMS_OF_USE.md',
        ],
    )

    for sample_index in range(5):
      sample_dir = os.path.join(output_dir, f'{prefix}_sample-{sample_index}')
      sample_prefix = (
          f'{fold_input.sanitised_name()}_seed-{seed}_sample-{sample_index}'
      )
      self.assertSameElements(
          os.listdir(sample_dir),
          [
              f'{sample_prefix}_confidences.json',
              f'{sample_prefix}_model.cif',
              f'{sample_prefix}_summary_confidences.json',
          ],
      )

    embeddings_dir = os.path.join(output_dir, f'{prefix}_embeddings')
    embeddings_filename = (
        f'{fold_input.sanitised_name()}_{prefix}_embeddings.npz'
    )
    self.assertSameElements(os.listdir(embeddings_dir), [embeddings_filename])

    with open(os.path.join(embeddings_dir, embeddings_filename), 'rb') as f:
      embeddings = np.load(f)
      self.assertSameElements(
          embeddings.keys(), ['single_embeddings', 'pair_embeddings']
      )
      # Ligand 7BU has 41 tokens.
      num_tokens = len(fold_input.protein_chains[0].sequence) + 41
      self.assertEqual(embeddings['single_embeddings'].shape, (num_tokens, 384))
      self.assertEqual(
          embeddings['pair_embeddings'].shape, (num_tokens, num_tokens, 128)
      )

    with open(os.path.join(output_dir, expected_data_json_filename), 'rt') as f:
      actual_input_json = json.load(f)

    self.assertEqual(
        actual_input_json['sequences'][0]['protein']['sequence'],
        fold_input.protein_chains[0].sequence,
    )
    self.assertSequenceEqual(
        actual_input_json['sequences'][1]['ligand']['ccdCodes'],
        fold_input.ligands[0].ccd_ids,
    )
    self.assertNotEmpty(
        actual_input_json['sequences'][0]['protein']['unpairedMsa']
    )
    self.assertNotEmpty(
        actual_input_json['sequences'][0]['protein']['pairedMsa']
    )
    self.assertIsNotNone(
        actual_input_json['sequences'][0]['protein']['templates']
    )

    ranking_scores_filename = (
        f'{fold_input.sanitised_name()}_ranking_scores.csv'
    )
    with open(os.path.join(output_dir, ranking_scores_filename), 'rt') as f:
      ranking_scores = list(csv.DictReader(f))

    self.assertLen(ranking_scores, 5)
    self.assertEqual([int(s['seed']) for s in ranking_scores], [seed] * 5)
    self.assertEqual(
        [int(s['sample']) for s in ranking_scores], [0, 1, 2, 3, 4]
    )

    # Ranking score should be between 0.66 and 0.76 for all samples.
    ranking_scores = [float(s['ranking_score']) for s in ranking_scores]
    scores_ok = [0.66 <= score <= 0.76 for score in ranking_scores]
    if not all(scores_ok):
      self.fail(f'{ranking_scores=} are not in expected range [0.66, 0.76]')

    with open(os.path.join(output_dir, 'TERMS_OF_USE.md'), 'rt') as f:
      actual_terms_of_use = f.read()
    self.assertStartsWith(
        actual_terms_of_use, '# ALPHAFOLD 3 OUTPUT TERMS OF USE'
    )

    bucket_label = 'default' if bucket is None else bucket
    output_filename = f'run_alphafold_test_output_bucket_{bucket_label}.pkl'

    # Convert to dict to enable simple serialization.
    actual_dict = [
        dict(
            seed=actual_inf.seed,
            inference_results=actual_inf.inference_results,
            full_fold_input=actual_inf.full_fold_input,
        )
        for actual_inf in actual
    ]
    with _output(output_filename) as (_, output):
      output.write(pickle.dumps(actual_dict))

    logging.info('Comparing inference results with expected values.')

    ### Assert that output is as expected.
    expected_dict = pickle.loads(
        (
            resources.ROOT
            / 'test_data'
            / 'alphafold_run_outputs'
            / output_filename
        ).read_bytes()
    )
    expected = [
        run_alphafold.ResultsForSeed(**expected_inf)
        for expected_inf in expected_dict
    ]

    actual_rmsds = []
    mask_proportions = []
    actual_masked_rmsds = []
    for actual_inf, expected_inf in zip(actual, expected, strict=True):
      for actual_inf, expected_inf in zip(
          actual_inf.inference_results,
          expected_inf.inference_results,
          strict=True,
      ):
        # Make sure the token chain IDs are the same as the input chain IDs.
        self.assertEqual(
            actual_inf.metadata['token_chain_ids'],
            ['P'] * len(fold_input.protein_chains[0].sequence) + ['LL'] * 41,
        )
        # All atom occupancies should be 1.0.
        np.testing.assert_array_equal(
            actual_inf.predicted_structure.atom_occupancy,
            [1.0] * actual_inf.predicted_structure.num_atoms,
        )
        actual_rmsds.append(
            alignment.rmsd_from_coords(
                decoy_coords=actual_inf.predicted_structure.coords,
                gt_coords=expected_inf.predicted_structure.coords,
            )
        )
        # Mask out atoms with b_factor < 80.0 (i.e. lower confidence regions).
        mask = actual_inf.predicted_structure.atom_b_factor > 80.0
        mask_proportions.append(
            np.sum(mask) / actual_inf.predicted_structure.num_atoms
        )
        actual_masked_rmsds.append(
            alignment.rmsd_from_coords(
                decoy_coords=actual_inf.predicted_structure.coords,
                gt_coords=expected_inf.predicted_structure.coords,
                include_idxs=mask,
            )
        )
    # 5tgy is stably predicted, samples should be all within 3.0 RMSD
    # regardless of seed, bucket, device type, etc.
    if any(rmsd > 3.0 for rmsd in actual_rmsds):
      self.fail(f'Full RMSD too high: {actual_rmsds=}')
    # Check proportion of atoms with b_factor > 80 is at least 70%.
    if any(prop < 0.7 for prop in mask_proportions):
      self.fail(f'Too many residues with low pLDDT: {mask_proportions=}')
    # Check masked RMSD is within tolerance (lower than full RMSD due to masking
    # of lower confidence regions).
    if any(rmsd > 1.4 for rmsd in actual_masked_rmsds):
      self.fail(f'Masked RMSD too high: {actual_masked_rmsds=}')


if __name__ == '__main__':
  absltest.main()
