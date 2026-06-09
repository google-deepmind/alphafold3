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

"""Tests the AlphaFold 3 data pipeline."""

import contextlib
import datetime
import difflib
import functools
import hashlib
import json
import os
import pathlib
import pickle
from typing import Any

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3 import structure
from alphafold3.common import folding_input
from alphafold3.common import resources
from alphafold3.common.testing import data as testing_data
from alphafold3.constants import chemical_components
from alphafold3.data import featurisation
from alphafold3.data import pipeline
from alphafold3.model.atom_layout import atom_layout
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


@functools.singledispatch
def _hash_data(x: Any, /) -> str:
  if x is None:
    return '<<None>>'
  return _hash_data(json.dumps(x).encode('utf-8'))


@_hash_data.register
def _(x: bytes, /) -> str:
  return hashlib.sha256(x).hexdigest()


@_hash_data.register
def _(x: jax.Array) -> str:
  return _hash_data(jax.device_get(x))


@_hash_data.register
def _(x: np.ndarray) -> str:
  if x.dtype == object:
    return ';'.join(map(_hash_data, x.ravel().tolist()))
  return _hash_data(x.tobytes())


@_hash_data.register
def _(_: structure.Structure) -> str:
  return '<<structure>>'


@_hash_data.register
def _(_: atom_layout.AtomLayout) -> str:
  return '<<atom-layout>>'


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


class DataPipelineTest(parameterized.TestCase):
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
                    'sequence': (
                        'SEFEKLRQTGDELVQAFQRLREIFDKGDDDSLEQVLEEIEELIQKHRQLFDNRQEAADTEAAKQGDQWVQLFQRFREAIDKGDKDSLEQLLEELEQALQKIRELAEKKN'
                    ),
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

  def compare_golden(self, result_path: str) -> None:
    filename = os.path.split(result_path)[1]
    golden_path = testing_data.Data(
        resources.ROOT / f'test_data/{filename}'
    ).path()
    with open(golden_path, 'r') as golden_file:
      golden_text = golden_file.read()
    with open(result_path, 'r') as result_file:
      result_text = result_file.read()

    diff = _generate_diff(result_text, golden_text)

    self.assertEqual(diff, "", f"Result differs from golden:\n{diff}")

  def test_config(self):
    model_config = run_alphafold.make_model_config()
    model_config_as_str = json.dumps(
        model_config.as_dict(), sort_keys=True, indent=2
    )
    with _output('model_config.json') as (result_path, output):
      output.write(model_config_as_str.encode('utf-8'))
    self.compare_golden(result_path)

  def test_featurisation(self):
    """Run featurisation and assert that the output is as expected."""
    fold_input = folding_input.Input.from_json(self._test_input_json)
    data_pipeline = pipeline.DataPipeline(self._data_pipeline_config)
    full_fold_input = data_pipeline.process(fold_input)
    featurised_example = featurisation.featurise_input(
        full_fold_input,
        ccd=chemical_components.Ccd(),
        buckets=None,
    )
    del featurised_example[0]['ref_pos']  # Depends on specific RDKit version.

    with _output('featurised_example.pkl') as (_, output):
      output.write(pickle.dumps(featurised_example))
    featurised_example = jax.tree_util.tree_map(_hash_data, featurised_example)
    with _output('featurised_example.json') as (result_path, output):
      output.write(
          json.dumps(featurised_example, sort_keys=True, indent=2).encode(
              'utf-8'
          )
      )
    self.compare_golden(result_path)

  def test_write_input_json(self):
    fold_input = folding_input.Input.from_json(self._test_input_json)
    output_dir = self.create_tempdir().full_path
    run_alphafold.write_fold_input_json(fold_input, output_dir)
    with open(
        os.path.join(output_dir, f'{fold_input.sanitised_name()}_data.json'),
        'rt',
    ) as f:
      actual_fold_input = folding_input.Input.from_json(f.read())

    self.assertEqual(actual_fold_input, fold_input)

  def test_process_fold_input_runs_only_data_pipeline(self):
    fold_input = folding_input.Input.from_json(self._test_input_json)
    output_dir = self.create_tempdir().full_path
    run_alphafold.process_fold_input(
        fold_input=fold_input,
        data_pipeline_config=self._data_pipeline_config,
        model_runner=None,
        output_dir=output_dir,
    )
    with open(
        os.path.join(output_dir, f'{fold_input.sanitised_name()}_data.json'),
        'rt',
    ) as f:
      actual_fold_input = folding_input.Input.from_json(f.read())

    featurisation.validate_fold_input(actual_fold_input)

  @parameterized.product(num_db_dirs=tuple(range(1, 3)))
  def test_replace_db_dir(self, num_db_dirs: int) -> None:
    """Test that the db_dir is replaced correctly."""
    db_dirs = [pathlib.Path(self.create_tempdir()) for _ in range(num_db_dirs)]
    db_dirs_posix = [db_dir.as_posix() for db_dir in db_dirs]

    for i, db_dir in enumerate(db_dirs):
      for j in range(i + 1):
        (db_dir / f'filename{j}.txt').write_text(f'hello world {i}')

    for i in range(num_db_dirs):
      self.assertEqual(
          pathlib.Path(
              run_alphafold.replace_db_dir(
                  f'${{DB_DIR}}/filename{i}.txt', db_dirs_posix
              )
          ).read_text(),
          f'hello world {i}',
      )
    with self.assertRaises(FileNotFoundError):
      run_alphafold.replace_db_dir(
          f'${{DB_DIR}}/filename{num_db_dirs}.txt', db_dirs_posix
      )


if __name__ == '__main__':
  absltest.main()
