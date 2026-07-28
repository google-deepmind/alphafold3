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

"""Unit tests for data-pipeline reliability fixes."""

import io
import os
import pathlib
import tarfile
import tempfile
from unittest import mock

from absl.testing import absltest

from alphafold3.data import structure_stores
from alphafold3.data.tools import hmmsearch
from alphafold3.data.tools import subprocess_utils


class CheckBinaryExistsTest(absltest.TestCase):

  def test_rejects_none_path(self):
    with self.assertRaisesRegex(RuntimeError, 'not set'):
      subprocess_utils.check_binary_exists(None, 'Jackhmmer')

  def test_rejects_empty_path(self):
    with self.assertRaisesRegex(RuntimeError, 'not set'):
      subprocess_utils.check_binary_exists('', 'Jackhmmer')

  def test_rejects_missing_path(self):
    with self.assertRaisesRegex(RuntimeError, 'not found'):
      subprocess_utils.check_binary_exists(
          '/path/does/not/exist/jackhmmer', 'Jackhmmer'
      )

  def test_accepts_existing_path(self):
    with tempfile.NamedTemporaryFile() as tmp:
      subprocess_utils.check_binary_exists(tmp.name, 'FakeBinary')


class StructureStoreLifecycleTest(absltest.TestCase):

  def _make_cif_tar(self, tar_path: pathlib.Path) -> None:
    cif_bytes = b'data_test\n_entry.id test\n'
    with tarfile.open(tar_path, 'w') as tar:
      info = tarfile.TarInfo(name='abcd.cif')
      info.size = len(cif_bytes)
      tar.addfile(info, io.BytesIO(cif_bytes))

  def test_context_manager_closes_tar(self):
    with tempfile.TemporaryDirectory() as tmp_dir:
      tar_path = pathlib.Path(tmp_dir) / 'pdb.tar'
      self._make_cif_tar(tar_path)
      with structure_stores.StructureStore(tar_path) as store:
        self.assertEqual(store.get_mmcif_str('abcd'), 'data_test\n_entry.id test\n')
      with self.assertRaisesRegex(IOError, 'closed'):
        store.get_mmcif_str('abcd')

  def test_mapping_store_roundtrip(self):
    store = structure_stores.StructureStore({'abcd': 'mmcif-data'})
    self.assertEqual(store.get_mmcif_str('abcd'), 'mmcif-data')
    store.close()


class HmmsearchCpuFlagTest(absltest.TestCase):

  def test_uses_configured_n_cpu(self):
    with tempfile.TemporaryDirectory() as tmp_dir:
      db_path = os.path.join(tmp_dir, 'db.fasta')
      pathlib.Path(db_path).write_text('>db\nAAAA\n')
      binary_path = os.path.join(tmp_dir, 'hmmsearch')
      hmmbuild_path = os.path.join(tmp_dir, 'hmmbuild')
      pathlib.Path(binary_path).write_text('#!/bin/sh\n')
      pathlib.Path(hmmbuild_path).write_text('#!/bin/sh\n')
      os.chmod(binary_path, 0o755)
      os.chmod(hmmbuild_path, 0o755)

      runner = hmmsearch.Hmmsearch(
          binary_path=binary_path,
          hmmbuild_binary_path=hmmbuild_path,
          database_path=db_path,
          n_cpu=3,
      )

      def _fake_run(*, cmd, **kwargs):
        del kwargs
        # Create the Stockholm output path passed via -A.
        a_idx = cmd.index('-A')
        pathlib.Path(cmd[a_idx + 1]).write_text('# STOCKHOLM 1.0\n//\n')
        return mock.Mock()

      with mock.patch.object(
          hmmsearch.subprocess_utils, 'run', side_effect=_fake_run
      ) as mock_run:
        runner.query_with_hmm('HMM')

      cmd = mock_run.call_args.kwargs['cmd']
      cpu_idx = cmd.index('--cpu')
      self.assertEqual(cmd[cpu_idx + 1], '3')

  def test_rejects_invalid_n_cpu(self):
    with tempfile.TemporaryDirectory() as tmp_dir:
      db_path = os.path.join(tmp_dir, 'db.fasta')
      pathlib.Path(db_path).write_text('>db\nAAAA\n')
      binary_path = os.path.join(tmp_dir, 'hmmsearch')
      hmmbuild_path = os.path.join(tmp_dir, 'hmmbuild')
      pathlib.Path(binary_path).write_text('#!/bin/sh\n')
      pathlib.Path(hmmbuild_path).write_text('#!/bin/sh\n')
      with self.assertRaisesRegex(ValueError, 'n_cpu'):
        hmmsearch.Hmmsearch(
            binary_path=binary_path,
            hmmbuild_binary_path=hmmbuild_path,
            database_path=db_path,
            n_cpu=0,
        )


if __name__ == '__main__':
  absltest.main()
