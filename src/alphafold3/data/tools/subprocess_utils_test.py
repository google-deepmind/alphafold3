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
import subprocess

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.data.tools import subprocess_utils


import shutil

_JACKHMMER_BINARY_PATH = shutil.which('jackhmmer')


class LaunchUtilsTest(absltest.TestCase):

  def test_create_query_fasta_file(self):
    tmp_dir = self.create_tempdir().full_path
    query_file_path = os.path.join(tmp_dir, 'query.fasta')
    subprocess_utils.create_query_fasta_file(
        sequence='AAABBB', path=query_file_path, linewidth=3
    )

    with open(query_file_path) as f:
      self.assertEqual(f.read(), '>query\nAAA\nBBB\n')


class RunSubprocessTest(parameterized.TestCase):

  def test_successful_subprocess(self):
    process = subprocess_utils.run(cmd=['echo', '0'], cmd_name='ok')
    self.assertEqual(process.returncode, 0)
    self.assertEqual(process.stdout, '0\n')
    self.assertEmpty(process.stderr)

  def test_subprocess_timeout(self):
    with self.assertRaises(subprocess.TimeoutExpired):
      subprocess_utils.run(cmd=['sleep', '3'], cmd_name='timeout', timeout=1)

  @parameterized.product(log_stdout=(False, True), log_stderr=(False, True))
  def test_logging(self, log_stdout, log_stderr):
    cmd = 'printf "THIS_IS_STDOUT"; printf "THIS_IS_STDERR" >&2; exit 0'
    tmp_dir = self.create_tempdir().full_path
    script_path = os.path.join(tmp_dir, 'script.sh')
    with open(script_path, 'w') as f:
      f.write(cmd)

    with self.assertLogs() as log_output:
      subprocess_utils.run(
          cmd=['bash', script_path],
          cmd_name='ok',
          log_stdout=log_stdout,
          log_stderr=log_stderr,
      )
    self.assertEqual(
        log_output[0][0].message, f'Launching subprocess "bash {script_path}"'
    )
    self.assertStartsWith(log_output[0][1].message, 'Finished ok in ')
    self.assertEndsWith(log_output[0][1].message, 'seconds')
    if log_stdout:
      self.assertEqual(log_output[0][2].message, 'ok stdout:\nTHIS_IS_STDOUT')
    if log_stderr:
      stderr_msg_index = 3 if log_stdout else 2
      self.assertEqual(
          log_output[0][stderr_msg_index].message, 'ok stderr:\nTHIS_IS_STDERR'
      )

  def test_failed_subprocess(self):
    cmd = 'printf "OKAY"; printf "FAILED" >&2; exit 1'
    tmp_dir = self.create_tempdir().full_path
    script_path = os.path.join(tmp_dir, 'script.sh')
    with open(script_path, 'w') as f:
      f.write(cmd)
    with self.assertRaisesRegex(
        RuntimeError, 'failing_cmd failed\nstdout:\nOKAY\n\nstderr:\nFAILED'
    ):
      subprocess_utils.run(cmd=['bash', script_path], cmd_name='failing_cmd')

  def test_jackhmmer_seq_limit_unsupported(self):
    # Test with the ls command which does not support the --seq-limit flag.
    self.assertFalse(subprocess_utils.jackhmmer_seq_limit_supported('ls'))

  def test_jackhmmer_seq_limit_bad_path(self):
    with self.assertRaises(FileNotFoundError):
      subprocess_utils.jackhmmer_seq_limit_supported('/non_existent_binary')


if __name__ == '__main__':
  absltest.main()
