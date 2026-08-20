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

"""Tests for data parsing helpers and the subprocess runner."""

from absl.testing import absltest
from alphafold3.data.tools import subprocess_utils
import sys
import time


class SubprocessUtilsTest(absltest.TestCase):

  def test_run_times_out_and_kills_hung_process(self):
    start = time.time()
    with self.assertRaises(RuntimeError) as cm:
      subprocess_utils.run(
          [sys.executable, '-c', 'import time; time.sleep(30)'],
          cmd_name='sleep',
          timeout=1,
      )
    elapsed = time.time() - start
    self.assertLess(elapsed, 10, 'process was not killed promptly')
    self.assertIn('timed out after 1 seconds', str(cm.exception))


if __name__ == '__main__':
  absltest.main()