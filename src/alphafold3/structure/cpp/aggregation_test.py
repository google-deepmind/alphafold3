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

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.cpp import aggregation
import numpy as np


class GroupbyTest(parameterized.TestCase):

  @parameterized.parameters([
      ([], {}),
      ([1], {1: [0]}),
      ([-1], {-1: [0]}),
      ([10, 10, 10], {10: [0, 1, 2]}),
      ([1, 1, 2, 3, 3, 1, 1], {1: [0, 1, 5, 6], 2: [2], 3: [3, 4]}),
  ])
  def test_get_group_indices(self, array, expected):
    self.assertEqual(
        aggregation.indices_grouped_by_value(np.array(array, dtype=np.int64)),
        expected,
    )


if __name__ == '__main__':
  absltest.main()
