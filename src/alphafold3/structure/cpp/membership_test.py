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
from alphafold3.cpp import membership
import numpy as np


class MembershipTest(parameterized.TestCase):

  @parameterized.parameters(
      ([], set(), False, []),
      ([], set([1, 2]), False, []),
      ([1, 2, 3], set([]), False, [False, False, False]),
      ([1, 2, 3], set([2]), False, [False, True, False]),
      ([1, 2, 3], set([1, 2]), False, [True, True, False]),
      ([1, 2, 3], set([1, 2, 3]), False, [True, True, True]),
      ([1, 2, 2], set([1, 2, 3]), False, [True, True, True]),
      ([2, 1, 3], set([1, 2, 3]), False, [True, True, True]),
      ([], set(), True, []),
      ([], set([1, 2]), True, []),
      ([1, 2, 3], set([]), True, [True, True, True]),
      ([1, 2, 3], set([2]), True, [True, False, True]),
      ([1, 2, 3], set([1, 2]), True, [False, False, True]),
      ([1, 2, 3], set([1, 2, 3]), True, [False, False, False]),
      ([1, 2, 2], set([1, 2, 3]), True, [False, False, False]),
      ([2, 1, 3], set([1, 2, 3]), True, [False, False, False]),
  )
  def test_isin(self, array, test_elements, invert, exp_result):
    result = membership.isin(
        np.array(array, dtype=np.int64), test_elements, invert=invert
    )
    self.assertSequenceEqual(result.tolist(), exp_result)

  def test_isin_multi_dim(self):
    result = membership.isin(
        np.array([[1, 2, 0], [3, 4, 0]], dtype=np.int64), set([1, 4])
    )
    self.assertSequenceEqual(
        result.tolist(), [[True, False, False], [False, True, False]]
    )


if __name__ == '__main__':
  absltest.main()
