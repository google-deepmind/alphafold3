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
from alphafold3.constants import residue_names
from alphafold3.model import data3
import numpy as np


class MsaTest(parameterized.TestCase):

  def test_get_profile_features(self):
    num_seqs = 12
    num_res = 3
    num_restypes = residue_names.POLYMER_TYPES_NUM_WITH_UNKNOWN_AND_GAP

    # Residue 0 is a 50-50 split between 0 and 1
    # Residue 1 is a 75-25 split between 0 and num_restypes - 1
    msa = np.zeros((num_seqs, num_res), dtype=np.int32)
    msa[:6, 0] = 1
    msa[:3, 1] = num_restypes - 1

    # Residue 0 has a deletion half the time. There are no other deletions.
    deletion_matrix = np.zeros((num_seqs, num_res), dtype=np.int32)
    deletion_matrix[:6, 0] = 1

    feats = data3.get_profile_features(msa, deletion_matrix)

    self.assertEqual((num_res, num_restypes), feats['profile'].shape)
    self.assertEqual((num_res,), feats['deletion_mean'].shape)

    expected_profile = np.zeros((num_res, num_restypes))
    expected_profile[0][0] = 0.5
    expected_profile[0][1] = 0.5
    expected_profile[1][0] = 0.75
    expected_profile[1][-1] = 0.25
    expected_profile[2][0] = 1.0
    np.testing.assert_allclose(expected_profile, feats['profile'])

    expected_deletion_m = np.zeros((num_res,))
    expected_deletion_m[0] = 0.5
    np.testing.assert_allclose(expected_deletion_m, feats['deletion_mean'])


if __name__ == '__main__':
  absltest.main()
