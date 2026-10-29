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
from alphafold3.cpp import msa_profile
import numpy as np


_RNG = np.random.default_rng(seed=42)


def _msa_profile(msa, num_residue_types):
  """The original Python implementation we are replacing by the C++ version."""
  _, batch = msa.shape
  profile = np.empty((batch, num_residue_types), dtype=np.float32)
  for i in range(batch):
    profile[i] = np.bincount(msa[:, i], minlength=num_residue_types).astype(
        np.float32
    )
    profile[i] /= np.sum(profile[i])
  return profile


class MsaProfileTest(parameterized.TestCase):

  @parameterized.parameters(
      (np.ones((1000, 500), dtype=np.int32), 10),
      (np.zeros((1000, 500), dtype=np.int32), 10),
      # Protein + DNA + RNA.
      (_RNG.integers(low=0, high=30, size=(1000, 500), dtype=np.int32), 30),
      # Protein only.
      (_RNG.integers(low=0, high=22, size=(1000, 500), dtype=np.int32), 30),
      # DNA + RNA only.
      (_RNG.integers(low=22, high=30, size=(1000, 500), dtype=np.int32), 30),
      # MSAs with depth 1.
      (np.ones((1, 500), dtype=np.int32), 30),
      (_RNG.integers(low=0, high=30, size=(1, 500), dtype=np.int32), 30),
  )
  def test_msa_profile_against_python(self, msa, num_residue_types):
    cpp_profile = msa_profile.compute_msa_profile(msa, num_residue_types)
    python_profile = _msa_profile(msa, num_residue_types)
    np.testing.assert_array_almost_equal(python_profile, cpp_profile)

  def test_msa_profile_value(self):
    msa = np.array(
        [[1, 2, 3, 4], [0, 2, 3, 4], [0, 0, 3, 4], [0, 0, 0, 4]], dtype=np.int32
    )
    expected_profile = np.array(
        [
            [0.75, 0.25, 0.0, 0.0, 0.0],
            [0.5, 0.0, 0.5, 0.0, 0.0],
            [0.25, 0.0, 0.0, 0.75, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )
    profile = msa_profile.compute_msa_profile(msa, num_residue_types=5)
    np.testing.assert_array_almost_equal(expected_profile, profile)

  def test_bad_shape(self):
    with self.assertRaisesRegex(
        ValueError,
        r'The MSA must be rectangular, got 1-dimensional MSA array\.',
    ):
      msa_profile.compute_msa_profile(np.array([0], dtype=np.int32), 5)

  def test_empty(self):
    with self.assertRaisesRegex(ValueError, r'The MSA must be non-empty\.'):
      msa_profile.compute_msa_profile(np.array([], dtype=np.int32), 5)

  @parameterized.parameters((5,), (10,))
  def test_bad_num_residue_types(self, num_residue_types):
    arr = np.array([[10]], dtype=np.int32)
    with self.assertRaises(ValueError):
      _msa_profile(arr, num_residue_types)

    with self.assertRaises(ValueError):
      msa_profile.compute_msa_profile(arr, num_residue_types)

  def test_negative_residue_codes(self):
    arr = np.array([[-1]], dtype=np.int32)
    with self.assertRaises(ValueError):
      _msa_profile(arr, num_residue_types=2)

    with self.assertRaises(ValueError):
      msa_profile.compute_msa_profile(arr, num_residue_types=2)


if __name__ == '__main__':
  absltest.main()
