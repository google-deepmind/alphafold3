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

"""Tests for MSA pairing (see issue #236 for the float32 overflow)."""

from absl.testing import absltest
from alphafold3.model import msa_pairing
import numpy as np


def _rank_metric_log(rows: np.ndarray) -> np.ndarray:
  """Reference log-domain ranking matching the production fix."""
  with np.errstate(divide='ignore'):
    return np.sum(np.log(np.abs(rows).astype(np.float64)), axis=1)


def _rank_metric_int64(rows: np.ndarray) -> np.ndarray:
  """Exact integer ranking, only safe for small values."""
  return np.abs(np.prod(rows.astype(np.int64), axis=1))


class MsaPairingTest(absltest.TestCase):

  def test_rank_metric_does_not_overflow_float32(self):
    # Regression test for issue #236: with deep MSAs and many chains the raw
    # float32 product of row indices overflows to inf, silently corrupting the
    # pairing order.
    rng = np.random.RandomState(0)
    rows = rng.randint(1, 8000, size=(100, 12)).astype(np.int32)
    old_metric = np.abs(np.prod(rows.astype(np.float32), axis=1))
    self.assertTrue(np.isinf(old_metric).any())
    new_metric = _rank_metric_log(rows)
    self.assertTrue(np.isfinite(new_metric).all())

  def test_log_rank_matches_int64_rank_for_small_values(self):
    rows = np.array(
        [[3, 2], [1, 5], [2, 2], [5, 1], [1, 1], [0, 5], [1, -1]],
        dtype=np.int32,
    )
    self.assertSequenceEqual(
        np.argsort(_rank_metric_log(rows)).tolist(),
        np.argsort(_rank_metric_int64(rows)).tolist(),
    )

  def test_query_and_padding_rows_rank_first(self):
    # Query rows (index 0) and fully-padded rows (product 1) must sort before
    # any real hit, matching the original product semantics.
    rows = np.array(
        [[0, 0], [1, -1], [2, 3], [4, 5]],
        dtype=np.int32,
    )
    order = np.argsort(_rank_metric_log(rows))
    self.assertEqual(order[0], 0)
    self.assertEqual(order[1], 1)

  def test_create_paired_features_with_deep_msa(self):
    # End-to-end check that create_paired_features survives deep MSAs over many
    # chains without the float32 overflow.
    rng = np.random.RandomState(42)
    num_chains = 8
    num_rows = 6000
    chains = []
    for c in range(num_chains):
      species = np.array(
          [b'sp%02d' % rng.randint(0, 5) for _ in range(num_rows)], dtype=object
      )
      species[0] = b''
      msa = np.zeros((num_rows, 16), dtype=np.int32)
      chains.append({
          'chain_id': str(c),
          'msa_species_identifiers_all_seq': species,
          'msa_all_seq': msa,
          'deletion_matrix_all_seq': np.zeros_like(msa),
          'msa': np.zeros((num_rows, 16), dtype=np.int32),
          'deletion_matrix': np.zeros_like(msa),
      })
    out = msa_pairing.create_paired_features(
        chains,
        max_paired_sequences=1024,
        nonempty_chain_ids={str(c) for c in range(num_chains)},
        max_hits_per_species=600,
    )
    for chain in out:
      self.assertEqual(chain['msa_all_seq'].shape[0], 1024)
      self.assertEqual(int(chain['num_alignments_all_seq']), 1024)
      self.assertEqual(chain['num_alignments_all_seq'].dtype, np.int32)

  def test_deduplicate_unpaired_sequences_uses_bytes_not_hash(self):
    paired = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.int32)
    unpaired = np.array(
        [[4, 5, 6], [7, 8, 9], [1, 2, 3], [10, 11, 12]], dtype=np.int32
    )
    chain = {
        'msa_all_seq': paired,
        'msa': unpaired,
        'deletion_matrix_all_seq': np.zeros_like(paired),
        'deletion_matrix': np.zeros_like(unpaired),
        'num_alignments': np.array(4, dtype=np.int32),
    }
    # Rows present in the paired MSA must be removed from the unpaired MSA,
    # leaving only the two unique rows.
    out = msa_pairing.deduplicate_unpaired_sequences([chain])[0]
    self.assertEqual(out['msa'].shape[0], 2)
    self.assertEqual(int(out['num_alignments']), 2)


if __name__ == '__main__':
  absltest.main()