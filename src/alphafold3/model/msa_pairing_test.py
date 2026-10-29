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

import copy
from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.model import data_constants
from alphafold3.model import msa_pairing
import numpy as np


class PairedAlignmentTest(parameterized.TestCase):

  def test_deduplicate_unpaired_sequences(self):
    chains = [
        {
            'msa_all_seq': np.ones((1, 5), dtype=np.int32),
            'msa': np.concatenate((
                np.ones((2, 5), dtype=np.int32),
                2 * np.ones((1, 5), dtype=np.int32),
            )),
        },
    ]
    result = msa_pairing.deduplicate_unpaired_sequences(chains)
    np.testing.assert_equal(
        result[0]['msa'], 2 * np.ones((1, 5), dtype=np.int32)
    )

    np.testing.assert_equal(
        result[0]['msa_all_seq'], np.ones((1, 5), dtype=np.int32)
    )

  def test_create_paired_features(self):
    np_chains = [
        {
            'chain_id': 'A',
            'msa_species_identifiers_all_seq': np.array([b'', b'S2', b'S3']),
            'num_alignments_all_seq': 3,
            'msa_all_seq': np.array([[1, 1], [2, 2], [3, 3]]),
            'deletion_matrix_all_seq': np.array([[0, 0], [0, 0], [0, 0]]),
        },
        {
            'chain_id': 'B',
            'msa_species_identifiers_all_seq': np.array([b'', b'S3']),
            'num_alignments_all_seq': 2,
            'msa_all_seq': np.array([[4, 4, 4], [5, 5, 5]]),
            'deletion_matrix_all_seq': np.array([[0, 0, 0], [0, 0, 0]]),
        },
    ]

    np_chains_paired = msa_pairing.create_paired_features(
        chains=np_chains,  # pyrefly: ignore[bad-argument-type]
        max_paired_sequences=128,
        nonempty_chain_ids={'A', 'B'},
        max_hits_per_species=600,
    )
    self.assertEqual(np_chains_paired[0]['num_alignments_all_seq'], 2)
    self.assertEqual(np_chains_paired[1]['num_alignments_all_seq'], 2)
    np.testing.assert_array_equal(
        np_chains_paired[0]['msa_all_seq'], np.array([[1, 1], [3, 3]])
    )
    np.testing.assert_array_equal(
        np_chains_paired[1]['msa_all_seq'], np.array([[4, 4, 4], [5, 5, 5]])
    )

    paired_feats = (
        ['num_alignments']
        + list(data_constants.NUM_SEQ_NUM_RES_MSA_FEATURES)
        + list(data_constants.NUM_SEQ_MSA_FEATURES)
    )
    for k in np_chains_paired[0]:
      k_split = k.split('_all_seq')[0]
      if k_split not in paired_feats:
        np.testing.assert_equal(
            np_chains_paired[0][k], np_chains[0][k], err_msg=k
        )
        np.testing.assert_equal(
            np_chains_paired[1][k], np_chains[1][k], err_msg=k
        )

  def test_create_paired_features_large(self):
    np_chains = [
        {
            'chain_id': 'A',
            'msa_species_identifiers_all_seq': np.array(
                [b'', b'S2'] + [b'S3'] * 900
            ),
            'num_alignments_all_seq': 3,
            'msa_all_seq': np.array([[1, 1]] * 902),
            'deletion_matrix_all_seq': np.array([[0, 0]] * 902),
        },
        {
            'chain_id': 'B',
            'msa_species_identifiers_all_seq': np.array([b''] + [b'S3'] * 800),
            'num_alignments_all_seq': 2,
            'msa_all_seq': np.array([[4, 4, 4]] * 801),
            'deletion_matrix_all_seq': np.array([[0, 0, 0]] * 801),
        },
    ]

    np_chains_paired = msa_pairing.create_paired_features(
        chains=np_chains,  # pyrefly: ignore[bad-argument-type]
        max_paired_sequences=128,
        max_hits_per_species=20,
        nonempty_chain_ids={'A', 'B'},
    )
    self.assertEqual(np_chains_paired[0]['num_alignments_all_seq'], 21)

  def test_create_paired_features_multichain(self):
    np_chains = [
        {
            'chain_id': 'A',
            'msa_species_identifiers_all_seq': np.array(
                [b'', b'S2', b'S3', b'S4', b'S2']
            ),
            'num_alignments_all_seq': 5,
            'msa_all_seq': np.array([[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]]),
            'deletion_matrix_all_seq': np.array([[0, 0]] * 5),
        },
        {
            'chain_id': 'B',
            'msa_species_identifiers_all_seq': np.array([b'', b'S4']),
            'num_alignments_all_seq': 2,
            'msa_all_seq': np.array([[6, 6, 6], [7, 7, 7]]),
            'deletion_matrix_all_seq': np.array([[0, 0, 0]] * 5),
        },
        {
            'chain_id': 'C',
            'msa_species_identifiers_all_seq': np.array(
                [b'', b'S2', b'S2', b'S4']
            ),
            'num_alignments_all_seq': 4,
            'msa_all_seq': np.array([[8, 8], [9, 9], [10, 10], [11, 11]]),
            'deletion_matrix_all_seq': np.array([[0, 0]] * 4),
        },
    ]

    np_chains_paired = msa_pairing.create_paired_features(
        chains=np_chains,  # pyrefly: ignore[bad-argument-type]
        max_paired_sequences=128,
        max_hits_per_species=20,
        nonempty_chain_ids={'A', 'B', 'C'},
    )
    self.assertEqual(np_chains_paired[0]['num_alignments_all_seq'], 4)
    np.testing.assert_array_equal(
        np_chains_paired[0]['msa_all_seq'],
        np.array([[1, 1], [4, 4], [2, 2], [5, 5]]),
    )

    gap_val = data_constants.MSA_GAP_IDX
    np.testing.assert_array_equal(
        np_chains_paired[1]['msa_all_seq'],
        np.array([
            [6, 6, 6],
            [7, 7, 7],
            [gap_val, gap_val, gap_val],
            [gap_val, gap_val, gap_val],
        ]),
    )
    np.testing.assert_array_equal(
        np_chains_paired[2]['msa_all_seq'],
        np.array([
            [8, 8],
            [11, 11],
            [9, 9],
            [10, 10],
        ]),
    )

  def test_create_paired_features_empty(self):
    np_chains = [
        {
            'chain_id': 'A',
            'msa_species_identifiers_all_seq': np.array([b'']),
            'num_alignments_all_seq': 1,
            'msa_all_seq': np.array([[1, 1]]),
            'deletion_matrix_all_seq': np.array([[0, 0]]),
        },
        {
            'chain_id': 'B',
            'msa_species_identifiers_all_seq': np.array([b'']),
            'num_alignments_all_seq': 1,
            'msa_all_seq': np.array([[6, 6, 6]]),
            'deletion_matrix_all_seq': np.array([[0, 0, 0]]),
        },
    ]

    np_chains_paired = msa_pairing.create_paired_features(
        chains=np_chains,  # pyrefly: ignore[bad-argument-type]
        max_paired_sequences=128,
        max_hits_per_species=20,
        nonempty_chain_ids={'A', 'B', 'C'},
    )
    self.assertEqual(np_chains_paired[0]['num_alignments_all_seq'], 1)

  def test_make_msa_crop_sizes_without_pairing(self):
    # Fits.
    crop_sizes = msa_pairing.choose_paired_unpaired_msa_crop_sizes(
        paired_msa=None,
        unpaired_msa=np.ones((5, 13)),
        total_msa_crop_size=11,
        max_paired_sequences=5,
    )
    self.assertTupleEqual(crop_sizes, (5, None))

    # Doesn't fit.
    crop_sizes = msa_pairing.choose_paired_unpaired_msa_crop_sizes(
        paired_msa=None,
        unpaired_msa=np.ones((11, 13)),
        total_msa_crop_size=7,
        max_paired_sequences=3,
    )
    self.assertTupleEqual(crop_sizes, (7, None))

  @parameterized.parameters(
      [11, 2, 1, (5, 2)],
      [11, 4, 1, (4, 3)],
      [2, 2, 1, (2, 2)],
      [2, 4, 1, (2, 3)],
  )
  def test_make_msa_crop_sizes_with_pairing(
      self,
      num_alignments_reduced,
      num_alignments_unreduced,
      num_non_gapped,
      expected,
  ):
    num_gapped = num_alignments_unreduced - num_non_gapped
    paired_msa = np.array(
        [[1, 2, 3]] * num_non_gapped
        + [[data_constants.MSA_GAP_IDX] * 3] * num_gapped
    )
    unpaired_msa = np.ones((num_alignments_reduced, 13))

    crop_sizes = msa_pairing.choose_paired_unpaired_msa_crop_sizes(
        paired_msa=paired_msa,
        unpaired_msa=unpaired_msa,
        total_msa_crop_size=7,
        max_paired_sequences=3,
    )
    self.assertTupleEqual(crop_sizes, expected)

  def test_remove_all_gapped_rows_from_all_seqs(self):
    """Tests that all gapped rows are correctly removed from all_seq feats."""

    num_seq = 10
    num_res = 4
    row_indices = np.arange(num_seq)
    chain1 = {
        'msa_all_seq': (
            np.ones((num_seq, num_res), dtype=np.int32) * row_indices[:, None]
        ),
        'asym_id': np.ones(num_res, dtype=np.int32),
    }
    chain2 = {
        'msa_all_seq': (
            np.ones((num_seq, num_res), dtype=np.int32) * row_indices[:, None]
        ),
        'asym_id': 2 * np.ones(num_res, dtype=np.int32),
    }
    chain3 = {
        'msa_all_seq': (
            np.ones((num_seq, num_res), dtype=np.int32) * row_indices[:, None]
        ),
        'asym_id': 3 * np.ones(num_res, dtype=np.int32),
    }

    # Set all gaps for chains 2 and 3 row index 3.
    chain2['msa_all_seq'][[1, 3, 4]] = data_constants.MSA_GAP_IDX
    chain3['msa_all_seq'][[3, 5]] = data_constants.MSA_GAP_IDX
    orig_chains_list = [chain1, chain2, chain3]

    with self.subTest('NoDeduplicationWhenNoAllGapsRows'):
      chains_list = copy.deepcopy(orig_chains_list)
      chains_list = msa_pairing.remove_all_gapped_rows_from_all_seqs(
          chains_list, [1, 2, 3]
      )
      for orig_chain, chain in zip(orig_chains_list, chains_list):
        exp_msa_all_seq = orig_chain['msa_all_seq']
        np.testing.assert_array_equal(exp_msa_all_seq, chain['msa_all_seq'])

    with self.subTest('DeduplicationWhenAllGapsRows'):
      chains_list = copy.deepcopy(orig_chains_list)
      chains_list = msa_pairing.remove_all_gapped_rows_from_all_seqs(
          chains_list, [2, 3]
      )
      # row index 3 should have been removed from all chains
      for orig_chain, chain in zip(orig_chains_list, chains_list):
        exp_msa_all_seq = orig_chain['msa_all_seq'][row_indices != 3]
        np.testing.assert_array_equal(exp_msa_all_seq, chain['msa_all_seq'])

  def test_pairing_preserves_within_species_order(self):
    # NumPy's default argsort uses insertion sort up to 16 elements, so 17 rows
    # is the minimum size where an unstable sort reorders equal species IDs.
    chains = [
        {
            'chain_id': 'A',
            'msa_species_identifiers_all_seq': np.array([b''] + [b'S'] * 16),
            'msa_all_seq': np.arange(17, dtype=np.int32)[:, None],
            'deletion_matrix_all_seq': np.zeros((17, 1), dtype=np.int32),
        },
        {
            'chain_id': 'B',
            'msa_species_identifiers_all_seq': np.array([b'', b'S', b'S']),
            'msa_all_seq': np.arange(3, dtype=np.int32)[:, None],
            'deletion_matrix_all_seq': np.zeros((3, 1), dtype=np.int32),
        },
    ]

    paired = msa_pairing.create_paired_features(
        chains=chains,  # pyrefly: ignore[bad-argument-type]
        max_paired_sequences=10,
        nonempty_chain_ids={'A', 'B'},
        max_hits_per_species=10,
    )

    np.testing.assert_array_equal(paired[0]['msa_all_seq'], [[0], [1], [2]])
    np.testing.assert_array_equal(paired[1]['msa_all_seq'], [[0], [1], [2]])


if __name__ == '__main__':
  absltest.main()
