# Copyright 2026 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for pairing MSA hits in their original within-species order."""

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.model import msa_pairing
import numpy as np


class MsaPairingTest(parameterized.TestCase):

  @parameterized.named_parameters(
      ('unequal_hit_counts', 20),
      ('per_species_cap', 3),
  )
  def test_pairing_preserves_within_species_order(self, max_hits_per_species):
    chains = []
    for chain_id, species in (
        ('A', [b''] + [b'species_a', b'species_b'] * 20),
        ('B', [b''] + [b'species_b', b'species_a'] * 15),
    ):
      chains.append({
          'chain_id': chain_id,
          'msa_species_identifiers_all_seq': np.array(species, dtype=object),
          'msa_all_seq': np.zeros((len(species), 1), dtype=np.int32),
          # Use deletion counts to identify the original rows after pairing.
          'deletion_matrix_all_seq': np.arange(
              len(species), dtype=np.int32
          )[:, None],
      })

    paired_chains = msa_pairing.create_paired_features(
        chains,
        max_paired_sequences=100,
        nonempty_chain_ids={'A', 'B'},
        max_hits_per_species=max_hits_per_species,
    )
    actual_pairs = np.column_stack([
        chain['deletion_matrix_all_seq'][:, 0] for chain in paired_chains
    ])
    expected_pairs = [(0, 0)]
    for i in range(min(15, max_hits_per_species)):
      expected_pairs.extend([
          (2 * i + 1, 2 * i + 2), (2 * i + 2, 2 * i + 1)
      ])

    np.testing.assert_array_equal(actual_pairs[0], [0, 0])
    # Global ranking may reorder pairs with equal rank, so ignore their order.
    self.assertCountEqual(map(tuple, actual_pairs), expected_pairs)


if __name__ == '__main__':
  absltest.main()
