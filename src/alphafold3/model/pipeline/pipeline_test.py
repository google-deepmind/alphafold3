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

from collections.abc import Sequence

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.common import folding_input
from alphafold3.constants import chemical_components
from alphafold3.model.pipeline import pipeline
import numpy as np


class PipelineTest(parameterized.TestCase):

  @parameterized.product(
      [
          dict(buckets=None, num_tokens=100, exp_bucket=100),
          dict(buckets=[100], num_tokens=100, exp_bucket=100),
          dict(buckets=[100], num_tokens=25, exp_bucket=100),
          dict(buckets=[100, 200, 300], num_tokens=25, exp_bucket=100),
          dict(buckets=[100, 200, 300], num_tokens=100, exp_bucket=100),
          dict(buckets=[100, 200, 300], num_tokens=101, exp_bucket=200),
          dict(buckets=[100, 200, 300], num_tokens=250, exp_bucket=300),
          dict(buckets=[100, 200, 300], num_tokens=300, exp_bucket=300),
          dict(buckets=[100, 200, 300], num_tokens=301, exp_bucket=301),
          dict(buckets=[100, 200, 300], num_tokens=501, exp_bucket=501),
      ],
      buckets_dtype=[list, tuple],
  )
  def test_calculate_bucket_size(
      self,
      buckets: Sequence[int] | None,
      num_tokens: int,
      exp_bucket: int,
      buckets_dtype,
  ):
    if buckets is not None:
      buckets = buckets_dtype(buckets)
    self.assertEqual(
        pipeline.calculate_bucket_size(num_tokens, buckets),
        exp_bucket,
    )

  @parameterized.product(
      [
          dict(buckets=[], num_tokens=10, error_message='non-empty.'),
          dict(
              buckets=[100, 300, 200],
              num_tokens=10,
              error_message='strictly increasing order.',
          ),
          dict(
              buckets=[100, 200, 200],
              num_tokens=10,
              error_message='strictly increasing order.',
          ),
      ],
      buckets_dtype=[list, tuple],
  )
  def test_calculate_bucket_size_error(
      self,
      buckets: Sequence[int] | None,
      num_tokens: int,
      error_message: str,
      buckets_dtype,
  ):
    if buckets is not None:
      buckets = buckets_dtype(buckets)
    with self.assertRaisesRegex(ValueError, error_message):
      pipeline.calculate_bucket_size(num_tokens, buckets)


class WholePdbPipelineTest(parameterized.TestCase):

  def test_process_item_simple(self):
    bucket = 6
    cfg = pipeline.WholePdbPipeline.Config(buckets=[bucket])
    data_pipeline = pipeline.WholePdbPipeline(config=cfg)

    fold_input = folding_input.Input(
        name='test',
        chains=[
            folding_input.ProteinChain(
                id='A',
                sequence='RRRRR',
                ptms=[],
                unpaired_msa='',
                paired_msa='',
                templates=[],
            ),
        ],
        rng_seeds=[0],
    )

    feats = data_pipeline.process_item(
        fold_input=fold_input,
        random_state=np.random.RandomState(0),
        ccd=chemical_components.Ccd(),
    )

    self.assertEqual(feats['aatype'].shape, (bucket,))
    np.testing.assert_array_equal(feats['aatype'], [1] * 5 + [0] * 1)

    self.assertEqual(feats['residue_index'].shape, (bucket,))
    np.testing.assert_array_equal(feats['residue_index'], [1, 2, 3, 4, 5, 0])

    self.assertEqual(feats['seq_mask'].shape, (bucket,))
    np.testing.assert_array_equal(feats['seq_mask'], [1, 1, 1, 1, 1, 0])

    self.assertEqual(feats['num_alignments'], 2)
    self.assertEqual(feats['seq_length'], 5)

    msa = feats['msa']
    self.assertEqual(msa.shape, (16384, bucket))
    np.testing.assert_array_equal(msa[0], [1, 1, 1, 1, 1, 0])
    np.testing.assert_array_equal(msa[1], [1, 1, 1, 1, 1, 0])
    np.testing.assert_array_equal(msa[2:], 0)

    msa_mask = feats['msa_mask']
    self.assertEqual(msa_mask.shape, (16384, bucket))
    np.testing.assert_array_equal(msa_mask[0], [1, 1, 1, 1, 1, 0])
    np.testing.assert_array_equal(msa_mask[1], [1, 1, 1, 1, 1, 0])
    np.testing.assert_array_equal(msa_mask[2:], 0)

    np.testing.assert_array_equal(feats['is_protein'], [1, 1, 1, 1, 1, 0])
    np.testing.assert_array_equal(feats['is_rna'], [0] * bucket)
    np.testing.assert_array_equal(feats['is_dna'], [0] * bucket)
    np.testing.assert_array_equal(feats['is_ligand'], [0] * bucket)

    self.assertEqual(feats['ref_pos'].shape, (bucket, 24, 3))
    np.testing.assert_array_equal(feats['ref_mask'][:5, 0], [True] * 5)
    np.testing.assert_array_equal(feats['ref_mask'][5:], False)

    self.assertEqual(feats['token_index'].shape, (bucket,))
    np.testing.assert_array_equal(feats['token_index'], [1, 2, 3, 4, 5, 0])

    self.assertEqual(feats['asym_id'].shape, (bucket,))
    np.testing.assert_array_equal(feats['asym_id'], [1, 1, 1, 1, 1, 0])

    self.assertEqual(feats['sym_id'].shape, (bucket,))
    np.testing.assert_array_equal(feats['sym_id'], [1, 1, 1, 1, 1, 0])

    self.assertEqual(feats['entity_id'].shape, (bucket,))
    np.testing.assert_array_equal(feats['entity_id'], [1, 1, 1, 1, 1, 0])

    self.assertEqual(feats['frames_mask'].shape, (bucket,))
    np.testing.assert_array_equal(feats['frames_mask'], [1, 1, 1, 1, 1, 0])

  def test_process_item_mixed(self):
    bucket = 30
    cfg = pipeline.WholePdbPipeline.Config(buckets=[bucket])
    data_pipeline = pipeline.WholePdbPipeline(config=cfg)

    fold_input = folding_input.Input(
        name='test',
        chains=[
            folding_input.ProteinChain(
                id='A',
                sequence='RRRRR',
                ptms=[],
                unpaired_msa='',
                paired_msa='',
                templates=[],
            ),
            folding_input.RnaChain(
                id='B',
                sequence='UUUUU',
                modifications=[],
                unpaired_msa='',
            ),
            folding_input.DnaChain(
                id='C',
                sequence='TTTTT',
                modifications=[],
            ),
            folding_input.Ligand(id='D', ccd_ids=['NAG']),  # NAG has 14 atoms.
        ],
        rng_seeds=[0],
    )

    feats = data_pipeline.process_item(
        fold_input=fold_input,
        random_state=np.random.RandomState(0),
        ccd=chemical_components.Ccd(),
    )

    msa = feats['msa']
    self.assertEqual(msa.shape, (16384, bucket))
    # RRRR (ARG=1), UUUU (U=25), TTTT (DT=29), gap (21), NAG (skipped in MSA)
    np.testing.assert_array_equal(
        msa[0], [1] * 5 + [25] * 5 + [29] * 5 + [21] * 14 + [0] * 1
    )
    np.testing.assert_array_equal(
        msa[1], [1] * 5 + [25] * 5 + [29] * 5 + [21] * 14 + [0] * 1
    )
    np.testing.assert_array_equal(msa[2:], 0)
    msa_mask = feats['msa_mask']
    self.assertEqual(msa_mask.shape, (16384, bucket))
    np.testing.assert_array_equal(msa_mask[0], [1] * 29 + [0] * 1)
    np.testing.assert_array_equal(msa_mask[1], [1] * 29 + [0] * 1)
    np.testing.assert_array_equal(msa_mask[2:], 0)
    self.assertEqual(feats['num_alignments'], 2)

    # ARG (5) + U (5) + DT (5) + NAG (14 atoms) = 29 tokens.
    self.assertEqual(feats['residue_index'].shape, (30,))
    # Protein A (5), RNA B (5), DNA C (5), Ligand D (14 tokens for 14 atoms)
    np.testing.assert_array_equal(
        feats['residue_index'], [1, 2, 3, 4, 5] * 3 + [1] * 14 + [0] * 1
    )
    self.assertEqual(feats['token_index'].shape, (30,))
    np.testing.assert_array_equal(feats['token_index'][:29], np.arange(1, 30))
    self.assertEqual(feats['seq_mask'].shape, (30,))
    np.testing.assert_array_equal(feats['seq_mask'], [True] * 29 + [False] * 1)

    # aatype: ARG=1, U=25, DT=29, Ligand=20 (UNK)
    expected_aatype = [1] * 5 + [25] * 5 + [29] * 5 + [20] * 14 + [0]
    np.testing.assert_array_equal(feats['aatype'], expected_aatype)

    np.testing.assert_array_equal(feats['is_protein'], [1] * 5 + [0] * 25)
    np.testing.assert_array_equal(feats['is_rna'], [0] * 5 + [1] * 5 + [0] * 20)
    np.testing.assert_array_equal(
        feats['is_dna'], [0] * 10 + [1] * 5 + [0] * 15
    )
    np.testing.assert_array_equal(
        feats['is_ligand'], [0] * 15 + [1] * 14 + [0] * 1
    )
    np.testing.assert_array_equal(
        feats['asym_id'], [1] * 5 + [2] * 5 + [3] * 5 + [4] * 14 + [0]
    )
    np.testing.assert_array_equal(
        feats['entity_id'], [1] * 5 + [2] * 5 + [3] * 5 + [4] * 14 + [0]
    )
    np.testing.assert_array_equal(feats['sym_id'], [1] * 29 + [0])
    np.testing.assert_array_equal(feats['frames_mask'], [1] * 29 + [0] * 1)

  def test_process_item_msa(self):
    bucket = 10
    cfg = pipeline.WholePdbPipeline.Config(buckets=[bucket])
    data_pipeline = pipeline.WholePdbPipeline(config=cfg)

    # Chain A: Protein 'RRRRR' (1), Chain B: Protein 'MMMMM' (12)
    fold_input = folding_input.Input(
        name='test',
        chains=[
            folding_input.ProteinChain(
                id='A',
                sequence='RRRRR',
                ptms=[],
                unpaired_msa='>query\nRRRRR\n>unp_a1\nWWWWW\n>unp_a2\nSSSSS',
                paired_msa='>query\nRRRRR\n>tr|12345678|ABC_SPEC1\nGGGGG\n>tr|87654321|XYZ_SPEC2\nDDDDD',
                templates=[],
            ),
            folding_input.ProteinChain(
                id='B',
                sequence='MMMMM',
                ptms=[],
                unpaired_msa='',
                paired_msa='>query\nMMMMM\n>tr|12345678|DEF_SPEC1\nLLLLL\n>tr|87654321|UVW_SPEC2\nFFFFF',
                templates=[],
            ),
        ],
        rng_seeds=[0],
    )

    random_state = np.random.RandomState(0)

    feats = data_pipeline.process_item(
        fold_input=fold_input,
        random_state=random_state,
        ccd=chemical_components.Ccd(),
    )

    self.assertEqual(feats['num_alignments'], 5)

    msa = feats['msa']
    self.assertEqual(msa.shape, (16384, bucket))
    np.testing.assert_array_equal(msa[0], [1] * 5 + [12] * 5)  #  RRRRRMMMMM
    np.testing.assert_array_equal(msa[1], [7] * 5 + [10] * 5)  #  GGGGGLLLLL
    np.testing.assert_array_equal(msa[2], [3] * 5 + [13] * 5)  #  DDDDDFFFFF
    np.testing.assert_array_equal(msa[3], [17] * 5 + [21] * 5)  # WWWWW-----
    np.testing.assert_array_equal(msa[4], [15] * 5 + [21] * 5)  # SSSSS-----
    np.testing.assert_array_equal(msa[5:], 0)

    msa_mask = feats['msa_mask']
    self.assertEqual(msa_mask.shape, (16384, bucket))
    np.testing.assert_array_equal(msa_mask[:5], True)
    np.testing.assert_array_equal(msa_mask[5:], False)

    np.testing.assert_array_equal(feats['residue_index'], [1, 2, 3, 4, 5] * 2)
    np.testing.assert_array_equal(feats['token_index'], np.arange(1, 11))
    np.testing.assert_array_equal(feats['seq_mask'], [True] * bucket)
    np.testing.assert_array_equal(feats['aatype'], [1] * 5 + [12] * 5)
    np.testing.assert_array_equal(feats['is_protein'], [True] * bucket)
    np.testing.assert_array_equal(feats['is_rna'], [False] * bucket)
    np.testing.assert_array_equal(feats['is_dna'], [False] * bucket)
    np.testing.assert_array_equal(feats['is_ligand'], [False] * bucket)
    np.testing.assert_array_equal(feats['asym_id'], [1] * 5 + [2] * 5)
    np.testing.assert_array_equal(feats['entity_id'], [1] * 5 + [2] * 5)
    np.testing.assert_array_equal(feats['sym_id'], [1] * bucket)
    np.testing.assert_array_equal(feats['frames_mask'], [True] * bucket)


if __name__ == '__main__':
  absltest.main()
