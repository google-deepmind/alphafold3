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

from unittest import mock

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.common import resources
from alphafold3.common.testing import data
from alphafold3.constants import mmcif_names
from alphafold3.data import msa
from alphafold3.data import msa_config
import numpy as np

import shutil

_JACKHMMER_BINARY_PATH = shutil.which('jackhmmer')
_NHMMER_BINARY_PATH = shutil.which('nhmmer')
_HMMALIGN_BINARY_PATH = shutil.which('hmmalign')
_HMMBUILD_BINARY_PATH = shutil.which('hmmbuild')

_SEQ_T1014 = (
    'SLAPVDIEGLLRQVAELMSPRAHEKGIEIAWAVSSPLPTILADEGRLRQILLNFAGNAVKFTEAGGVLLTASAI'
    'DGGRVRFSVADTGPGVAPDARARIFEAFVQTDVTHATQLGGAGLGLAIVSRLSAAMGGAVGVGGELGQGAEFWF'
    'EAPFATAAAPLRAAPLEGRNVAIASPNAIVRAATARQIEAAGGRAYAAVDIASALAGAPADAVLLIDAALSGPR'
    'GALKPPAGRRSVVLLTPEQRDRIDRLKAAGFSGYLIKPLRAASLVAQVLQAVTA'
)

_SEQ_5FK3 = (
    'GGCUUAUCAAGAGAGGGCGAGCGACUGGCGCGAAGACCCCCGGCAACCAGAAAUGGUGCCAAUUCCUGCAGCGG'
    'AAACGUUGAAAGAUGAGCCG'
)

_SEQ_T0760 = (
    'MYMKKLKFITLTLAIIITLPMLQSCLDDNDHQSRSLVISTINQISEDSKEFYFTLDNGKTMFPSNSQAWGGEKF'
    'ENGQRAFVIFNELEQPVNGYDYNIQVRDITKVLTKEIVTMDDEENTEEKIGDDKINATYMWISKDKKYLTIEFQ'
    'YYSTHSEDKKHFLNLVINNKEADSAAENEDNTDDEYINLEFRHNSERDSPDHLGEGYVSFKLDKIEEQIEGKKG'
    'LNIRVRTLYDGIKNYKVQFP'
)

_RNA_TEST_SEQ = 'CCAAAAATACGAAAAAGTAGCCAGGCATAATGGCACACATGTGTGGTCCCGCTGCTTGGG'


def _msa_from_seqs(seqs):
  return msa.Msa(
      query_sequence=seqs[0],
      chain_poly_type=mmcif_names.PROTEIN_CHAIN,
      sequences=seqs,
      descriptions=[''] * len(seqs),
  )


class MsaTest(parameterized.TestCase):

  def setUp(self):
    super().setUp()
    self.jackhmmer_test_data = data.Data(
        resources.ROOT / 'data/tools/testdata/jackhmmer'
    )
    self.nhmmer_test_data = data.Data(
        resources.ROOT / 'data/tools/testdata/nhmmer'
    )

  def _msa_test_query(
      self, chain_poly_type: str, crop_size: int | None
  ) -> msa_config.RunConfig:
    match chain_poly_type:
      case mmcif_names.PROTEIN_CHAIN:
        msa_tool_config = msa_config.JackhmmerConfig(
            database_config=msa_config.DatabaseConfig(
                name='test_db',
                path=self.jackhmmer_test_data.path('test_db.fasta'),
            ),
            binary_path=_JACKHMMER_BINARY_PATH,
            n_cpu=1,
            n_iter=3,
            e_value=1e-3,
            z_value=None,
            dom_z_value=None,
            max_sequences=5000,
        )
      case mmcif_names.RNA_CHAIN:
        msa_tool_config = msa_config.NhmmerConfig(
            database_config=msa_config.DatabaseConfig(
                name='rfam',
                path=self.nhmmer_test_data.path('dna_target.fasta'),
            ),
            binary_path=_NHMMER_BINARY_PATH,
            hmmalign_binary_path=_HMMALIGN_BINARY_PATH,
            hmmbuild_binary_path=_HMMBUILD_BINARY_PATH,
            n_cpu=1,
            e_value=1e1,
            z_value=None,
            alphabet=None,
            max_sequences=32,
        )
      case _:
        raise ValueError(f'Unsupported {chain_poly_type=}.')

    return msa_config.RunConfig(
        chain_poly_type=chain_poly_type,
        config=msa_tool_config,
        crop_size=crop_size,
    )

  @parameterized.named_parameters(
      ('max_depth', list('ABC'), list('123'), list('AB'), list('12'), 2, 2),
      ('max_depth_zero', list('AB'), list('12'), list('AB'), list('12'), 0, 2),
      ('max_depth_equal', list('AB'), list('12'), list('AB'), list('12'), 2, 2),
      ('max_depth_more', list('AB'), list('12'), list('AB'), list('12'), 3, 2),
      ('empty', [], [], ['A'], ['Original query'], None, 1),
  )
  def test_single_msa(
      self,
      seqs,
      descs,
      exp_seqs,
      exp_descs,
      max_depth,
      exp_depth,
  ):
    a3m = '\n'.join([f'>{desc}\n{seq}' for seq, desc in zip(seqs, descs)])
    single_msa = msa.Msa.from_a3m(
        query_sequence='A',
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
        a3m=a3m,
        max_depth=max_depth,
    )
    self.assertEqual(single_msa.depth, exp_depth)
    self.assertEqual(single_msa.sequences, exp_seqs)
    self.assertEqual(single_msa.descriptions, exp_descs)
    self.assertEqual(
        str(single_msa),
        f'Msa({exp_depth} sequences, {mmcif_names.PROTEIN_CHAIN})',
    )

  def test_single_msa_from_a3m(self):
    single_msa = msa.Msa.from_a3m(
        query_sequence='ABJOUXZ',
        a3m='>seq1\nADXXCXE\n>seq2\nDEF',
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
    )
    self.assertEqual(single_msa.depth, 2)
    self.assertEqual(single_msa.sequences, ['ADXXCXE', 'DEF'])
    self.assertEqual(single_msa.descriptions, ['seq1', 'seq2'])

  @parameterized.named_parameters(
      ('bad_lengths', 'A', ['A', 'B'], ['1', '2', '3']),
      ('bad_first_seq', 'AGA', ['A', 'B'], ['1', '2']),
  )
  def test_single_msa_error(self, query_sequence, seqs, descs):
    with self.assertRaises(ValueError):
      msa.Msa(
          query_sequence=query_sequence,
          chain_poly_type=mmcif_names.PROTEIN_CHAIN,
          sequences=seqs,
          descriptions=descs,
      )

  @parameterized.named_parameters(
      ('deduplicated', True, 4),
      ('not_deduplicated', False, 6),
  )
  def test_from_multiple_a3ms(self, deduplicate, exp_depth):
    single_a3m_1 = _msa_from_seqs(['A', 'B', 'Caa']).to_a3m()
    single_a3m_2 = _msa_from_seqs(['A', 'Coo', 'D']).to_a3m()

    merged_msa = msa.Msa.from_multiple_a3ms(
        a3ms=[single_a3m_1, single_a3m_2],
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
        deduplicate=deduplicate,
    )
    self.assertEqual(merged_msa.depth, exp_depth)

  @parameterized.named_parameters(
      ('deduplicated', True, 4),
      ('not_deduplicated', False, 6),
  )
  def test_from_multiple_msas(self, deduplicate, exp_depth):
    single_msa_1 = _msa_from_seqs(['A', 'B', 'Caa'])
    single_msa_2 = _msa_from_seqs(['A', 'Coo', 'D'])

    merged_msa = msa.Msa.from_multiple_msas(
        msas=[single_msa_1, single_msa_2], deduplicate=deduplicate
    )
    self.assertEqual(merged_msa.depth, exp_depth)

  def test_empty_msa(self):
    empty_msa = msa.Msa.from_empty(
        query_sequence='DEEPMIND', chain_poly_type=mmcif_names.PROTEIN_CHAIN
    )
    self.assertEqual(empty_msa.depth, 1)
    self.assertEqual(empty_msa.to_a3m(), '>Original query\nDEEPMIND\n')

  @parameterized.named_parameters(
      ('empty', [], []),
      ('queries', [['A', 'B'], ['C', 'D']], [mmcif_names.PROTEIN_CHAIN] * 2),
      ('type', [['A']] * 2, [mmcif_names.RNA_CHAIN, mmcif_names.PROTEIN_CHAIN]),
  )
  def test_merged_msa_error(self, msa_seqs, poly_types):
    msas = [_msa_from_seqs(seqs) for seqs in msa_seqs] if msa_seqs else []

    for single_msa, poly_type in zip(msas, poly_types, strict=True):
      single_msa.chain_poly_type = poly_type
      single_msa.tool_config = None

    with self.assertRaises(ValueError):
      msa.Msa.from_multiple_msas(msas=msas, deduplicate=True)

  def test_get_alignment_features(self):
    msa_sequences = [
        'ABCDEF',
        'e-BCDEFss',
        'Ai-CiDiEF',
        'ABCD-o-o',
    ]
    msa_descriptions = [
        'tr|AAAAA1|AAAAA1_SSS1 OS=Homo sapiens OX=9606',
        'sp|AAAAA2|PPP2_SSS2 OS=E. hepatopenaei CLUSTER:cluster_3:CLUSTER',
        'seq_3 OX=646526 CLUSTER:cluster_2:CLUSTER',
        'seq_4 CLUSTER:cluster_4:CLUSTER',
    ]
    merged_msa = msa.Msa(
        query_sequence='ABCDEF',
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
        sequences=msa_sequences,
        descriptions=msa_descriptions,
    )

    expected_msa = np.array(
        [
            [0, 3, 4, 3, 6, 13],
            [21, 3, 4, 3, 6, 13],
            [0, 21, 4, 3, 6, 13],
            [0, 3, 4, 3, 21, 21],
        ],
        dtype=np.int32,
    )
    expected_deletion_matrix = np.array(
        [
            [0, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0],
            [0, 1, 0, 1, 1, 0],
            [0, 0, 0, 0, 0, 1],
        ],
        dtype=np.int32,
    )
    expected_species_identifiers = np.array(
        [b'SSS1', b'SSS2', b'', b''], dtype=object
    )

    features = merged_msa.featurize()
    np.testing.assert_array_equal(features['num_alignments'], 4)
    np.testing.assert_array_equal(features['msa'], expected_msa, strict=True)
    deletion_matrix = features['deletion_matrix']
    np.testing.assert_array_equal(
        deletion_matrix, expected_deletion_matrix, strict=True
    )
    species_identifiers = features['msa_species_identifiers']
    np.testing.assert_array_equal(
        species_identifiers, expected_species_identifiers, strict=True
    )

  def test_bad_amino_acids(self):
    msa_sequences = [
        'ABCDEF',
        '> Actually a comment',
    ]
    err = (
        'Error extracting MSA or deletion features: Unknown residues in MSA: '
        '(1, 0):>, (1, 1): , (1, 10): , (1, 12): . target_sequence: ABCDEF'
    )
    merged_msa = msa.Msa.from_multiple_msas(
        msas=[_msa_from_seqs(msa_sequences)], deduplicate=False
    )
    with self.assertRaisesWithLiteralMatch(msa.Error, err):
      merged_msa.featurize()

  def test_bad_alignment(self):
    msa_sequences = [
        'ABCDEF',
        'GH-JKLb',
        'MNOPQ-R',
    ]
    err = (
        'Error extracting MSA or deletion features: Invalid shape all strings '
        'must have the same number of non-lowercase characters; First string '
        "has 6 non-lowercase characters but 'MNOPQ-R' has 7. target_sequence: "
        'ABCDEF'
    )
    merged_msa = msa.Msa.from_multiple_msas(
        msas=[_msa_from_seqs(msa_sequences)], deduplicate=False
    )
    with self.assertRaisesWithLiteralMatch(msa.Error, err):
      merged_msa.featurize()

  @parameterized.parameters(
      (_SEQ_T0760, mmcif_names.PROTEIN_CHAIN, None, 2),
      (_RNA_TEST_SEQ, mmcif_names.RNA_CHAIN, None, 9),
      # MSA depth is less than crop size due to deduplication occurring after
      # cropping.
      (_SEQ_T0760, mmcif_names.PROTEIN_CHAIN, 2, 1),
      (_RNA_TEST_SEQ, mmcif_names.RNA_CHAIN, 7, 6),
  )
  def test_get_msa(self, sequence, chain_poly_type, crop_size, exp_msa_depth):
    msa_query = self._msa_test_query(
        chain_poly_type=chain_poly_type, crop_size=crop_size
    )
    merged_msa = msa.get_msa(
        target_sequence=sequence,
        run_config=msa_query,
        chain_poly_type=chain_poly_type,
        deduplicate=True,
    )
    self.assertEqual(merged_msa.depth, exp_msa_depth)
    self.assertEqual(sequence, merged_msa.sequences[0])

  @mock.patch.object(msa.jackhmmer, 'Jackhmmer', autospec=True)
  def test_get_msa_tool_jackhmmer(self, mock_jackhmmer_cls):
    """Tests that all JackhmmerConfig fields are forwarded to Jackhmmer."""
    n_cpu = 2
    n_iter = 1
    e_value = 1e-4
    z_value = 200
    dom_z_value = 400
    max_sequences = 5432
    max_parallel_shards = 32
    config = msa_config.JackhmmerConfig(
        binary_path='jackhmmer',
        database_config=msa_config.DatabaseConfig(name='tst', path='db.fasta'),
        n_cpu=n_cpu,
        n_iter=n_iter,
        e_value=e_value,
        z_value=z_value,
        dom_z_value=dom_z_value,
        max_sequences=max_sequences,
        max_parallel_shards=max_parallel_shards,
    )

    msa.get_msa_tool(config)

    mock_jackhmmer_cls.assert_called_once_with(
        binary_path='jackhmmer',
        database_path='db.fasta',
        n_cpu=n_cpu,
        n_iter=n_iter,
        e_value=e_value,
        z_value=z_value,
        dom_z_value=dom_z_value,
        max_sequences=max_sequences,
        max_threads=max_parallel_shards,
    )

  @mock.patch.object(msa.nhmmer, 'Nhmmer', autospec=True)
  def test_get_msa_tool_nhmmer(self, mock_nhmmer_cls):
    """Tests that all NhmmerConfig fields are forwarded to Nhmmer."""
    n_cpu = 2
    e_value = 1e-7
    z_value = 500.0
    max_sequences = 4321
    alphabet = 'dna'
    max_parallel_shards = 16
    config = msa_config.NhmmerConfig(
        binary_path='nhmmer',
        hmmalign_binary_path='hmmalign',
        hmmbuild_binary_path='hmmbuild',
        database_config=msa_config.DatabaseConfig(name='tst', path='db.fasta'),
        n_cpu=n_cpu,
        e_value=e_value,
        z_value=z_value,
        max_sequences=max_sequences,
        alphabet=alphabet,
        max_parallel_shards=max_parallel_shards,
    )

    msa.get_msa_tool(config)

    mock_nhmmer_cls.assert_called_once_with(
        binary_path='nhmmer',
        hmmalign_binary_path='hmmalign',
        hmmbuild_binary_path='hmmbuild',
        database_path='db.fasta',
        n_cpu=n_cpu,
        e_value=e_value,
        z_value=z_value,
        max_sequences=max_sequences,
        alphabet=alphabet,
        max_threads=max_parallel_shards,
    )


if __name__ == '__main__':
  absltest.main()
