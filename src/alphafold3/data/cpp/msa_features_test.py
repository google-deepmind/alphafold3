# Copyright 2026 DeepMind Technologies Limited
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
import sys

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.constants import mmcif_names
from alphafold3.constants import residue_names
from alphafold3.cpp import msa_features
import numpy as np


_INVALID_MSA_SHORT_SEQUENCES = [
    'ABCDEF',
    'e-BCDEFss',
    'Ai-CiDiEF',
    'ABCD-oo',
]
_SHORT_ERR_MSG = (
    'Invalid shape all strings must have the same number of non-lowercase '
    "characters; First string has 6 non-lowercase characters but 'ABCD-oo' "
    'has 5. target_sequence: ABCDEF'
)

_INVALID_MSA_LONG_SEQUENCES = [
    'ABCDEF',
    'GH-JKLb',
    'MNOPQ-R',
]
_LONG_ERR_MSG = (
    'Invalid shape all strings must have the same number of non-lowercase'
    " characters; First string has 6 non-lowercase characters but 'MNOPQ-R' has"
    ' 7. target_sequence: ABCDEF'
)


class MsaOpsTest(parameterized.TestCase):

  @parameterized.named_parameters(
      ('protein', mmcif_names.PROTEIN_CHAIN, 0),
      ('rna', mmcif_names.RNA_CHAIN, 22),
      ('dna', mmcif_names.DNA_CHAIN, 26),
  )
  def test_extract_msa_features_basic(self, chain_poly_type, a_id):
    # Test basic features, the deletion matrix and correct reference counts.
    msa_sequences = [
        'AAAAAA',
        'deepmind-AAAAAss',
        'Ai-AiAiAA',
        'AAAA-o-o',
    ]
    expected_msa = np.array([
        [a_id, a_id, a_id, a_id, a_id, a_id],
        [21, a_id, a_id, a_id, a_id, a_id],
        [a_id, 21, a_id, a_id, a_id, a_id],
        [a_id, a_id, a_id, a_id, 21, 21],
    ])
    expected_deletion_matrix = np.array([
        [0, 0, 0, 0, 0, 0],
        [8, 0, 0, 0, 0, 0],
        [0, 1, 0, 1, 1, 0],
        [0, 0, 0, 0, 0, 1],
    ])

    msa, deletion_matrix = msa_features.extract_msa_features(
        msa_sequences, chain_poly_type
    )

    self.assertEqual(msa.dtype, np.int32)
    self.assertEqual(deletion_matrix.dtype, np.int32)

    self.assertEqual(sys.getrefcount(msa), sys.getrefcount(expected_msa))
    self.assertEqual(
        sys.getrefcount(deletion_matrix),
        sys.getrefcount(expected_deletion_matrix),
    )

    np.testing.assert_array_equal(deletion_matrix, expected_deletion_matrix)
    np.testing.assert_array_equal(msa, expected_msa)

  @parameterized.named_parameters(
      (
          'protein_short',
          mmcif_names.PROTEIN_CHAIN,
          _INVALID_MSA_SHORT_SEQUENCES,
          _SHORT_ERR_MSG,
      ),
      (
          'rna_short',
          mmcif_names.RNA_CHAIN,
          _INVALID_MSA_SHORT_SEQUENCES,
          _SHORT_ERR_MSG,
      ),
      (
          'dna_short',
          mmcif_names.DNA_CHAIN,
          _INVALID_MSA_SHORT_SEQUENCES,
          _SHORT_ERR_MSG,
      ),
      (
          'protein_long',
          mmcif_names.PROTEIN_CHAIN,
          _INVALID_MSA_LONG_SEQUENCES,
          _LONG_ERR_MSG,
      ),
      (
          'rna_long',
          mmcif_names.RNA_CHAIN,
          _INVALID_MSA_LONG_SEQUENCES,
          _LONG_ERR_MSG,
      ),
      (
          'dna_long',
          mmcif_names.DNA_CHAIN,
          _INVALID_MSA_LONG_SEQUENCES,
          _LONG_ERR_MSG,
      ),
  )
  def test_invalid_shape_throws(
      self, chain_poly_type: str, msa_sequences: Sequence[str], err_msg: str
  ):
    with self.assertRaisesWithLiteralMatch(ValueError, err_msg):
      msa_features.extract_msa_features(msa_sequences, chain_poly_type)

  @parameterized.named_parameters(
      ('protein', mmcif_names.PROTEIN_CHAIN),
      ('rna', mmcif_names.RNA_CHAIN),
      ('dna', mmcif_names.DNA_CHAIN),
  )
  def test_invalid_character_throws(self, chain_poly_type):
    invalid_msa_sequences = [
        '>ABCDEF',
        '>e-BCDEFss',
        '>Ai-CiDiEF',
        '>ABCD-o-o',
    ]
    with self.assertRaises(ValueError):
      msa_features.extract_msa_features(invalid_msa_sequences, chain_poly_type)

  def test_invalid_chain_poly_type(self):
    with self.assertRaises(ValueError):
      msa_features.extract_msa_features(['ABCDEF'], 'INVALID_CHAIN')

  @parameterized.named_parameters(
      (
          'protein',
          mmcif_names.PROTEIN_CHAIN,
          residue_names.PROTEIN_TYPES_WITH_UNKNOWN,
      ),
      (
          'rna',
          mmcif_names.RNA_CHAIN,
          residue_names.RNA_TYPES + (residue_names.UNK_NUCLEIC_ONE_LETTER,),
      ),
      (
          # Note UNK DNA and UNK RNA are mapped to the same residue: "N".
          'dna',
          mmcif_names.DNA_CHAIN,
          residue_names.DNA_TYPES + (residue_names.UNK_NUCLEIC_ONE_LETTER,),
      ),
  )
  def test_residue_encoding(self, chain_poly_type, residues):
    sequence = ''.join(
        [residue_names.CCD_NAME_TO_ONE_LETTER.get(c, '-') for c in residues]
    )
    msa, _ = msa_features.extract_msa_features([sequence], chain_poly_type)
    encoded = msa[0].tolist()
    exp_encoded = [
        residue_names.POLYMER_TYPES_ORDER_WITH_UNKNOWN_AND_GAP[r]
        for r in residues
    ]
    self.assertSequenceEqual(encoded, exp_encoded)

  @parameterized.named_parameters(
      ('protein', mmcif_names.PROTEIN_CHAIN),
      ('rna', mmcif_names.RNA_CHAIN),
      ('dna', mmcif_names.DNA_CHAIN),
  )
  def test_extract_msa_features_gap_encoding(self, chain_poly_type):
    msa, _ = msa_features.extract_msa_features(
        [residue_names.GAP], chain_poly_type
    )
    exp_encoded = residue_names.POLYMER_TYPES_ORDER_WITH_UNKNOWN_AND_GAP[
        residue_names.GAP
    ]
    self.assertEqual(msa[0][0], exp_encoded)

  def test_extract_msa_features_non_standard_protein_residues(self):
    msa, _ = msa_features.extract_msa_features(
        ['BJOUXZ'], mmcif_names.PROTEIN_CHAIN
    )
    exp_msa, _ = msa_features.extract_msa_features(
        ['DXXCXE'], mmcif_names.PROTEIN_CHAIN
    )
    self.assertSequenceEqual(msa[0].tolist(), exp_msa[0].tolist())

  def test_empty(self):
    msa_sequences = []
    msa, deletion_matrix = msa_features.extract_msa_features(
        msa_sequences, mmcif_names.PROTEIN_CHAIN
    )
    self.assertTupleEqual(msa.shape, (0, 0))
    self.assertTupleEqual(deletion_matrix.shape, (0, 0))

  @parameterized.named_parameters(
      ('empty', [''], [b'']),
      ('non_match', ['foo bar'], [b'']),
      ('bad_start', ['foo sp|P0C2L1|A3X1_LOXLA'], [b'']),
      ('swissprot', ['sp|P0C2L1|A3X1_LOXLA'], [b'LOXLA']),
      ('trembl', ['tr|A0A146SKV9|A0A146SKV9_FUNHE'], [b'FUNHE']),
      ('match_non_match', ['sp|P0C2L1|A3X1_LOXLA', 'foo'], [b'LOXLA', b'']),
      ('double_match', ['sp|P0C2L1|A3X1_A', 'sp|P0C2L1|A3X1_B'], [b'A', b'B']),
      ('sp_with_text', ['sp|P0C2L1|A3X1_LOXLA foo bar'], [b'LOXLA']),
      ('sp_with_space_after', ['sp|P0C2L1|A3X1_LOXLA '], [b'LOXLA']),
      ('sp_with_space_before', [' sp|P0C2L1|A3X1_LOXLA'], [b'LOXLA']),
      ('isoform', ['sp|P0C2L1_2|A3X1_LOXLA'], [b'LOXLA']),
      ('alt_isoform_excluded', ['sp|O88974-2|SETB1_MOUSE'], [b'']),
      (
          'subseq_range',
          ['sp|Q32KD2|SETB1_DROME/444-700 [subseq from] Histone'],
          [b'DROME'],
      ),
      (
          'non_uniprot_contig',
          ['MDSW01.2.fsa_nt_gb|MDSW01133750.1|_1 # 131'],
          [b''],
      ),
  )
  def test_extract_species_ids(self, descriptions, exp_species_ids):
    species_ids = msa_features.extract_species_ids(descriptions)
    self.assertEqual(species_ids, exp_species_ids)


if __name__ == '__main__':
  absltest.main()
