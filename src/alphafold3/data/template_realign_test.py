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

from collections.abc import Mapping, Sequence
import csv

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.common import resources
from alphafold3.data import template_realign


def _align_hit(
    hit_length: int, alignment: Mapping[int, int], struc_sequence: str
) -> str:
  aligned_hit = ['-'] * hit_length
  for query_index in range(hit_length):
    hit_index = alignment.get(query_index)
    if hit_index is not None:
      aligned_hit[query_index] = struc_sequence[hit_index]
  return ''.join(aligned_hit)


def _load_unequal_pdb_seqs() -> Sequence[Mapping[str, str | int]]:
  unequal_sequences = []
  with resources.open_resource(
      resources.ROOT
      / 'data/testdata/template_realign/mismatch_seqres_vs_mmcif_seqs.csv',
      mode='rt',
  ) as mismatch_seqs_csv:
    for row in csv.DictReader(mismatch_seqs_csv):
      max_mismatch_columns = 1 if row['target_name'].startswith('4ryw') else 0
      unequal_sequences.append({
          'testcase_name': row['target_name'],
          'seqres_seq': row['seqres_seq'],
          'struc_seq': row['struc_seq'],
          'max_mismatch_columns': max_mismatch_columns,
      })
  return unequal_sequences


_UNEQUAL_SEQUENCES: Sequence[Mapping[str, str | int]] = _load_unequal_pdb_seqs()


class TemplateRealignTest(parameterized.TestCase):

  @parameterized.parameters(
      ('ABC', 0, 'A', 'A--'),
      ('ABC', 0, 'C', '--C'),
      ('ABC', 0, 'X', 'X--'),
      ('ABCD', 0, 'AD', 'A--D'),
      ('ABCD', 0, 'BD', '-B-D'),
      ('ABCD', 0, 'CD', '--CD'),
      ('ABCD', 0, 'XD', 'X--D'),
      ('ABCD', 0, 'ABX', 'ABX-'),
      ('ABBD', 0, 'ABD', 'AB-D'),
      ('ABBBBC', 0, 'ABBC', 'ABB--C'),
      ('ABBBD', 0, 'ABD', 'AB--D'),
      ('ABCCD', 0, 'ABD', 'AB--D'),
      ('ABCDEFG', 0, 'ABFG', 'AB---FG'),
      ('ABCDEFG', 0, 'ABCFG', 'ABC--FG'),
      ('ABCDEFG', 0, 'ABCDFG', 'ABCD-FG'),
      ('ABCDEFGH', 0, 'XDEG', 'X--DE-G-'),
      ('ABCDEFGH', 0, 'XDEXGH', 'X--DEXGH'),
      ('ACDEIKFNNGH', 0, 'ACFWNGH', 'AC----FWNGH'),
      ('GSMVSLPTPWPTLDDXLMCF', 0, 'XVSLFTDWPTLXLMCF', 'X--VSLFTDWPTL--XLMCF'),
      ('ABCDEFGH', 0, 'ABCDEFGH', 'ABCDEFGH'),
      ('ABCDEFGH', 1, 'BCDEFGH', '-BCDEFGH'),
      ('ABCDEFGH', 2, 'CDEFGH', '--CDEFGH'),
      ('ABCDEFGH', 3, 'CDEFGH', '---DEFGH'),
      ('AABCDEFGH', 1, 'XBCDEFGH', '-XBCDEFGH'),
      ('AABCDEFGH', 2, 'XXCDEFGH', '--XCDEFGH'),
      ('AABCDEFGH', 1, 'BCDEFGH', '--BCDEFGH'),
      ('ABCDEFGH', 0, 'ABCDEGH', 'ABCDE-GH'),
      ('ABCDEFGH', 1, 'BCDEGH', '-BCDE-GH'),
      ('ABCDEFGH', 1, 'XBCDEGH', '-BCDE-GH'),
      ('ABCDEFGH', 2, 'CDEGH', '--CDE-GH'),
      ('ABCDEFGH', 2, 'XCDEGH', '--CDE-GH'),
      ('ABCDEFGH', 2, 'XXCDEGH', '--CDE-GH'),
      ('DEEPMATCHMIND', 4, 'CATCHMIND', '----CATCHMIND'),
      ('DEEPMATCHMIND', 4, 'CATCHMIN', '----CATCHMIN-'),
      ('DEEPMATCHMIND', 4, 'CATCHMID', '----CATCHMI-D'),
      ('DEEPMATCHMIND', 4, 'CATCHID', '----CATCH-I-D'),
  )
  def test_realign_hit_to_structure(
      self,
      hit_sequence: str,
      hit_start_index: int,
      struc_sequence: str,
      expected_aligned_hit: str,
  ):
    # We assume that query and the hmm hit are equal so there is identity map.
    query_to_hit_mapping = {
        i: i - hit_start_index
        for i in range(hit_start_index, len(hit_sequence))
    }
    new_alignment = template_realign.realign_hit_to_structure(
        hit_sequence=hit_sequence[hit_start_index:],
        hit_start_index=hit_start_index,
        hit_end_index=len(hit_sequence),
        full_length=len(hit_sequence),
        structure_sequence=struc_sequence,
        query_to_hit_mapping=query_to_hit_mapping,
    )

    aligned_hit = _align_hit(len(hit_sequence), new_alignment, struc_sequence)
    self.assertSequenceEqual(aligned_hit, expected_aligned_hit)
    self.assertLen(new_alignment, len(expected_aligned_hit.replace('-', '')))

  @parameterized.parameters(
      dict(
          query_sequence='ABCDEFGHIJKL',
          hit_sequence='---DEFGHIJK-',
          struc_sequence='DEFGHKL',
          expected_aligned_hit='---DEFGH--K-',
          start_index=3,
          end_index=11,
          full_length=12,
          query_to_hit_mapping={
              3: 0, 4: 1, 5: 2, 6: 3, 7: 4, 8: 5, 9: 6, 10: 7
          },
      ),
      dict(
          query_sequence='ABCDEFGHIJKL',
          hit_sequence='---DEFGHIJK-',
          struc_sequence='XDEFGHKL',
          expected_aligned_hit='---DEFGH--K-',
          start_index=3,
          end_index=11,
          full_length=12,
          query_to_hit_mapping={
              3: 0, 4: 1, 5: 2, 6: 3, 7: 4, 8: 5, 9: 6, 10: 7
          },
      ),
      dict(
          query_sequence='ABCDEFGHIJKL',
          hit_sequence='--CDEFGHIJK-',
          struc_sequence='XDEFGHKL',
          expected_aligned_hit='--XDEFGH--K-',
          start_index=2,
          end_index=11,
          full_length=12,
          query_to_hit_mapping={
              2: 0, 3: 1, 4: 2, 5: 3, 6: 4, 7: 5, 8: 6, 9: 7, 10: 8
          },
      ),
      dict(
          query_sequence='ABCDEFGHIJKL',
          hit_sequence='---DEFGHIJKL',
          struc_sequence='XDEFGHKL',
          expected_aligned_hit='---DEFGH--KL',
          start_index=3,
          end_index=12,
          full_length=12,
          query_to_hit_mapping={
              3: 0, 4: 1, 5: 2, 6: 3, 7: 4, 8: 5, 9: 6, 10: 7, 11: 8
          },
      ),
      dict(
          query_sequence='ABCDEFGHIJKL',
          hit_sequence='---DEFGHIJ--',
          struc_sequence='XDEFGHKL',
          expected_aligned_hit='---DEFGHKL--',
          start_index=3,
          end_index=10,
          full_length=12,
          query_to_hit_mapping={
              3: 0, 4: 1, 5: 2, 6: 3, 7: 4, 8: 5, 9: 6
          },
      ),
  )  # pyformat: disable
  def test_realign_hit_to_structure_complex(
      self,
      query_sequence: str,
      hit_sequence: str,
      struc_sequence: str,
      expected_aligned_hit: str,
      start_index: int,
      end_index: int,
      full_length: int,
      query_to_hit_mapping: Mapping[int, int],
  ):
    # Test the case where there are gaps before the aligned region.
    new_alignment = template_realign.realign_hit_to_structure(
        hit_sequence=hit_sequence.upper().replace('-', ''),
        hit_start_index=start_index,
        hit_end_index=end_index,
        full_length=full_length,
        structure_sequence=struc_sequence,
        query_to_hit_mapping=query_to_hit_mapping,
    )

    aligned_hit = _align_hit(len(query_sequence), new_alignment, struc_sequence)
    self.assertSequenceEqual(''.join(aligned_hit), expected_aligned_hit)
    self.assertLen(new_alignment, len(expected_aligned_hit.replace('-', '')))

  @parameterized.named_parameters(_UNEQUAL_SEQUENCES)
  def test_realign_all_unequal_pdb_sequences(
      self,
      seqres_seq: str,
      struc_seq: str,
      max_mismatch_columns: int,
  ):
    alignment = template_realign.realign_hit_to_structure(
        hit_sequence=seqres_seq,
        hit_start_index=0,
        hit_end_index=len(seqres_seq),
        full_length=len(seqres_seq),
        structure_sequence=struc_seq,
        query_to_hit_mapping={i: i for i in range(len(seqres_seq))},
    )

    aligned_seq = _align_hit(len(seqres_seq), alignment, struc_seq)

    num_aligned_columns_before = sum(
        r1 == r2 or r2 in {'X', '-'}
        for r1, r2 in zip(seqres_seq, struc_seq, strict=False)
    )
    num_aligned_columns_after = sum(
        r1 == r2 or r2 in {'X', '-'}
        for r1, r2 in zip(seqres_seq, aligned_seq, strict=True)
    )
    # The re-alignment improves number of aligned columns.
    self.assertLess(num_aligned_columns_before, num_aligned_columns_after)
    # The number of unaligned non-gap + non-unknown columns is at most X.
    self.assertLessEqual(
        len(aligned_seq) - num_aligned_columns_after, max_mismatch_columns
    )
    # The number of inserted gaps equals the difference in lengths.
    self.assertEqual(aligned_seq.count('-'), len(seqres_seq) - len(struc_seq))
    self.assertLen(alignment, len(struc_seq))


if __name__ == '__main__':
  absltest.main()
