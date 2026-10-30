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

import io

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.common import resources
from alphafold3.common.testing import data
from alphafold3.data import parsers


class ParsersTest(parameterized.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.test_data = data.Data(resources.ROOT / 'data/testdata/parsers')

  def test_lazy_parse_fasta_string(self):
    fasta_str = self.test_data.load('test_db.fasta')
    fasta_iterator = parsers.lazy_parse_fasta_string(fasta_str)
    sequences = [s[0] for s in fasta_iterator]
    self.assertLen(sequences, 34)
    self.assertEqual(
        sequences[26],
        (
            'MEQVPIKEMRLSDLRPNNKSIDTDLGGTKLVVIGKPGSGKSTLIKALLDSKRHIIPCAVVISGSEE'
            'ANGFYKGVVPDLFIYHQFSPSIIDRIHRRQVKAKAEMGSKKSWLLVVIDDCMDNAKMFNDKEVRAL'
            'FKNGRHWNVLVVIANQYVMDLTPDLRSSVDGVFLFRENNVTYRDKTYANFASVVPKKLYPTVMETV'
            'CQNYRCMFIDNTKATDNWHDSVFWYKAPYSKSAVAPFGARSYWKYACSKTGEEMPAVFDNVKILGD'
            'LLLKELPEAGEALVTYGGKDGPSDNEDGPSDDEDGPSDDEEGLSKDGVSEYYQSDLDD'
        ),
    )
    for seq in sequences:
      self.assertNotIn('\n', seq)

  def test_lazy_parse_fasta_string_error(self):
    fasta_iterator = parsers.lazy_parse_fasta_string(fasta_string='')

    with self.assertRaises(ValueError):
      list(fasta_iterator)

  def test_parse_fasta(self):
    fasta_str = self.test_data.load('test_db.fasta')
    seqs, descs = parsers.parse_fasta(fasta_str)
    self.assertLen(seqs, 34)
    self.assertLen(descs, 34)
    self.assertEqual(
        seqs[30],
        (
            'MLQNYAIVLGMAVAVAIWYFFKIEEEAPPGPNPPKPDPPKPDPPKMHMPKKKPHWMDPHL'
            'TGSQTVQYSRNRSMGDPIRGDLPIIPRDDGWFSTAANPAHTLHAGALSMIAPASTGGGLT'
            'VNKLISAYADKGNAMSGRHNSPSYYGSS'
        ),
    )
    self.assertEqual(
        descs[30],
        (
            'sp|Q6GZV5|020R_FRG3G Uncharacterized protein 020R OS=Frog virus 3 '
            '(isolate Goorha) GN=FV3-020R PE=4 SV=1'
        ),
    )
    for seq in seqs:
      self.assertNotRegex(seq, r'\s')

  @parameterized.parameters(
      ('> comment\nAAAAAA', ['AAAAAA']),
      ('> comment\nAAAAAA\n> comment 2\nBBBBBB', ['AAAAAA', 'BBBBBB']),
      ('> comment\nA\n> comment 2\nB\nC', ['A', 'BC']),
      ('> comment\n\nA\n\n> comment 2\nB\nC', ['A', 'BC']),
      ('> comment\n> comment 2\nA', ['', 'A']),
      ('> comment\nAAAAAA\nBBBBBB', ['AAAAAABBBBBB']),
      ('> comment\nAAAbbb', ['AAAbbb']),
      ('> comment\nAAA---', ['AAA---']),
      ('> comment\nAAA---bCc', ['AAA---bCc']),
      ('> comment\n\nAAA---bCc', ['AAA---bCc']),
      ('', []),
      ('\n', []),
  )
  def test_parse_a3m(self, a3m_file_string, expected_sequences):
    sequences, _ = parsers.parse_fasta(a3m_file_string)
    for exp_seq, seq in zip(expected_sequences, sequences, strict=True):
      self.assertEqual(exp_seq, seq)

  def test_a3m_to_sequences_with_comment(self):
    sequences, descriptions = parsers.parse_fasta(
        '> tst\nAAA\n> tst 2| bla hello\nB\n'
    )
    self.assertEqual(sequences, ['AAA', 'B'])
    self.assertEqual(descriptions, [' tst', ' tst 2| bla hello'])

  @parameterized.named_parameters(
      ('single_line_alignments', 'test_input.sto', 'expected_result.a3m', None),
      (
          'multi_line_alignments',
          'multiline_stockholm.sto',
          'multiline_a3m.a3m',
          None,
      ),
      (
          'single_line_alignments_truncated',
          'test_input.sto',
          'expected_result.a3m',
          4,
      ),
      (
          'multi_line_alignments_truncated',
          'multiline_stockholm.sto',
          'multiline_a3m.a3m',
          4,
      ),
  )
  def test_convert_stockholm_to_a3m(
      self, input_filename, output_filename, max_sequences
  ):
    with open(self.test_data.path(input_filename)) as f:
      a3m = parsers.convert_stockholm_to_a3m(
          stockholm=f, max_sequences=max_sequences
      )
    expected_a3m = self.test_data.load(output_filename)
    if max_sequences:
      expected_a3m = (
          '\n'.join(expected_a3m.splitlines()[: max_sequences * 2]) + '\n'
      )
    self.assertEqual(a3m, expected_a3m)

  def test_convert_stockholm_to_a3m_keep_gaps(self):
    sto = """# STOCKHOLM 1.0
#=GS 9nh9_Z/1-39 DE [subseq from] mol:unk length:40 foo bar baz
9nh9_Z/1-39         AAAAACGGGGUUU--UUUUUUUAAAAAAGGGGGGGGGGGGG-
//"""
    a3m = parsers.convert_stockholm_to_a3m(
        stockholm=io.StringIO(sto), remove_first_row_gaps=False
    )
    expected_a3m = """>9nh9_Z/1-39 [subseq from] mol:unk length:40 foo bar baz
AAAAACGGGGUUU--UUUUUUUAAAAAAGGGGGGGGGGGGG-
"""
    self.assertEqual(a3m, expected_a3m)

  def test_convert_stockholm_to_a3m_removes_dots(self):
    sto = """# STOCKHOLM 1.0
#=GS humanchr1_frag/1-60          DE [subseq from] humanchr1_frag
humanchr1_frag/1-60                  CCAAAAATACGAAAAA...GTAGCCAGGCATAATGGCACACATGTGTGGTC.CC.GCTGCTTGGG
//"""
    a3m = parsers.convert_stockholm_to_a3m(stockholm=io.StringIO(sto))
    expected_a3m = """>humanchr1_frag/1-60 [subseq from] humanchr1_frag
CCAAAAATACGAAAAAGTAGCCAGGCATAATGGCACACATGTGTGGTCCCGCTGCTTGGG
"""
    self.assertEqual(a3m, expected_a3m)

  def test_convert_stockholm_to_a3m_linewidth(self):
    sto = """# STOCKHOLM 1.0
#=GS 6v0c_A/2-68    DE [subseq from] mol:protein length:832  Transport protein YebT
6v0c_A/2-68            -------------------------------------NTVTIDFMSADGIVPGrTPVRYQGVEVGTVQDISLSDDLRKIEVKVSIKSDMkD.A.LREETQFWLVT...........P----------------....-------------------------
//"""
    a3m = parsers.convert_stockholm_to_a3m(
        stockholm=io.StringIO(sto), remove_first_row_gaps=False, linewidth=60
    )
    expected_a3m = """>6v0c_A/2-68 [subseq from] mol:protein length:832  Transport protein YebT
-------------------------------------NTVTIDFMSADGIVPGrTPVRYQ
GVEVGTVQDISLSDDLRKIEVKVSIKSDMkDALREETQFWLVTP----------------
-------------------------
"""
    self.assertEqual(a3m, expected_a3m)

  def test_convert_stockholm_to_a3m_empty(self):
    self.assertEmpty(parsers.convert_stockholm_to_a3m(io.StringIO('')))

  def test_convert_a3m_to_stockholm(self):
    a3m = self.test_data.load('test_input.a3m')
    expected_stockholm = self.test_data.load('expected_result.sto')

    stockholm = parsers.convert_a3m_to_stockholm(a3m=a3m)
    self.assertEqual(stockholm, expected_stockholm)

  def test_convert_stockholm_to_a3m_zero_linewidth(self):
    with self.assertRaisesRegex(ValueError, 'linewidth must be > 0 or None'):
      parsers.convert_stockholm_to_a3m(stockholm=io.StringIO('aa'), linewidth=0)

  def test_convert_a3m_to_stockholm_with_tabs(self):
    a3m = """
    >101
LDSPDRPWNPPTFSPALLVVTE
>UniRef100_UPI001E69966C\t167\t0.705\t5.170E-44\t0\t145\t146\t24\t169\t233
LDSPDRPWSPLTFSPAQLTVHE
"""
    stockholm = parsers.convert_a3m_to_stockholm(a3m=a3m)
    exp_stockholm = """# STOCKHOLM 1.0

#=GS 101_0 DE <EMPTY>
#=GS UniRef100_UPI001E69966C_1 DE 167 0.705 5.170E-44 0 145 146 24 169 233

101_0                     LDSPDRPWNPPTFSPALLVVTE
UniRef100_UPI001E69966C_1 LDSPDRPWSPLTFSPAQLTVHE
#=GC RF                   xxxxxxxxxxxxxxxxxxxxxx
//"""

    self.assertEqual(stockholm, exp_stockholm)


if __name__ == '__main__':
  absltest.main()
