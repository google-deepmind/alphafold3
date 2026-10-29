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

import tempfile

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.cpp import fasta_iterator


TEST_CASES = (
    ('>name\nAAAAAA', ['AAAAAA'], ['name']),
    ('>name \t\r\nAAAAAA \n\r BB ', ['AAAAAABB'], ['name']),
    ('>\nAAAA\n>\nbbbb', ['AAAA', 'bbbb'], ['', '']),
    ('>name', [''], ['name']),
    ('>name\n', [''], ['name']),
    ('>name\n\n', [''], ['name']),
    ('>tst\nAAA\n>tst 2| bla hi\nB\n', ['AAA', 'B'], ['tst', 'tst 2| bla hi']),
    ('>name\nAAAAAA\n>name 2\nBBBBB', ['AAAAAA', 'BBBBB'], ['name', 'name 2']),
    ('>name\nA\n>name 2\nB\nC', ['A', 'BC'], ['name', 'name 2']),
    ('>name\n\nA\n\n>name 2\nB\nC', ['A', 'BC'], ['name', 'name 2']),
    ('>name\n>name 2\nA', ['', 'A'], ['name', 'name 2']),
    ('>name\nAAAAAA\nBBBBBB', ['AAAAAABBBBBB'], ['name']),
    ('>name\nAAAbbb', ['AAAbbb'], ['name']),
    ('>name\nAAA---', ['AAA---'], ['name']),
    ('>name\nAAA---bCc', ['AAA---bCc'], ['name']),
    ('>name\n\nAAA---bCc', ['AAA---bCc'], ['name']),
)

BAD_TEST_CASES = (
    ('NOHEADER', [], []),
    ('NOHEADER\nSTILLNONE', [], []),
    ('', [], []),
    ('\n', [], []),
)


class FastaCppTest(parameterized.TestCase):

  @parameterized.parameters(*TEST_CASES, *BAD_TEST_CASES)
  def test_parse_fasta(self, fasta_file_string, expected_sequences, _):
    sequences = fasta_iterator.parse_fasta(fasta_file_string)
    self.assertSequenceEqual(sequences, expected_sequences)

  @parameterized.parameters(*TEST_CASES, *BAD_TEST_CASES)
  def test_parse_fasta_with_description(
      self, fasta_file_string, expected_sequences, expected_descriptions
  ):
    sequences, descriptions = fasta_iterator.parse_fasta_include_descriptions(
        fasta_file_string
    )
    self.assertSequenceEqual(sequences, expected_sequences)
    self.assertSequenceEqual(descriptions, expected_descriptions)

  @parameterized.parameters(*TEST_CASES)
  def test_lazy_parse_fasta_file(
      self, fasta_file_string, expected_sequences, expected_descriptions
  ):
    with tempfile.NamedTemporaryFile('wt') as f:
      f.write(fasta_file_string)
      f.flush()
      iterator = fasta_iterator.FastaFileIterator(f.name)

    sequences = []
    descriptions = []
    for seq, descr in iterator:
      sequences.append(seq)
      descriptions.append(descr)
    self.assertSequenceEqual(sequences, expected_sequences)
    self.assertSequenceEqual(descriptions, expected_descriptions)

  @parameterized.parameters(*BAD_TEST_CASES)
  def test_lazy_parse_fasta_file_error(self, fasta_file_string, *args):
    del args  # Unused.
    with tempfile.NamedTemporaryFile('wt') as f:
      f.write(fasta_file_string)
      f.flush()
      iterator = fasta_iterator.FastaFileIterator(f.name)

    with self.assertRaises(ValueError):
      for _ in iterator:
        pass

  def test_lazy_parse_fasta_bad_path(self):
    iterator = fasta_iterator.FastaFileIterator('/bad/path')

    with self.assertRaisesRegex(ValueError, 'Invalid FASTA file: /bad/path'):
      for _ in iterator:
        pass

  @parameterized.parameters(*TEST_CASES)
  def test_lazy_parse_fasta_string(
      self, fasta_string, expected_sequences, expected_descriptions
  ):
    sequences = []
    descriptions = []
    for seq, descr in fasta_iterator.FastaStringIterator(fasta_string):
      sequences.append(seq)
      descriptions.append(descr)

    self.assertSequenceEqual(sequences, expected_sequences)
    self.assertSequenceEqual(descriptions, expected_descriptions)

  @parameterized.parameters(*BAD_TEST_CASES)
  def test_lazy_parse_fasta_string_error(self, fasta_string, *args):
    del args  # Unused.
    iterator = fasta_iterator.FastaStringIterator(fasta_string)

    with self.assertRaises(ValueError):
      for _ in iterator:
        pass


if __name__ == '__main__':
  absltest.main()
