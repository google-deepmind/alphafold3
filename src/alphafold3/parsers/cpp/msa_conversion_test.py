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
from alphafold3.cpp import msa_conversion


class MsaConversionTest(parameterized.TestCase):

  @parameterized.named_parameters(
      dict(testcase_name='empty', a3m=[], expected_sto=[]),
      dict(testcase_name='no_insertions', a3m=['AAAA'], expected_sto=['AAAA']),
      dict(testcase_name='all_insertions', a3m=['aaaa'], expected_sto=['AAAA']),
      dict(testcase_name='single_row', a3m=['-aA-'], expected_sto=['-AA-']),
      dict(
          testcase_name='no_overlap',
          a3m=[
              'aaa---',
              '-aaa--',
              '--aaa-',
              '---aaa',
          ],
          expected_sto=[
              'AAA------------',
              '----AAA--------',
              '--------AAA----',
              '------------AAA',
          ]),
      dict(
          testcase_name='small_example',
          a3m=[
              'abCD',
              'CgD',
              'fCDa',
          ],
          expected_sto=[
              'ABC-D-',
              '--CGD-',
              'F-C-DA',
          ]),
      dict(
          testcase_name='complete_overlap',
          a3m=['AAAA', 'AAAA'],
          expected_sto=['AAAA', 'AAAA']))
  def test_convert_a3m_to_stockholm(self, a3m, expected_sto):
    observed_sto = msa_conversion.convert_a3m_to_stockholm(a3m_sequences=a3m)
    self.assertListEqual(expected_sto, observed_sto)

  def test_convert_a3m_to_stockholm_bad_a3m(self):
    a3m = ['AAAA', 'BBB']
    with self.assertRaises(ValueError):
      msa_conversion.convert_a3m_to_stockholm(a3m)

  def test_align_sequence_to_gapless_query(self):
    aligned = msa_conversion.align_sequence_to_gapless_query(
        sequence='AB--E', query_sequence='A--DE')
    self.assertEqual(aligned, 'Ab-E')

  def test_align_sequence_to_gapless_query_bad_length(self):
    with self.assertRaises(ValueError):
      msa_conversion.align_sequence_to_gapless_query(
          sequence='AAA', query_sequence='AAAA')


if __name__ == '__main__':
  absltest.main()
