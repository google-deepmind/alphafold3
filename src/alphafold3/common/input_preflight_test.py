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

from alphafold3.common import input_preflight


class ProteinInputStatusTest(parameterized.TestCase):

  def test_unset_msa_and_templates_search_both(self):
    status = input_preflight.protein_input_status('A', None, None, None)
    self.assertIs(
        status, input_preflight.ProteinInputStatus.SEARCH_MSA_AND_TEMPLATES
    )

  def test_empty_templates_search_msa_only(self):
    status = input_preflight.protein_input_status('A', None, None, [])
    self.assertIs(status, input_preflight.ProteinInputStatus.SEARCH_MSA)

  def test_supplied_msa_searches_templates(self):
    status = input_preflight.protein_input_status('A', '>q\nA', '>q\nA', None)
    self.assertIs(status, input_preflight.ProteinInputStatus.SEARCH_TEMPLATES)

  def test_supplied_msa_and_templates(self):
    status = input_preflight.protein_input_status('A', '', '', [])
    self.assertIs(status, input_preflight.ProteinInputStatus.SUPPLIED)

  @parameterized.parameters(
      ('>q\nA', None, None),
      (None, '>q\nA', None),
      ('>q\nA', None, []),
      (None, '>q\nA', []),
      (None, None, ['template']),
  )
  def test_partial_fields_use_pipeline_error(self, unpaired, paired, templates):
    with self.assertRaisesRegex(ValueError, 'set only partially'):
      input_preflight.protein_input_status('A', unpaired, paired, templates)

  def test_descriptions(self):
    self.assertEqual(
        input_preflight.describe_protein_input('A', None, None, None),
        'protein A: MSA and templates require the data pipeline',
    )
    self.assertEqual(
        input_preflight.describe_protein_input('B', None, None, []),
        'protein B: MSA requires the data pipeline; templates disabled',
    )
    self.assertEqual(
        input_preflight.describe_protein_input('C', '>q\nA', '', None),
        'protein C: MSA supplied; templates require the data pipeline',
    )
    self.assertEqual(
        input_preflight.describe_protein_input('D', '', '', []),
        'protein D: MSA empty; templates disabled',
    )
    self.assertEqual(
        input_preflight.describe_protein_input('E', '>q\nA', '', ['t']),
        'protein E: MSA custom; templates 1 supplied',
    )
    self.assertEqual(
        input_preflight.describe_rna_msa('R', None),
        'rna R: MSA requires the data pipeline',
    )
    self.assertEqual(
        input_preflight.describe_rna_msa('R', ''),
        'rna R: MSA supplied empty',
    )
    self.assertEqual(
        input_preflight.describe_rna_msa('R', '>q\nA'),
        'rna R: MSA supplied',
    )


if __name__ == '__main__':
  absltest.main()
