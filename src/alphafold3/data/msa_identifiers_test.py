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
from alphafold3.data import msa_identifiers
from alphafold3.data import parsers


EMPTY_IDENTIFIER = msa_identifiers.Identifiers(species_id='')

BFD_IDENTIFIERS = [
    EMPTY_IDENTIFIER,
    msa_identifiers.Identifiers('9RHOB'),
    msa_identifiers.Identifiers('9GAMM'),
    msa_identifiers.Identifiers('9BACL'),
    msa_identifiers.Identifiers('9BACL'),
    msa_identifiers.Identifiers('9BACL'),
    msa_identifiers.Identifiers('9BACL'),
    msa_identifiers.Identifiers('9FIRM'),
    msa_identifiers.Identifiers('9BACL'),
    EMPTY_IDENTIFIER,
]


class MsaIdentifiersTest(parameterized.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.test_data = data.Data(resources.ROOT / 'data/testdata')

  @parameterized.named_parameters(
      dict(
          testcase_name='uniprot_trembl',
          identifier_str='tr|A0A146SKV9|A0A146SKV9_FUNHE',
          species_id='FUNHE',
      ),
      dict(
          testcase_name='uniprot_swissprot',
          identifier_str='sp|P0C2L1|A3X1_LOXLA',
          species_id='LOXLA',
      ),
  )
  def test_get_identifiers(self, identifier_str, species_id):
    identifiers = msa_identifiers.get_identifiers(identifier_str)
    self.assertEqual(
        identifiers, msa_identifiers.Identifiers(species_id=species_id)
    )

  @parameterized.named_parameters(
      dict(
          testcase_name='bfd',
          filename='bfd_uniclust',
          ext='a3m',
          expected_identifiers=BFD_IDENTIFIERS,
      ),
      dict(
          testcase_name='small_bfd',
          filename='small_bfd',
          ext='sto',
          expected_identifiers=BFD_IDENTIFIERS,
      ),
      dict(
          testcase_name='uniref',
          filename='uniref',
          ext='sto',
          expected_identifiers=[EMPTY_IDENTIFIER] * 10,
      ),
      dict(
          testcase_name='mgnify',
          filename='mgnify',
          ext='sto',
          expected_identifiers=[EMPTY_IDENTIFIER] * 10,
      ),
  )
  def test_real_msas(self, filename, ext, expected_identifiers):
    msa_string = self.test_data.load(
        f'casp14_data_outputs/T0949/{filename}.{ext}'
    )
    if ext == 'sto':
      _, descriptions = parsers.parse_fasta(
          parsers.convert_stockholm_to_a3m(io.StringIO(msa_string))
      )
    elif ext == 'a3m':
      _, descriptions = parsers.parse_fasta(msa_string)
    else:
      self.fail(f'extension {ext} not recognised.')
    identifiers = [
        msa_identifiers.get_identifiers(desc) for desc in descriptions[:10]
    ]
    for identifier, expected in zip(identifiers, expected_identifiers):
      self.assertEqual(identifier, expected)


if __name__ == '__main__':
  absltest.main()
