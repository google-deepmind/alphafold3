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
from alphafold3.constants import mmcif_names
from alphafold3.data import msa_config


class MsaConfigTest(parameterized.TestCase):

  def test_invalid_crop_size(self):
    with self.assertRaisesWithLiteralMatch(
        ValueError, 'crop_size must be None or >= 2: 1'
    ):
      msa_config.RunConfig(
          config=msa_config.JackhmmerConfig(
              binary_path='/path/to/jackhmmer',
              database_config=msa_config.DatabaseConfig(
                  name='test',
                  path='/path/to/test',
              ),
              n_cpu=1,
              n_iter=1,
              e_value=1e-3,
              z_value=None,
              dom_z_value=None,
              max_sequences=5000,
          ),
          chain_poly_type=mmcif_names.PROTEIN_CHAIN,
          crop_size=1,
      )

  def test_run_config_invalid_chain_poly_type(self):
    with self.assertRaisesWithLiteralMatch(
        ValueError,
        'chain_poly_type must be one of'
        f' {mmcif_names.STANDARD_POLYMER_CHAIN_TYPES}: invalid_chain_poly_type',
    ):
      msa_config.RunConfig(
          config=msa_config.JackhmmerConfig(
              binary_path='/path/to/jackhmmer',
              database_config=msa_config.DatabaseConfig(
                  name='test',
                  path='/path/to/test',
              ),
              n_cpu=1,
              n_iter=1,
              e_value=1e-3,
              z_value=None,
              dom_z_value=None,
              max_sequences=5000,
          ),
          chain_poly_type='invalid_chain_poly_type',
          crop_size=2,
      )

  def test_template_tool_config_invalid_chain_poly_type(self):
    with self.assertRaisesWithLiteralMatch(
        ValueError,
        'chain_poly_type must be one of'
        f' {mmcif_names.STANDARD_POLYMER_CHAIN_TYPES}: invalid_chain_poly_type',
    ):
      msa_config.TemplateToolConfig(
          database_path='/path/to/database',
          chain_poly_type='invalid_chain_poly_type',
          hmmsearch_config=msa_config.HmmsearchConfig(
              hmmsearch_binary_path='/path/to/hmmsearch',
              hmmbuild_binary_path='/path/to/hmmbuild',
              n_cpu=8,
              e_value=1e-3,
              inc_e=1e-3,
              dom_e=1e-3,
              incdom_e=1e-3,
              alphabet='amino',
              filter_f1=0.1,
              filter_f2=0.1,
              filter_f3=0.1,
              filter_max=False,
          ),
      )


if __name__ == '__main__':
  absltest.main()
