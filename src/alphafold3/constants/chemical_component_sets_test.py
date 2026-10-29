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
from alphafold3.constants import chemical_component_sets


class ChemicalComponentSetsTest(parameterized.TestCase):

  @parameterized.product(glycal_other=('FFX', 'LAK', 'CR1'))
  def test_glycal_other_in(self, glycal_other):
    self.assertIn(glycal_other, chemical_component_sets.GLYCAN_OTHER_LIGANDS)

  @parameterized.product(glycal_linking=('G6P', '46Z', '1GN'))
  def test_glycal_linking_in(self, glycal_linking):
    self.assertIn(
        glycal_linking, chemical_component_sets.GLYCAN_LINKING_LIGANDS
    )

  @parameterized.product(glycal_other=('G6P', '46Z', '1GN'))
  def test_glycal_other_not_in(self, glycal_other):
    self.assertNotIn(glycal_other, chemical_component_sets.GLYCAN_OTHER_LIGANDS)

  @parameterized.product(glycal_linking=('FFX', 'LAK', 'CR1'))
  def test_glycal_linking_not_in(self, glycal_linking):
    self.assertNotIn(
        glycal_linking, chemical_component_sets.GLYCAN_LINKING_LIGANDS
    )

if __name__ == '__main__':
  absltest.main()
