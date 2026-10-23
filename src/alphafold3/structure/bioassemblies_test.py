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
from alphafold3.structure import bioassemblies
from alphafold3.structure import mmcif
from alphafold3.structure import test_utils
import numpy as np


class BioassembliesTest(parameterized.TestCase):

  def test_oper(self):
    oper = bioassemblies.Operation(
        trans=np.array([5.0, 6.0, 7.0]),
        rot=np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    )
    coords = np.eye(3)
    transformed = oper.apply_to_coords(coords)
    expected = np.array([[5.0, 7.0, 7.0], [4.0, 6.0, 7.0], [5.0, 6.0, 8.0]])
    np.testing.assert_array_equal(transformed, expected)

  def test_transform(self):
    oper1 = bioassemblies.Operation(
        trans=np.array([5.0, 6.0, 7.0]),
        rot=np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    )
    oper2 = bioassemblies.Operation(
        trans=np.array([0.0, -6.0, 0.0]), rot=np.eye(3)
    )
    transform = bioassemblies.Transform(
        operations=[oper2, oper1],
        chain_ids=['A'],
        chain_id_rename_map={'A': 'A'},
    )
    coords = np.eye(3)
    transformed = transform.apply_to_coords(coords)
    expected = np.array([[5.0, 1.0, 7.0], [4.0, 0.0, 7.0], [5.0, 0.0, 8.0]])
    np.testing.assert_array_equal(transformed, expected)

  def test_bioassembly_data_from_cif(self):
    cif_str = test_utils.read_test_cif_file('2bfu')
    cif = mmcif.from_string(cif_str)
    bd = bioassemblies.BioassemblyData.from_mmcif(cif)
    self.assertLen(bd._operations, 62)
    np.testing.assert_array_equal(bd._operations['X0'].trans, np.zeros((3,)))
    np.testing.assert_array_equal(bd._operations['X0'].rot, np.eye(3))

    np.testing.assert_array_equal(bd._operations['1'].trans, np.zeros((3,)))
    np.testing.assert_array_equal(bd._operations['1'].rot, np.eye(3))

    np.testing.assert_array_equal(bd._operations['P'].trans, np.zeros((3,)))
    np.testing.assert_array_equal(
        bd._operations['P'].rot,
        np.array(
            [
                [0.30901699, -0.80901699, 0.50000000],
                [0.80901699, 0.50000000, 0.30901699],
                [-0.50000000, 0.30901699, 0.80901699],
            ],
            dtype=np.float32,
        ),
    )

  @parameterized.parameters(
      ('2bfu', list('123456')),
      ('7a9y', list('12')),
      ('6tug_trimmed', list('1234')),
  )
  def test_bioassembly_data_assembly_ids(self, pdb_id, expected_ids):
    cif_str = test_utils.read_test_cif_file(pdb_id)
    cif = mmcif.from_string(cif_str)
    bd = bioassemblies.BioassemblyData.from_mmcif(cif)
    self.assertCountEqual(bd.assembly_ids, expected_ids)

  def test_asym_id_by_assembly_chain_id(self):
    cif_str = test_utils.read_test_cif_file('2bfu')
    cif = mmcif.from_string(cif_str)
    bd = bioassemblies.BioassemblyData.from_mmcif(cif)
    # Assembly 3 has 5 copies of chains A and B.
    asym_id_by_assembly_chain_id = bd.asym_id_by_assembly_chain_id('3')
    self.assertDictEqual(
        asym_id_by_assembly_chain_id,
        {
            'A': 'A',
            'B': 'B',
            'C': 'A',
            'D': 'B',
            'E': 'A',
            'F': 'B',
            'G': 'A',
            'H': 'B',
            'I': 'A',
            'J': 'B',
        },
    )

  def test_assembly_chain_ids_by_asym_id(self):
    cif_str = test_utils.read_test_cif_file('2bfu')
    cif = mmcif.from_string(cif_str)
    bd = bioassemblies.BioassemblyData.from_mmcif(cif)
    # Assembly 3 has 5 copies of chains A and B.
    assembly_chain_ids_by_asym_id = bd.assembly_chain_ids_by_asym_id('3')
    self.assertDictEqual(
        assembly_chain_ids_by_asym_id,
        {
            'A': {'A', 'C', 'E', 'G', 'I'},
            'B': {'B', 'D', 'F', 'H', 'J'},
        },
    )

  @parameterized.parameters(('2bfu', '1'), ('7a9y', '1'), ('6tug_trimmed', '1'))
  def test_bioassembly_data_get_default_assembly_id(self, pdb_id, expected_id):
    cif_str = test_utils.read_test_cif_file(pdb_id)
    cif = mmcif.from_string(cif_str)
    bd = bioassemblies.BioassemblyData.from_mmcif(cif)
    self.assertEqual(bd.get_default_assembly_id(), expected_id)

  def test_bioassembly_data_get_assembly_info(self):
    cif_str = test_utils.read_test_cif_file('2bfu')
    cif = mmcif.from_string(cif_str)
    bd = bioassemblies.BioassemblyData.from_mmcif(cif)
    info = bd.get_assembly_info('4')
    expected = {
        'id': '4',
        'details': 'icosahedral 23 hexamer',
        'method_details': '?',
        'oligomeric_details': 'dodecameric',
        'oligomeric_count': '12',
    }
    self.assertDictEqual(info, expected)

  def test_bioassembly_data_get_transforms_2bfu(self):
    cif_str = test_utils.read_test_cif_file('2bfu')
    cif = mmcif.from_string(cif_str)
    bd = bioassemblies.BioassemblyData.from_mmcif(cif)
    transforms = bd.get_transforms('6')
    self.assertLen(transforms, 5)
    seen_chain_ids = set()
    for transform in transforms:
      self.assertLen(transform.operations, 2)
      self.assertCountEqual(transform.chain_ids, ['A', 'B'])
      renamed_chains = [
          transform.chain_id_rename_map[c] for c in transform.chain_ids
      ]
      self.assertLen(set(renamed_chains), len(renamed_chains))
      self.assertEmpty(seen_chain_ids & set(renamed_chains))
      seen_chain_ids |= set(renamed_chains)

  def test_bioassembly_data_get_transforms_7a9y(self):
    cif_str = test_utils.read_test_cif_file('7a9y')
    cif = mmcif.from_string(cif_str)
    bd = bioassemblies.BioassemblyData.from_mmcif(cif)
    transforms = bd.get_transforms('1')
    self.assertLen(transforms, 2)
    seen_chain_ids = set()
    all_expected_chain_ids = [['A', 'C', 'D', 'F'], ['B', 'E', 'G']]
    for transform, expected_chain_ids in zip(
        sorted(transforms, key=lambda t: t.chain_ids[0]),
        all_expected_chain_ids,
        strict=True,
    ):
      self.assertLen(transform.operations, 1)
      self.assertCountEqual(transform.chain_ids, expected_chain_ids)
      renamed_chains = [
          transform.chain_id_rename_map[c] for c in transform.chain_ids
      ]
      self.assertLen(set(renamed_chains), len(renamed_chains))
      self.assertEmpty(seen_chain_ids & set(renamed_chains))
      seen_chain_ids |= set(renamed_chains)

  def test_bioassembly_data_get_transforms_6tug(self):
    cif_str = test_utils.read_test_cif_file('6tug_trimmed')
    cif = mmcif.from_string(cif_str)
    bd = bioassemblies.BioassemblyData.from_mmcif(cif)
    transforms = bd.get_transforms('1')
    self.assertLen(transforms, 1)
    transform = transforms[0]
    self.assertLen(transform.operations, 1)
    self.assertCountEqual(
        transform.chain_ids,
        ['A', 'B', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'KA', 'LA'],
    )
    renamed_chains = [
        transform.chain_id_rename_map[c] for c in transform.chain_ids
    ]
    self.assertLen(set(renamed_chains), len(renamed_chains))

  def test_bioassembly_data_to_mmcif_dict(self):
    cif_str = test_utils.read_test_cif_file('2bfu')
    cif = mmcif.from_string(cif_str)
    bd = bioassemblies.BioassemblyData.from_mmcif(cif)
    cif_dict = bd.to_mmcif_dict()

    for k, v in cif_dict.items():
      # Rows might be reordered so we just check the element counts
      # are equal to be order-invariant.
      self.assertCountEqual(v, cif[k], msg=k)

  @parameterized.parameters(
      ('2bfu', {'A': 'Z'}, {'A', 'B'}, ['A,B'] * 6, ['Z,B'] * 6),
      ('2bfu', {}, {'A', 'B'}, ['A,B'] * 6, ['A,B'] * 6),
      (
          '7a9y',
          {'A': 'X', 'B': 'ZZ'},
          {'A', 'B', 'C', 'D', 'E', 'F', 'G'},
          ['B,E,G', 'A,C,D,F', 'A,C,D,F', 'B,E,G'],
          ['ZZ,E,G', 'X,C,D,F', 'X,C,D,F', 'ZZ,E,G'],
      ),
  )
  def test_rename_label_asym_ids(
      self,
      target_name,
      mapping,
      present_chains,
      exp_old_asym_ids,
      exp_new_asym_ids,
  ):
    cif_str = test_utils.read_test_cif_file(target_name)
    cif = mmcif.from_string(cif_str)
    old_bd = bioassemblies.BioassemblyData.from_mmcif(cif)
    old_cif_dict = old_bd.to_mmcif_dict()
    self.assertSequenceEqual(
        old_cif_dict['_pdbx_struct_assembly_gen.asym_id_list'], exp_old_asym_ids
    )

    new_bd = old_bd.rename_label_asym_ids(
        mapping, present_chains=present_chains
    )
    new_cif_dict = new_bd.to_mmcif_dict()
    self.assertSequenceEqual(
        new_cif_dict['_pdbx_struct_assembly_gen.asym_id_list'], exp_new_asym_ids
    )

  def test_rename_label_asym_ids_error(self):
    cif_str = test_utils.read_test_cif_file('2bfu')
    cif = mmcif.from_string(cif_str)
    old_bd = bioassemblies.BioassemblyData.from_mmcif(cif)
    with self.assertRaises(ValueError):
      # Chain B already exists in the BioassemblyData.
      old_bd.rename_label_asym_ids({'A': 'B'}, present_chains={'A', 'B'})


if __name__ == '__main__':
  absltest.main()
