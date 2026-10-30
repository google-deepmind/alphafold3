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
import datetime
from typing import Any

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3 import structure
from alphafold3.constants import atom_types
from alphafold3.constants import chemical_components
from alphafold3.constants import mmcif_names
from alphafold3.constants import residue_names
from alphafold3.cpp import cif_dict
from alphafold3.structure import structure as structure_lib
from alphafold3.structure import structure_tables
from alphafold3.structure import test_utils
import numpy as np


class StructureTest(test_utils.StructureTestCase):

  def test_rename_auth_asym_ids_fails_if_not_unique(self):
    struc = structure.from_mmcif(
        test_utils.read_test_cif_file('2pv7'), include_bonds=False
    )
    with self.assertRaisesRegex(
        ValueError, 'The new polymer auth_asym_ids are not unique'
    ):
      # This should fail because it maps two polymers to the same auth ID.
      _ = struc.rename_auth_asym_ids({'A': 'B'})

  def test_split_by_chain(self):
    original_struc = structure.from_mmcif(
        test_utils.read_test_cif_file('2pv7'), include_bonds=False
    )
    chain_strucs = original_struc.split_by_chain()
    concat_struc = structure.concat(chain_strucs, assign_unique_chain_ids=False)

    self.assertLen(chain_strucs, original_struc.num_chains)
    self.assertEqual(
        sum(s.num_residues(count_unresolved=True) for s in chain_strucs),
        original_struc.num_residues(count_unresolved=True),
    )
    self.assertEqual(
        sum(s.num_residues(count_unresolved=False) for s in chain_strucs),
        original_struc.num_residues(count_unresolved=False),
    )
    self.assertEqual(
        sum(s.num_atoms for s in chain_strucs), original_struc.num_atoms
    )
    # Ignore all global fields that concat removes.
    ignore_fields = (
        'name',
        'release_date',
        'resolution',
        'structure_method',
        'bioassembly_data',
    )
    self.assertStructuresEqual(
        concat_struc, original_struc, ignore_fields=ignore_fields
    )

  def make_test_structure(
      self,
      res_id: np.ndarray,
      name: str | None = None,
      release_date: datetime.date | None = None,
      all_residues: Any = None,
      bonds: structure.Bonds | None = None,
      chain_id: np.ndarray | None = None,
      chain_type: np.ndarray | None = None,
      res_name: np.ndarray | None = None,
      atom_key: np.ndarray | None = None,
      atom_name: np.ndarray | None = None,
      atom_element: np.ndarray | None = None,
      atom_x: np.ndarray | None = None,
      atom_y: np.ndarray | None = None,
      atom_z: np.ndarray | None = None,
      atom_b_factor: np.ndarray | None = None,
      atom_occupancy: np.ndarray | None = None,
      chemical_components_data: structure.ChemicalComponentsData | None = None,
      shuffle_atom_keys: bool = True,
  ) -> structure.Structure:
    """Function to create a test structure."""

    if name is None:
      name = 'None'

    # Allow users to specify only atom_x for tests.
    if atom_y is None:
      atom_y = atom_x
    if atom_z is None:
      atom_z = atom_x

    atoms, residues, chains = structure_tables.tables_from_atom_arrays(
        res_id=res_id,
        all_residues=all_residues,
        chain_id=chain_id,
        chain_type=chain_type,
        res_name=res_name,
        atom_key=atom_key,
        atom_name=atom_name,
        atom_element=atom_element,
        atom_x=atom_x,
        atom_y=atom_y,
        atom_z=atom_z,
        atom_b_factor=atom_b_factor,
        atom_occupancy=atom_occupancy,
    )

    # Randomise the keys to make sure we are not relying on the default ones.
    np.random.seed(42)
    new_chain_key = np.array(chains.key)
    np.random.shuffle(new_chain_key)
    new_residue_key = np.array(residues.key)
    np.random.shuffle(new_residue_key)

    residue_key_map = dict(zip(residues.key, new_residue_key))
    chain_key_map = dict(zip(chains.key, new_chain_key))

    new_residue_chain_key = np.vectorize(chain_key_map.get)(residues.chain_key)
    new_atom_chain_key = np.vectorize(chain_key_map.get)(atoms.chain_key)
    new_atom_res_key = np.vectorize(residue_key_map.get)(atoms.res_key)

    chains = chains.copy_and_update(key=new_chain_key)
    residues = residues.copy_and_update(
        key=new_residue_key, chain_key=new_residue_chain_key
    )
    if not shuffle_atom_keys:
      atoms = atoms.copy_and_update(
          chain_key=new_atom_chain_key,
          res_key=new_atom_res_key,
      )
    else:
      new_atom_key = np.array(atoms.key)
      np.random.shuffle(new_atom_key)
      atom_key_map = dict(zip(atoms.key, new_atom_key))

      atoms = atoms.copy_and_update(
          chain_key=new_atom_chain_key,
          res_key=new_atom_res_key,
          key=new_atom_key,
      )
      if bonds is not None:
        new_bonds_from_key = np.vectorize(atom_key_map.get)(bonds.from_atom_key)
        new_bonds_dest_key = np.vectorize(atom_key_map.get)(bonds.dest_atom_key)
        bonds = bonds.copy_and_update(
            from_atom_key=new_bonds_from_key,
            dest_atom_key=new_bonds_dest_key,
        )

    return structure.Structure(
        name=name,
        release_date=release_date,
        chemical_components_data=chemical_components_data,
        atoms=atoms,
        chains=chains,
        residues=residues,
        bonds=bonds if bonds is not None else structure.Bonds.make_empty(),
    )

  def test_init_without_string_array(self):
    res_name = np.array(['ALA', 'GLY', 'ALA'], dtype=object)
    chain_id = np.array(['A', 'A', 'B'], dtype=object)
    chain_type = np.array(['polypeptide(L)'] * 3, dtype=object)
    struc = self.make_test_structure(
        name='test_structure',
        res_name=res_name,
        res_id=np.array([1, 2, 1], dtype=np.int32),
        chain_id=chain_id,
        chain_type=chain_type,
    )
    self.assertEqual(struc.name, 'test_structure')
    np.testing.assert_array_equal(struc.res_name, res_name)
    np.testing.assert_array_equal(struc.chain_id, chain_id)
    np.testing.assert_array_equal(struc.chain_type, chain_type)

  def test_get_item(self):
    res_id = np.array([0, 1, 2], dtype=np.int32)
    struc = self.make_test_structure(name='test_structure', res_id=res_id)
    self.assertEqual(struc['name'], 'test_structure')
    np.testing.assert_array_equal(struc['res_id'], res_id)
    with self.assertRaises((KeyError, AttributeError)):
      _ = struc['invalid_field']

  def test_num_elements(self):
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(['A', 'A', 'A', 'A', 'A', 'B'], dtype=object),
        res_id=np.array([0, 0, 1, 2, 3, 1], dtype=np.int32),
        atom_x=np.array([0.0, 0.0, 1.0, 1.5, 2.0, 3.0], dtype=np.float32),
    )
    self.assertEqual(struc.num_atoms, 6)
    self.assertEqual(struc.num_residues(count_unresolved=False), 5)
    self.assertEqual(struc.num_chains, 2)

  def test_iter_atoms(self):
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(['A', 'A', 'A'], dtype=object),
        res_id=np.array([0, 1, 2], dtype=np.int32),
        atom_x=np.array([[[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]]], dtype=np.float32),
    )  # 1 x 2 x 3
    expected = [
        {'chain_id': 'A', 'res_id': 0, 'atom_x': np.array([[0.0, 3.0]])},
        {'chain_id': 'A', 'res_id': 1, 'atom_x': np.array([[1.0, 4.0]])},
        {'chain_id': 'A', 'res_id': 2, 'atom_x': np.array([[2.0, 5.0]])},
    ]
    for atom, expected_atom in zip(struc.iter_atoms(), expected, strict=True):
      for field in expected_atom:
        if field == 'atom_x':
          np.testing.assert_array_equal(atom[field], expected_atom[field])
        else:
          self.assertEqual(atom[field], expected_atom[field])

  def test_iter_residues(self):
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(['A', 'A', 'A', 'A', 'A', 'B'], dtype=object),
        res_id=np.array([0, 0, 1, 1, 2, 1], dtype=np.int32),
        atom_x=np.array([0.0, 0.0, 1.0, 1.5, 2.0, 3.0], dtype=np.float32),
    )
    expected = [
        {'chain_id': 'A', 'res_id': 0},
        {'chain_id': 'A', 'res_id': 1},
        {'chain_id': 'A', 'res_id': 2},
        {'chain_id': 'B', 'res_id': 1},
    ]
    for atom, expected_atom in zip(
        struc.iter_residues(), expected, strict=True
    ):
      for field in expected_atom:
        self.assertEqual(atom[field], expected_atom[field])

  def test_iter_chains(self):
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(['A', 'A', 'A', 'A', 'A', 'B'], dtype=object),
        res_id=np.array([0, 0, 1, 1, 2, 1], dtype=np.int32),
        atom_x=np.array([0.0, 0.0, 1.0, 1.5, 2.0, 3.0], dtype=np.float32),
    )
    expected = [{'chain_id': 'A'}, {'chain_id': 'B'}]
    for chain, expected_chain in zip(
        struc.iter_chains(), expected, strict=True
    ):
      for field in expected_chain:
        self.assertEqual(chain[field], expected_chain[field])

  def test_iter_chain_ranges(self):
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(['A', 'A', 'A', 'C', 'B', 'B'], dtype=object),
        res_id=np.array([4, 5, 6, 1, 2, 5], dtype=np.int32),
    )
    expected = [(0, 3), (3, 4), (4, 6)]
    self.assertEqual(list(struc.iter_chain_ranges()), expected)

  def test_iter_residue_ranges(self):
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(['A', 'A', 'A', 'B', 'B', 'B', 'B'], dtype=object),
        res_id=np.array([1, 2, 2, 2, 2, 3, 4], dtype=np.int32),
    )
    expected = [(0, 1), (1, 3), (3, 5), (5, 6), (6, 7)]
    self.assertEqual(list(struc.iter_residue_ranges()), expected)

  def test_iter_ranges_empty(self):
    struc = structure.make_empty_structure()
    self.assertEmpty(list(struc.iter_residue_ranges()))
    self.assertEmpty(list(struc.iter_chain_ranges()))

  def test_iter_ranges_single_residue(self):
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(['A', 'A', 'A'], dtype=object),
        res_id=np.array([1, 1, 1], dtype=np.int32),
        atom_name=np.array(['C', 'CA', 'N'], dtype=object),
    )
    self.assertEqual(list(struc.iter_residue_ranges()), [(0, 3)])
    self.assertEqual(list(struc.iter_chain_ranges()), [(0, 3)])

  def test_iter_bonds(self):
    struc = self.make_test_structure(
        name='test_structure',
        bonds=structure.Bonds(
            key=np.array([1, 2]),
            type=np.array(['type1', 'type2'], dtype=object),
            role=np.array(['role1', 'role2'], dtype=object),
            from_atom_key=np.array([11, 2]),
            dest_atom_key=np.array([2, 3]),
        ),
        chain_id=np.array(['A', 'A', 'A', 'A', 'A', 'B'], dtype=object),
        res_id=np.array([0, 0, 0, 0, 1, 0], dtype=np.int32),
        atom_key=np.array([2, 11, 1, 0, 4, 3]),
        atom_name=np.array(['C', 'CA', 'N', 'O', 'C', 'C'], dtype=object),
        atom_x=np.array(
            [
                [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                [10.0, 20.0, 30.0, 40.0, 50.0, 60.0],
            ],
            dtype=np.float32,
        ),
        atom_b_factor=np.array(
            [0.5, 1.0, 1.5, 2.0, 2.5, 3.0], dtype=np.float32
        ),
        shuffle_atom_keys=False,
    )
    from_atoms, dest_atoms, bond_info_dicts = zip(
        *struc.iter_bonds(), strict=True
    )
    exp_from_atoms = [
        {
            'chain_id': 'A',
            'res_id': 0,
            'atom_key': 11,
            'atom_name': 'CA',
            'atom_x': np.array([2.0, 20.0]),
            'atom_b_factor': 1.0,
        },
        {
            'chain_id': 'A',
            'res_id': 0,
            'atom_key': 2,
            'atom_name': 'C',
            'atom_x': np.array([1.0, 10.0]),
            'atom_b_factor': 0.5,
        },
    ]
    for from_atom, exp_from_atoms in zip(
        from_atoms, exp_from_atoms, strict=True
    ):
      for field in exp_from_atoms:
        if field == 'atom_x':
          np.testing.assert_array_equal(from_atom[field], exp_from_atoms[field])
        else:
          self.assertEqual(from_atom[field], exp_from_atoms[field])
    exp_dest_atoms = [
        {
            'chain_id': 'A',
            'res_id': 0,
            'atom_key': 2,
            'atom_name': 'C',
            'atom_x': np.array([1.0, 10.0]),
            'atom_b_factor': 0.5,
        },
        {
            'chain_id': 'B',
            'res_id': 0,
            'atom_key': 3,
            'atom_name': 'C',
            'atom_x': np.array([6.0, 60.0]),
            'atom_b_factor': 3.0,
        },
    ]
    for dest_atom, exp_dest_atoms in zip(
        dest_atoms, exp_dest_atoms, strict=True
    ):
      for field in exp_dest_atoms:
        if field == 'atom_x':
          np.testing.assert_array_equal(dest_atom[field], exp_dest_atoms[field])
        else:
          self.assertEqual(dest_atom[field], exp_dest_atoms[field])
    exp_bond_info_dicts = [
        {
            'key': 1,
            'type': 'type1',
            'role': 'role1',
            'from_atom_key': 11,
            'dest_atom_key': 2,
        },
        {
            'key': 2,
            'type': 'type2',
            'role': 'role2',
            'from_atom_key': 2,
            'dest_atom_key': 3,
        },
    ]
    self.assertSequenceEqual(bond_info_dicts, exp_bond_info_dicts)

  def test_apply_atom_index_array(self):
    struc = self.make_test_structure(
        name='test_structure',
        res_id=np.array([0, 1, 2, 3, 4], dtype=np.int32),
        atom_x=np.array([[0.0, 1.0, 2.0, 3.0, 4.0]], dtype=np.float32),
    )

    int_indexed = struc._apply_atom_index_array(np.array([0, 2, 4]))
    np.testing.assert_array_equal(int_indexed.res_id, [0, 2, 4])
    np.testing.assert_array_equal(
        int_indexed.atom_x, np.array([[0.0, 2.0, 4.0]])
    )

    bool_indexed = struc._apply_atom_index_array(
        np.array([True, False, True, False, True])
    )
    np.testing.assert_array_equal(bool_indexed.res_id, [0, 2, 4])
    np.testing.assert_array_equal(
        bool_indexed.atom_x, np.array([[0.0, 2.0, 4.0]])
    )

  def test_by_residue(self):
    all_residues = {
        'A': [('ARG', 0), ('HIS', 1)],
        'B': [('HIS', 1), ('GLU', 2), ('SER', 3)],  # Missing (SER, 3).
    }
    struc = self.make_test_structure(
        chain_id=np.array(['A', 'A', 'A', 'B', 'B', 'B'], dtype=object),
        chain_type=np.array([mmcif_names.PROTEIN_CHAIN] * 6, dtype=object),
        res_id=np.array([0, 1, 1, 1, 2, 2], dtype=np.int32),
        res_name=np.array(
            ['ARG', 'HIS', 'HIS', 'HIS', 'GLU', 'GLU'], dtype=object
        ),
        atom_name=np.array(['C', 'C', 'CA', 'C', 'C', 'CA'], dtype=object),
        all_residues=all_residues,
    )
    by_res = struc.group_by_residue
    np.testing.assert_array_equal(by_res.chain_id, ['A', 'A', 'B', 'B'])
    np.testing.assert_array_equal(by_res.res_id, [0, 1, 1, 2])
    self.assertSequenceEqual(
        [r['res_id'] for r in by_res.iter_residues()], [0, 1, 1, 2]
    )

    self.assertDictEqual(
        by_res.chain_single_letter_sequence(include_missing_residues=False),
        {'A': 'RH', 'B': 'HE'},
    )
    self.assertDictEqual(
        by_res.chain_single_letter_sequence(include_missing_residues=True),
        {'A': 'RH', 'B': 'HES'},
    )
    self.assertDictEqual(by_res.all_residues, all_residues)

    np.testing.assert_array_equal(
        by_res.atom_name, np.array(['C', 'C', 'C', 'C'], dtype=object)
    )

  def test_by_chain(self):
    all_residues = {
        'A': [('ARG', 0), ('HIS', 1)],
        'B': [('HIS', 1), ('GLU', 2), ('SER', 3)],  # Missing (SER, 3).
    }
    struc = self.make_test_structure(
        chain_id=np.array(['A', 'A', 'A', 'B', 'B', 'B'], dtype=object),
        chain_type=np.array([mmcif_names.PROTEIN_CHAIN] * 6, dtype=object),
        res_id=np.array([0, 1, 1, 1, 2, 2], dtype=np.int32),
        res_name=np.array(
            ['ARG', 'HIS', 'HIS', 'HIS', 'GLU', 'GLU'], dtype=object
        ),
        atom_name=np.array(['C', 'C', 'CA', 'C', 'C', 'CA'], dtype=object),
        all_residues=all_residues,
    )
    by_chain = struc.group_by_chain
    np.testing.assert_array_equal(by_chain.chain_id, ['A', 'B'])

    self.assertDictEqual(by_chain.all_residues, all_residues)

    self.assertDictEqual(
        by_chain.chain_single_letter_sequence(include_missing_residues=False),
        {'A': 'R', 'B': 'H'},
    )
    self.assertDictEqual(
        by_chain.chain_single_letter_sequence(include_missing_residues=True),
        {'A': 'RH', 'B': 'HES'},
    )
    np.testing.assert_array_equal(
        by_chain.res_id, np.array([0, 1], dtype=np.int64)
    )
    np.testing.assert_array_equal(
        by_chain.atom_name, np.array(['C', 'C'], dtype=object)
    )

  def test_with_sorted_chains(self):
    chain_ids = ('A', 'C', 'AA', 'BA', 'CA', 'AB', 'BB', 'DDD')
    residue_ids = []
    strucs = []
    for chain_index, chain_id in enumerate(chain_ids):
      # Give residues unique IDs so we can test their ordering.
      res1_id = 10 * chain_index + 1
      res2_id = 10 * chain_index + 2
      res3_id = 10 * chain_index + 3
      residue_ids.extend((res1_id, res1_id, res2_id, res2_id))
      struc = self.make_test_structure(
          chain_id=np.array([chain_id] * 4, dtype=object),
          chain_type=np.array([mmcif_names.PROTEIN_CHAIN] * 4, dtype=object),
          res_id=np.array([res1_id, res1_id, res2_id, res2_id], dtype=np.int32),
          res_name=np.array(['ARG', 'ARG', 'HIS', 'HIS'], dtype=object),
          all_residues={
              chain_id: [('ARG', res1_id), ('HIS', res2_id), ('SER', res3_id)]
          },
          atom_name=np.array(['CA', 'N', 'CA', 'N'], dtype=object),
          atom_x=np.arange(4.0, dtype=np.float32) + 100 * chain_index,
      )
      strucs.append(struc)

    struc = structure.concat(strucs, assign_unique_chain_ids=False)
    reversed_struc = structure.concat(
        strucs[::-1], assign_unique_chain_ids=False
    )
    assert struc.chains[::-1] == reversed_struc.chains

    resorted = reversed_struc.with_sorted_chains

    with self.subTest('chain_fields'):
      self.assertSequenceEqual(struc.chains, resorted.chains)

    with self.subTest('residue_fields'):
      self.assertDictEqual(
          resorted.chain_single_letter_sequence(include_missing_residues=False),
          struc.chain_single_letter_sequence(include_missing_residues=False),
      )
      self.assertDictEqual(
          resorted.chain_single_letter_sequence(include_missing_residues=True),
          struc.chain_single_letter_sequence(include_missing_residues=True),
      )
      self.assertSequenceEqual(
          [r['res_id'] for r in resorted.iter_residues()],
          [r['res_id'] for r in struc.iter_residues()],
      )
      self.assertSequenceEqual(
          [
              r['res_id']
              for r in resorted.iter_residues(include_unresolved=True)
          ],
          [r['res_id'] for r in struc.iter_residues(include_unresolved=True)],
      )
      np.testing.assert_array_equal(resorted.res_id, struc.res_id)
      np.testing.assert_array_equal(resorted.res_name, struc.res_name)

    with self.subTest('atom_fields'):
      np.testing.assert_array_equal(struc.atom_x, resorted.atom_x)

  @parameterized.parameters([True, False])
  def test_filter(self, apply_per_element):
    struc = self.make_test_structure(
        name='test_filter',
        bonds=structure.Bonds(
            key=np.array([1, 2]),
            type=np.array(['?', '?'], dtype=object),
            role=np.array(['?', '?'], dtype=object),
            from_atom_key=np.array([11, 2]),
            dest_atom_key=np.array([2, 3]),
        ),
        chain_id=np.array(['A', 'A', 'A', 'A', 'A', 'B'], dtype=object),
        res_id=np.array([0, 0, 0, 0, 1, 0], dtype=np.int32),
        atom_key=np.array([2, 11, 1, 0, 4, 3]),
        atom_name=np.array(['C', 'CA', 'N', 'O', 'C', 'C'], dtype=object),
        atom_x=np.array(
            [
                [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                [10.0, 20.0, 30.0, 40.0, 50.0, 60.0],
            ],
            dtype=np.float32,
        ),
        atom_b_factor=np.array(
            [0.5, 1.0, 1.5, 2.0, 2.5, 3.0], dtype=np.float32
        ),
    )
    filtered = struc.filter(
        apply_per_element=apply_per_element,
        chain_id='A',
        res_id=0,
        atom_name=('CA', 'C'),
        atom_b_factor=lambda b: b <= 1.5,
    )
    self.assertEqual(filtered.name, 'test_filter')
    np.testing.assert_array_equal(filtered.chain_id, ['A', 'A'])
    np.testing.assert_array_equal(filtered.res_id, [0, 0])
    np.testing.assert_array_equal(filtered.atom_name, ['C', 'CA'])
    np.testing.assert_array_equal(
        filtered.atom_x, np.array([[1.0, 2.0], [10.0, 20.0]])
    )
    np.testing.assert_array_equal(filtered.atom_b_factor, [0.5, 1.0])
    np.testing.assert_array_equal(filtered.bonds.key, [1])

  def test_filter_and_restore_newly_unresolved_chain(self):
    struc = self.make_test_structure(
        name='test_filter',
        chain_id=np.array(['A', 'A', 'B'], dtype=object),
        chain_type=np.array(
            [
                mmcif_names.PROTEIN_CHAIN,
                mmcif_names.PROTEIN_CHAIN,
                mmcif_names.RNA_CHAIN,
            ],
            dtype=object,
        ),
        res_name=np.array(['ALA', 'GLY', 'A'], dtype=object),
        res_id=np.array([1, 2, 1], dtype=np.int32),
        atom_b_factor=np.array([2.0, 2.0, 1.0], dtype=np.float32),
        all_residues={
            'A': [('ALA', 1), ('GLY', 2)],
            'B': [('A', 1)],
        },
    )

    filtered = struc.filter(
        cascade_delete=structure.CascadeDelete.NONE,
        mask=struc.atom_b_factor > 1.5,
    )
    # Check that chain A shows up in all_residues but not chain_ids, res_names.
    self.assertCountEqual(filtered.all_residues, ['A', 'B'])
    np.testing.assert_array_equal(filtered.res_name, ['ALA', 'GLY'])
    self.assertSequenceEqual(
        filtered.all_residues['A'], [('ALA', 1), ('GLY', 2)]
    )
    np.testing.assert_array_equal(filtered.chain_id, ['A', 'A'])
    self.assertEqual(filtered.chain_type.shape, filtered.chain_id.shape)

  def test_filter_missing_residues(self):
    struc = self.make_test_structure(
        all_residues={
            'A': [('ALA', 1), ('ILE', 2), ('GLY', 3)],
            'B': [('ASN', 1), ('PRO', 2), ('MET', 3)],
            'C': [('UNK', 1)],
        },
        chain_id=np.array(['A', 'A', 'A', 'B', 'B', 'B', 'C'], dtype=object),
        res_id=np.array([1, 1, 3, 2, 3, 3, 1], dtype=np.int32),
        res_name=np.array(
            ['ALA', 'ALA', 'GLY', 'PRO', 'MET', 'MET', 'UNK'], dtype=object
        ),
        atom_name=np.array(['C', 'N', 'C', 'C', 'C', 'N', 'C'], dtype=object),
    )
    struc = struc.filter(
        np.array([True] * 6 + [False]),
        cascade_delete=structure.CascadeDelete.NONE,
    )
    with self.subTest('empty_chains_get_removed_atom_filter'):
      filtered = struc.filter(atom_name='C')
      self.assertEqual(filtered.chains, ('A', 'B'))
    with self.subTest('empty_chains_get_removed_chain_filter'):
      filtered = struc.filter(chain_id='A')
      self.assertEqual(filtered.chains, ('A',))
    with self.subTest('missing_residues_not_filtered_atom_filter'):
      filtered = struc.filter(atom_name='N')
      np.testing.assert_array_equal(filtered.res_id, [1, 3])
      self.assertEqual(filtered.num_residues(count_unresolved=False), 2)
      self.assertEqual(filtered.num_residues(count_unresolved=True), 6)
    with self.subTest('missing_residues_not_filtered_res_filter'):
      # NB: This behavior will change in structure v2, where directly filtering
      # residues will affect both resolved and unresolved residues.
      filtered = struc.filter(res_id=lambda i: i >= 2)
      np.testing.assert_array_equal(filtered.res_id, [3, 2, 3, 3])
      self.assertEqual(filtered.num_residues(count_unresolved=False), 3)
      self.assertEqual(filtered.num_residues(count_unresolved=True), 6)

  @parameterized.named_parameters(
      dict(
          testcase_name='none_no_filter',
          apply_filter_to_protein=False,
          cascade_delete=structure.CascadeDelete.NONE,
          expected_chains=['A', 'B', 'C', 'D'],
          expected_all_residues={
              'A': [('ALA', 1), ('ILE', 2), ('GLY', 3)],
              'B': [('ASN', 1), ('PRO', 2), ('MET', 3)],
              'C': [('DA', 1)],
              'D': [('UNK', 1), ('UNK', 2)],
          },
      ),
      dict(
          testcase_name='chains_no_filter',
          cascade_delete=structure.CascadeDelete.CHAINS,
          apply_filter_to_protein=False,
          expected_chains=['A', 'C', 'D'],
          expected_all_residues={
              'A': [('ALA', 1), ('ILE', 2), ('GLY', 3)],
              'C': [('DA', 1)],
              'D': [('UNK', 1), ('UNK', 2)],
          },
      ),
      dict(
          testcase_name='full_no_filter',
          cascade_delete=structure.CascadeDelete.FULL,
          apply_filter_to_protein=False,
          expected_chains=['A', 'C', 'D'],
          expected_all_residues={
              'A': [('ALA', 1), ('GLY', 3)],
              'C': [('DA', 1)],
              'D': [('UNK', 1), ('UNK', 2)],
          },
      ),
      dict(
          testcase_name='none_with_filter',
          apply_filter_to_protein=True,
          cascade_delete=structure.CascadeDelete.NONE,
          expected_chains=['A', 'B', 'C', 'D'],
          expected_all_residues={
              'A': [('ALA', 1), ('ILE', 2), ('GLY', 3)],
              'B': [('ASN', 1), ('PRO', 2), ('MET', 3)],
              'C': [('DA', 1)],
              'D': [('UNK', 1), ('UNK', 2)],
          },
      ),
      dict(
          testcase_name='chains_with_filter',
          cascade_delete=structure.CascadeDelete.CHAINS,
          apply_filter_to_protein=True,
          expected_chains=['A', 'D'],
          expected_all_residues={
              'A': [('ALA', 1), ('ILE', 2), ('GLY', 3)],
              'D': [('UNK', 1), ('UNK', 2)],
          },
      ),
      dict(
          testcase_name='full_with_filter',
          cascade_delete=structure.CascadeDelete.FULL,
          apply_filter_to_protein=True,
          expected_chains=['A', 'D'],
          expected_all_residues={
              'A': [('ALA', 1), ('GLY', 3)],
              'D': [('UNK', 1), ('UNK', 2)],
          },
      ),
  )
  def test_filter_cascade_delete(
      self,
      cascade_delete: structure.CascadeDelete,
      apply_filter_to_protein: bool,
      expected_chains: Sequence[str],
      expected_all_residues: Mapping[str, Sequence[tuple[str, int]]],
  ):
    chains = structure.Chains(
        key=np.array([0, 1, 2, 3], dtype=np.int64),
        id=np.array(['A', 'B', 'C', 'D'], dtype=object),
        type=np.array(
            [
                mmcif_names.PROTEIN_CHAIN,
                mmcif_names.PROTEIN_CHAIN,
                mmcif_names.DNA_CHAIN,
                mmcif_names.PROTEIN_CHAIN,
            ],
            dtype=object,
        ),
        auth_asym_id=np.array(['A', 'B', 'C', 'D'], dtype=object),
        entity_id=np.array(['1', '2', '3', '4'], dtype=object),
        entity_desc=np.array(['.', '.', '.', '.'], dtype=object),
    )
    residues = structure.Residues(
        key=np.array([0, 1, 2, 3, 4, 5, 6, 7, 8], dtype=np.int64),
        chain_key=np.array([0, 0, 0, 1, 1, 1, 2, 3, 3], dtype=np.int64),
        id=np.array([1, 2, 3, 1, 2, 3, 1, 1, 2], dtype=np.int64),
        name=np.array(
            ['ALA', 'ILE', 'GLY', 'ASN', 'PRO', 'MET', 'DA', 'UNK', 'UNK'],
            dtype=object,
        ),
        auth_seq_id=np.array(
            ['1', '2', '3', '1', '2', '3', '1', '1', '2'], dtype=object
        ),
        insertion_code=np.array(['?'] * 9, dtype=object),
    )
    atoms = structure.Atoms(
        key=np.array([0, 1, 2, 3, 4, 5, 6], dtype=np.int64),
        chain_key=np.array([0, 0, 0, 0, 2, 3, 3], dtype=np.int64),
        res_key=np.array([0, 0, 0, 2, 6, 7, 8], dtype=np.int64),
        name=np.array(['N', 'CA', 'C', 'CA', "C1'", 'CA', 'CA'], dtype=object),
        element=np.array(['N', 'C', 'C', 'C', 'C', 'C', 'C'], dtype=object),
        x=np.zeros(7, dtype=np.float32),
        y=np.zeros(7, dtype=np.float32),
        z=np.zeros(7, dtype=np.float32),
        b_factor=np.zeros(7, dtype=np.float32),
        occupancy=np.ones(7, dtype=np.float32),
    )

    struc = structure.Structure(
        chains=chains,
        residues=residues,
        atoms=atoms,
        bonds=structure.Bonds.make_empty(),
    )

    self.assertDictEqual(
        {
            'A': [('ALA', 1), ('ILE', 2), ('GLY', 3)],  # ILE is unresolved
            'B': [('ASN', 1), ('PRO', 2), ('MET', 3)],  # Chain B is unresolved
            'C': [('DA', 1)],
            'D': [('UNK', 1), ('UNK', 2)],  # Chain D is resolved but unk.
        },
        struc.all_residues,
    )

    if apply_filter_to_protein:
      struc = struc.filter(
          cascade_delete=cascade_delete, chain_type=mmcif_names.PROTEIN_CHAIN
      )
    else:
      struc = struc.filter(cascade_delete=cascade_delete)

    self.assertDictEqual(expected_all_residues, struc.all_residues)
    self.assertSequenceEqual(expected_chains, struc.chains)

  def test_filter_multidim(self):
    struc = self.make_test_structure(
        name='test_structure',
        res_id=np.array([0, 1, 2], dtype=np.int32),
        atom_b_factor=np.array(
            [[0.1, 0.2, 0.3], [1.0, 2.0, 3.0]], dtype=np.float32
        ),
    )
    with self.assertRaises(ValueError):
      struc.filter(atom_b_factor=lambda b: b < 1.5)

  def test_filter_with_mask(self):
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(['A', 'A', 'A', 'A', 'A', 'B'], dtype=object),
        res_id=np.array([0, 0, 0, 0, 1, 0], dtype=np.int32),
        atom_name=np.array(['C', 'CA', 'N', 'O', 'C', 'C'], dtype=object),
        atom_x=np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=np.float32),
        atom_b_factor=np.array(
            [0.5, 1.0, 1.5, 2.0, 2.5, 3.0], dtype=np.float32
        ),
    )
    filtered_mask_and_preds = struc.filter(
        np.array([True, False, True, False, False, True]), chain_id='A'
    )
    np.testing.assert_array_equal(
        filtered_mask_and_preds.atom_x, np.array([[1.0, 3.0]])
    )
    filtered_mask_only = struc.filter(
        np.array([True, False, True, False, False, True])
    )
    np.testing.assert_array_equal(
        filtered_mask_only.atom_x, np.array([[1.0, 3.0, 6.0]])
    )

  @parameterized.named_parameters(
      dict(
          testcase_name='single_chain',
          chain_pred='A',
          cascade_delete=structure_lib.CascadeDelete.CHAINS,
          exp_chains=['A'],
          exp_all_residues={'A': [('ALA', 1), ('CYS', 2)]},
          exp_atom_x=[[1.0], [11.0]],
          exp_bond_keys=[],
      ),
      dict(
          testcase_name='adjacent_chains_chains_cascade',
          chain_pred=('A', 'B'),
          cascade_delete=structure_lib.CascadeDelete.CHAINS,
          exp_chains=['A', 'B'],
          exp_all_residues={
              'A': [('ALA', 1), ('CYS', 2)],
              'B': [('GLY', 1), ('SER', 2)],
          },
          exp_atom_x=[[1.0, 2.0, 3.0], [11.0, 12.0, 13.0]],
          exp_bond_keys=[1],
      ),
      dict(
          testcase_name='adjacent_chains_full_cascade',
          chain_pred=('A', 'B'),
          cascade_delete=structure_lib.CascadeDelete.FULL,
          exp_chains=['A', 'B'],
          exp_all_residues={
              'A': [('ALA', 1)],
              'B': [('GLY', 1), ('SER', 2)],
          },
          exp_atom_x=[[1.0, 2.0, 3.0], [11.0, 12.0, 13.0]],
          exp_bond_keys=[1],
      ),
      dict(
          testcase_name='disjoint_chains',
          chain_pred=('A', 'D'),
          cascade_delete=structure_lib.CascadeDelete.CHAINS,
          exp_chains=['A', 'D'],
          exp_all_residues={
              'A': [('ALA', 1), ('CYS', 2)],
              'D': [('LEU', 1)],
          },
          exp_atom_x=[[1.0, 5.0], [11.0, 15.0]],
          exp_bond_keys=[3],
      ),
  )
  def test_filter_chain_spans(
      self,
      chain_pred,
      cascade_delete,
      exp_chains,
      exp_all_residues,
      exp_atom_x,
      exp_bond_keys,
  ):
    atom_x = np.array(
        [[1.0, 2.0, 3.0, 4.0, 5.0], [11.0, 12.0, 13.0, 14.0, 15.0]],
        dtype=np.float32,
    )
    struc = structure.from_atom_arrays(
        name='span_test',
        chain_id=np.array(['A', 'B', 'B', 'C', 'D'], dtype=object),
        res_id=np.array([1, 1, 2, 1, 1], dtype=np.int32),
        res_name=np.array(['ALA', 'GLY', 'SER', 'VAL', 'LEU'], dtype=object),
        all_residues={
            'A': [('ALA', 1), ('CYS', 2)],
            'B': [('GLY', 1), ('SER', 2)],
            'C': [('VAL', 1)],
            'D': [('LEU', 1)],
        },
        atom_key=np.array([10, 20, 30, 40, 50], dtype=np.int64),
        atom_name=np.array(['CA', 'CA', 'CA', 'CA', 'CA'], dtype=object),
        atom_x=atom_x,
        atom_y=atom_x,
        atom_z=atom_x,
        bond_table=structure.Bonds(
            key=np.array([1, 2, 3], dtype=np.int64),
            type=np.array(['covale', 'covale', 'covale'], dtype=object),
            role=np.array(['.', '.', '.'], dtype=object),
            from_atom_key=np.array([10, 20, 10], dtype=np.int64),
            dest_atom_key=np.array([20, 40, 50], dtype=np.int64),
        ),
    )
    filtered = struc.filter(chain_id=chain_pred, cascade_delete=cascade_delete)
    np.testing.assert_array_equal(
        filtered.chains_table.id,
        np.array(exp_chains, dtype=object),
        strict=True,
    )
    self.assertAllResiduesEqual(filtered.all_residues, exp_all_residues)
    np.testing.assert_array_equal(
        filtered.atom_x, np.array(exp_atom_x, dtype=np.float32), strict=True
    )
    np.testing.assert_array_equal(
        filtered.bonds_table.key,
        np.array(exp_bond_keys, dtype=np.int64),
        strict=True,
    )

  def test_filter_many_disjoint_chains(self):
    num_chains = 80
    coords = np.arange(num_chains, dtype=np.float32)
    struc = structure.from_atom_arrays(
        name='many_chains',
        chain_id=np.array([f'C{i}' for i in range(num_chains)], dtype=object),
        res_id=np.ones(num_chains, dtype=np.int32),
        atom_x=coords,
        atom_y=coords,
        atom_z=coords,
    )
    target_chains = tuple(f'C{i}' for i in range(0, num_chains, 2))
    filtered = struc.filter(chain_id=target_chains)
    np.testing.assert_array_equal(
        filtered.chains_table.id,
        np.array(target_chains, dtype=object),
        strict=True,
    )
    np.testing.assert_array_equal(
        filtered.atom_x,
        np.arange(0, num_chains, 2, dtype=np.float32),
        strict=True,
    )

  @parameterized.named_parameters(
      dict(
          testcase_name='mask_only',
          kwargs={},
          expected=[2.0],
      ),
      dict(
          testcase_name='mask_with_chain_predicate',
          kwargs={'chain_id': 'A'},
          expected=[2.0, 3.0],
      ),
  )
  def test_filter_does_not_mutate_input_mask(self, kwargs, expected):
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(['A', 'A', 'B'], dtype=object),
        res_id=np.array([1, 2, 1], dtype=np.int32),
        atom_x=np.array([1.0, 2.0, 3.0], dtype=np.float32),
    )
    mask = np.array([True, False, True], dtype=bool)
    mask.flags.writeable = False
    np.testing.assert_array_equal(
        struc.filter(mask, invert=True, **kwargs).atom_x,
        np.array(expected, dtype=np.float32),
        strict=True,
    )

  def test_filter_does_not_mutate_callable_mask(self):
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(['A', 'A', 'B'], dtype=object),
        res_id=np.array([1, 2, 1], dtype=np.int32),
        atom_x=np.array([1.0, 2.0, 3.0], dtype=np.float32),
    )
    callable_mask = np.array([True, False, True], dtype=bool)
    np.testing.assert_array_equal(
        struc.filter(invert=True, atom_x=lambda _: callable_mask).atom_x,
        np.array([2.0], dtype=np.float32),
        strict=True,
    )
    np.testing.assert_array_equal(
        callable_mask, np.array([True, False, True]), strict=True
    )

  def test_filter_empty_chain_preserves_multimodel_shape(self):
    coords = np.array([[1.0], [2.0]], dtype=np.float32)
    struc = structure_lib.Structure(
        name='unresolved_chain',
        chains=structure_tables.Chains.from_defaults(
            id=np.array(['A', 'D'], dtype=object),
        ),
        residues=structure_tables.Residues.from_defaults(
            id=np.array([1, 1], dtype=np.int32),
            chain_key=np.array([0, 1], dtype=np.int64),
        ),
        atoms=structure_tables.Atoms.from_defaults(
            chain_key=np.array([0], dtype=np.int64),
            res_key=np.array([0], dtype=np.int64),
            x=coords,
            y=coords,
            z=coords,
        ),
        bonds=structure.Bonds.make_empty(),
    )
    filtered_chains = struc.filter(
        chain_id='D', cascade_delete=structure_lib.CascadeDelete.CHAINS
    )
    self.assertEmpty(filtered_chains.chains_table.id)
    self.assertEqual(filtered_chains.atom_x.shape, (2, 0))
    filtered_none = struc.filter(
        chain_id='D', cascade_delete=structure_lib.CascadeDelete.NONE
    )
    np.testing.assert_array_equal(
        filtered_none.chains_table.id,
        np.array(['A', 'D'], dtype=object),
        strict=True,
    )
    self.assertEqual(filtered_none.atom_x.shape, (2, 0))

  def test_filter_non_monotonic_chain_keys(self):
    coords = np.array([1.0, 2.0], dtype=np.float32)
    struc = structure_lib.Structure(
        name='non_monotonic',
        chains=structure_tables.Chains.from_defaults(
            key=np.array([10, 5], dtype=np.int64),
            id=np.array(['A', 'B'], dtype=object),
        ),
        residues=structure_tables.Residues.from_defaults(
            key=np.array([100, 50], dtype=np.int64),
            id=np.array([1, 1], dtype=np.int32),
            chain_key=np.array([10, 5], dtype=np.int64),
        ),
        atoms=structure_tables.Atoms.from_defaults(
            key=np.array([1000, 500], dtype=np.int64),
            chain_key=np.array([10, 5], dtype=np.int64),
            res_key=np.array([100, 50], dtype=np.int64),
            x=coords,
            y=coords,
            z=coords,
        ),
        bonds=structure.Bonds.make_empty(),
    )
    np.testing.assert_array_equal(
        struc.filter(chain_id='A').atom_x,
        np.array([1.0], dtype=np.float32),
        strict=True,
    )
    np.testing.assert_array_equal(
        struc.filter(chain_id='B').atom_x,
        np.array([2.0], dtype=np.float32),
        strict=True,
    )

  def test_filter_single_chain_owns_sliced_columns(self):
    coords = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    struc = structure.from_atom_arrays(
        name='owns_slice',
        chain_id=np.array(['A', 'B', 'C'], dtype=object),
        res_id=np.array([1, 1, 1], dtype=np.int32),
        atom_x=coords,
        atom_y=coords,
        atom_z=coords,
    )
    self.assertIsNone(struc.filter(chain_id='B').atom_x.base)

  def test_filter_unsorted_chain_keys(self):
    coords = np.array([1.0, 2.0], dtype=np.float32)
    base = structure.from_atom_arrays(
        name='unsorted',
        chain_id=np.array(['A', 'B'], dtype=object),
        res_id=np.array([1, 1], dtype=np.int32),
        atom_x=coords,
        atom_y=coords,
        atom_z=coords,
    )
    unsorted = base.copy_and_update(
        chains=base.chains_table.copy_and_update(
            key=np.array([1, 0], dtype=np.int64)
        ),
        residues=base.residues_table.copy_and_update(
            chain_key=np.array([1, 0], dtype=np.int64)
        ),
        atoms=base.atoms_table.copy_and_update(
            chain_key=np.array([1, 0], dtype=np.int64)
        ),
    )
    filtered = unsorted.filter(chain_id='B')
    np.testing.assert_array_equal(filtered.chain_id, ['B'])
    np.testing.assert_array_equal(
        filtered.atom_x, np.array([2.0], dtype=np.float32), strict=True
    )

  def test_filter_after_concat_preserves_chain_order(self):
    coords = np.array([1.0, 2.0], dtype=np.float32)
    s1 = structure.from_atom_arrays(
        name='s1',
        chain_id=np.array(['D', 'B'], dtype=object),
        res_id=np.array([1, 1], dtype=np.int32),
        atom_x=coords,
        atom_y=coords,
        atom_z=coords,
    )
    s2 = structure.from_atom_arrays(
        name='s2',
        chain_id=np.array(['C', 'A'], dtype=object),
        res_id=np.array([1, 1], dtype=np.int32),
        atom_x=coords + 10.0,
        atom_y=coords,
        atom_z=coords,
    )
    merged = structure.concat([s1, s2], assign_unique_chain_ids=False)
    filtered = merged.filter(chain_id=('B', 'A'))
    np.testing.assert_array_equal(filtered.chain_id, ['B', 'A'])
    np.testing.assert_array_equal(
        filtered.atom_x, np.array([2.0, 12.0], dtype=np.float32), strict=True
    )

  def test_filter_residues(self):
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(
            ['A', 'A', 'A', 'B', 'B', 'B', 'C', 'C'], dtype=object
        ),
        res_id=np.array([3, 3, 4, 2, 5, 5, 5, 5], dtype=np.int32),
        atom_name=np.array(
            ['C', 'CA', 'C', 'C', 'CA', 'C', 'C', 'CA'], dtype=object
        ),
    )
    res_mask = np.array([True, False, True, False, True])
    filtered = struc.filter_residues(res_mask)
    np.testing.assert_array_equal(filtered.chain_id, ['A', 'A', 'B', 'C', 'C'])
    np.testing.assert_array_equal(filtered.res_id, [3, 3, 2, 5, 5])
    np.testing.assert_array_equal(
        filtered.atom_name, ['C', 'CA', 'C', 'C', 'CA']
    )

  def test_filter_coords(self):
    struc = self.make_test_structure(
        name='test_structure',
        res_id=np.array([1, 1, 1, 1], dtype=np.int32),
        atom_name=np.array(['N', 'CA', 'C', 'O'], dtype=object),
        atom_x=np.array([0.0, 0.0, 0.0, 0.8], dtype=np.float32),
        atom_y=np.array([0.0, 1.0, 0.0, 0.8], dtype=np.float32),
        atom_z=np.array([0.0, 0.0, 2.0, 0.8], dtype=np.float32),
    )
    filtered = struc.filter_coords(lambda p: np.linalg.norm(p) > 1)  # pyrefly: ignore[bad-argument-type]
    np.testing.assert_array_equal(filtered.atom_name, ['C', 'O'])
    np.testing.assert_array_equal(
        filtered.atom_x, np.array([0.0, 0.8], dtype=np.float32)
    )
    np.testing.assert_array_equal(
        filtered.atom_y, np.array([0.0, 0.8], dtype=np.float32)
    )
    np.testing.assert_array_equal(
        filtered.atom_z, np.array([2.0, 0.8], dtype=np.float32)
    )

  def test_filter_polymers_to_single_atom_per_res(self):
    chain_to_atoms_dict = {
        mmcif_names.PROTEIN_CHAIN: 'CA',
        mmcif_names.NON_POLYMER_CHAIN: 'CA',
        mmcif_names.DNA_CHAIN: "C1'",
    }
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(
            ['A', 'A', 'A', 'B', 'B', 'B', 'C', 'C', 'D', 'D'], dtype=object
        ),
        res_id=np.array([3, 3, 4, 2, 5, 5, 5, 5, 7, 8], dtype=np.int32),
        chain_type=np.array(
            [mmcif_names.DNA_CHAIN] * 3
            + [mmcif_names.PROTEIN_CHAIN] * 3
            + [mmcif_names.NON_POLYMER_CHAIN] * 2
            + [mmcif_names.RNA_CHAIN] * 2,
            dtype=object,
        ),
        atom_name=np.array(
            ['O', "C1'", 'P', 'C', 'CA', 'C', 'C', 'CA', "C2'", 'OP1'],
            dtype=object,
        ),
    )
    filtered = struc.filter_polymers_to_single_atom_per_res(chain_to_atoms_dict)
    np.testing.assert_array_equal(filtered.chain_id, ['A', 'B', 'C', 'D', 'D'])
    np.testing.assert_array_equal(filtered.res_id, [3, 5, 5, 7, 8])
    np.testing.assert_array_equal(
        filtered.atom_name, ["C1'", 'CA', 'CA', "C2'", 'OP1']
    )

  def test_filter_polymers_to_single_atom_per_res_repeated_chain_atoms(self):
    chain_to_atoms_dict = {
        mmcif_names.DNA_CHAIN: "C4'",
        mmcif_names.RNA_CHAIN: "C2'",
    }
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(
            ['A', 'A', 'A', 'B', 'B', 'B', 'C', 'C', 'D', 'D'], dtype=object
        ),
        res_id=np.array([3, 3, 4, 2, 5, 5, 5, 5, 7, 8], dtype=np.int32),
        chain_type=np.array(
            [mmcif_names.DNA_CHAIN] * 3
            + [mmcif_names.PROTEIN_CHAIN] * 3
            + [mmcif_names.NON_POLYMER_CHAIN] * 2
            + [mmcif_names.RNA_CHAIN] * 2,
            dtype=object,
        ),
        atom_name=np.array(
            ['O', "C4'", "C2'", 'C', 'CA', 'C', 'C', 'CA', "C2'", "C4'"],
            dtype=object,
        ),
    )
    filtered = struc.filter_polymers_to_single_atom_per_res(chain_to_atoms_dict)
    np.testing.assert_array_equal(
        filtered.chain_id, ['A', 'B', 'B', 'B', 'C', 'C', 'D']
    )
    np.testing.assert_array_equal(filtered.res_id, [3, 2, 5, 5, 5, 5, 7])
    np.testing.assert_array_equal(
        filtered.atom_name, ["C4'", 'C', 'CA', 'C', 'C', 'CA', "C2'"]
    )

  @parameterized.named_parameters(
      dict(
          testcase_name='empty_mapping_cascades_unresolved_chain',
          mapping={},
          exp_chains=['A', 'B', 'C'],
          exp_atoms=['N', 'CA', 'C', 'CA', "C2'", 'OP1'],
      ),
      dict(
          testcase_name='shared_target_atom_with_sparse_chain_keys',
          mapping={
              mmcif_names.PROTEIN_CHAIN: 'CA',
              mmcif_names.NON_POLYMER_CHAIN: 'CA',
          },
          exp_chains=['A', 'B', 'C'],
          exp_atoms=['CA', 'CA', "C2'", 'OP1'],
      ),
  )
  def test_filter_polymers_to_single_atom_per_res_edge_cases(
      self, mapping, exp_chains, exp_atoms
  ):
    struc = structure_lib.Structure(
        name='polymer_edge_cases',
        chains=structure_tables.Chains.from_defaults(
            key=np.array([100, 5000, 9000, 12000], dtype=np.int64),
            id=np.array(['A', 'B', 'C', 'D'], dtype=object),
            type=np.array(
                [
                    mmcif_names.PROTEIN_CHAIN,
                    mmcif_names.NON_POLYMER_CHAIN,
                    mmcif_names.RNA_CHAIN,
                    mmcif_names.PROTEIN_CHAIN,
                ],
                dtype=object,
            ),
        ),
        residues=structure_tables.Residues.from_defaults(
            id=np.array([1, 1, 1, 1], dtype=np.int32),
            chain_key=np.array([100, 5000, 9000, 12000], dtype=np.int64),
        ),
        atoms=structure_tables.Atoms.from_defaults(
            chain_key=np.array(
                [100, 100, 5000, 5000, 9000, 9000], dtype=np.int64
            ),
            res_key=np.array([0, 0, 1, 1, 2, 2], dtype=np.int64),
            name=np.array(['N', 'CA', 'C', 'CA', "C2'", 'OP1'], dtype=object),
        ),
        bonds=structure.Bonds.make_empty(),
    )
    filtered = struc.filter_polymers_to_single_atom_per_res(mapping)
    np.testing.assert_array_equal(
        filtered.chains_table.id,
        np.array(exp_chains, dtype=object),
        strict=True,
    )
    np.testing.assert_array_equal(
        filtered.atom_name, np.array(exp_atoms, dtype=object), strict=True
    )

  def test_atom_ids(self):
    expected = [
        ('A', '1', None, 'C'), ('A', '2', None, 'C'), ('A', '2', None, 'O'),
        ('B', '1', None, 'C'), ('B', '3', None, 'C'), ('B', '4', None, 'C'),
    ]  # pyformat: disable

    struc = self.make_test_structure(
        chain_id=np.array(['A', 'A', 'A', 'B', 'B', 'B'], dtype=object),
        res_id=np.array([1, 2, 2, 1, 3, 4], dtype=np.int32),
        atom_name=np.array(['C', 'C', 'O', 'C', 'C', 'C'], dtype=object),
    )
    self.assertEqual(struc.atom_ids, expected)

  @parameterized.named_parameters(
      [('auth_scheme', True), ('no_auth_scheme', False)]
  )
  def test_label_asym_id_to_entity_id(self, use_auth_scheme):
    struc = structure.Structure(
        chains=structure.Chains.from_defaults(
            key=np.array([0, 1, 2], dtype=np.int64),
            id=np.array(['A', 'B', 'C'], dtype=object),
            entity_id=np.array(['1', '3', '2'], dtype=object),
        ),
        residues=structure.Residues.from_defaults(
            key=np.array([0, 1, 2, 3, 4], dtype=np.int64),
            chain_key=np.array([0, 0, 1, 1, 2], dtype=np.int64),
            id=np.array([1, 2, 1, 3, 4], dtype=np.int32),
        ),
        atoms=structure.Atoms.from_defaults(
            chain_key=np.array([0, 0, 0, 1, 1, 2], dtype=np.int64),
            res_key=np.array([0, 1, 1, 2, 3, 4], dtype=np.int64),
            name=np.array([0, 1, 1, 2, 3, 4], dtype=np.int64),
        ),
        bonds=structure.Bonds.make_empty(),
    )
    if use_auth_scheme:
      expected = {'A': '1', 'B': '3', 'C': '2'}  # '3' and '2' reversed.
    else:
      expected = {'A': '1', 'B': '2', 'C': '3'}  # '3' and '2' in order.
      struc = struc.reset_author_naming_scheme()
    self.assertEqual(struc.label_asym_id_to_entity_id, expected)

  def test_order_and_drop_atoms_to_match(self):
    struc = self.make_test_structure(
        chain_id=np.array(['A', 'A', 'A', 'A', 'B', 'B', 'B'], dtype=object),
        res_id=np.array([1, 1, 1, 2, 10, 10, 10], dtype=np.int32),
        atom_name=np.array(['C', 'CA', 'N', 'C', 'C', 'CA', 'N'], dtype=object),
        atom_b_factor=np.array([0, 1, 2, 3, 4, 5, 6], dtype=np.float32),
    )
    other = self.make_test_structure(
        chain_id=np.array(['B', 'B', 'B', 'A', 'A', 'A'], dtype=object),
        res_id=np.array([10, 10, 10, 1, 1, 1], dtype=np.int32),
        atom_name=np.array(['N', 'C', 'CA', 'N', 'CA', 'C'], dtype=object),
        atom_b_factor=np.array([0, 1, 2, 3, 4, 5], dtype=np.float32),
    )

    with self.subTest('with_missing_residue'):
      aligned = struc.order_and_drop_atoms_to_match(other)
      np.testing.assert_array_equal(aligned.chains, other.chains)
      self.assertSequenceEqual(
          [r['res_id'] for r in aligned.iter_residues()],
          [r['res_id'] for r in other.iter_residues()],
      )
      np.testing.assert_array_equal(aligned.res_id, other.res_id)
      np.testing.assert_array_equal(aligned.atom_name, other.atom_name)
      np.testing.assert_array_equal(
          aligned.atom_b_factor, [6.0, 4.0, 5.0, 2.0, 1.0, 0.0]
      )
      self.assertDictEqual(aligned.all_residues, other.all_residues)

    with self.subTest('missing_atom_error'):
      with self.assertRaises(structure.MissingAtomError):
        other.order_and_drop_atoms_to_match(struc)

    with self.subTest('without_missing_residue'):
      struc = struc.filter(res_id=(1, 10))  # Remove res 2 missing from other.
      other = other.order_and_drop_atoms_to_match(struc)
      np.testing.assert_array_equal(other.chains, struc.chains)
      self.assertSequenceEqual(
          [r['res_id'] for r in other.iter_residues()],
          [r['res_id'] for r in struc.iter_residues()],
      )
      np.testing.assert_array_equal(other.res_id, struc.res_id)
      np.testing.assert_array_equal(other.atom_name, struc.atom_name)
      np.testing.assert_array_equal(
          other.atom_b_factor, [5.0, 4.0, 3.0, 1.0, 2.0, 0.0]
      )

  def test_order_and_drop_atoms_to_match_with_missing_atoms(self):
    struc = self.make_test_structure(
        chain_id=np.array(['A', 'A', 'A'], dtype=object),
        res_id=np.array([1, 2, 3], dtype=np.int32),
        atom_name=np.array(['C', 'C', 'C'], dtype=object),
    )
    other = self.make_test_structure(
        chain_id=np.array(['A', 'A', 'A'], dtype=object),
        res_id=np.array([2, 3, 4], dtype=np.int32),
        atom_name=np.array(['C', 'C', 'C'], dtype=object),
    )
    with self.assertRaises(structure.MissingAtomError):
      _ = struc.order_and_drop_atoms_to_match(other)
    struc = struc.order_and_drop_atoms_to_match(other, allow_missing_atoms=True)
    np.testing.assert_array_equal(struc.chains, other.chains)
    self.assertSequenceEqual(
        [r['res_id'] for r in struc.iter_residues()], [2, 3]
    )
    np.testing.assert_array_equal(struc.res_id, [2, 3])

  def test_order_and_drop_atoms_to_match_with_empty_residues(self):
    struc = self.make_test_structure(
        chain_id=np.array(['A', 'A'], dtype=object),
        res_id=np.array([1, 2], dtype=np.int32),
        res_name=np.array(['ARG', 'ARG'], dtype=object),
        atom_name=np.array(['C', 'C'], dtype=object),
        all_residues={'A': [('ARG', 1), ('ARG', 2), ('ARG', 3)]},
    )
    other = self.make_test_structure(
        chain_id=np.array(['A', 'A'], dtype=object),
        res_id=np.array([1, 2], dtype=np.int32),
        res_name=np.array(['ARG', 'ARG'], dtype=object),
        atom_name=np.array(['C', 'C'], dtype=object),
        all_residues={'A': [('ARG', 1), ('ARG', 2), ('ARG', 3)]},
    )
    struc = struc.order_and_drop_atoms_to_match(other)
    np.testing.assert_array_equal(struc.chains, other.chains)
    self.assertEqual(
        struc.all_residues, {'A': [('ARG', 1), ('ARG', 2), ('ARG', 3)]}
    )

  def test_order_and_drop_atoms_to_match_identical_reshuffled(self):
    mmcif_str = test_utils.read_test_cif_file('2pv7')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False).filter(
        chain_id='A'
    )
    struc = structure.concat([struc] * 3)
    reordered = struc.rename_chain_ids({'A': 'C', 'B': 'B', 'C': 'A'})

    struc = struc.filter(struc.res_id > 50)

    restored = reordered.order_and_drop_atoms_to_match(
        struc, allow_missing_atoms=True
    )
    np.testing.assert_array_equal(restored.chains, struc.chains)
    self.assertSequenceEqual(
        [r['res_id'] for r in restored.iter_residues()],
        [r['res_id'] for r in struc.iter_residues()],
    )
    self.assertStructuresEqual(
        struc,
        restored,
        ignore_fields=('author_naming_scheme', 'bioassembly_data'),
    )

  @parameterized.parameters(
      ('CA', ('A', 'P', 'B', 'X', 'CA')),
      ('H', ('A', 'H', 'P', 'B', 'X')),
      ('M', ('A', 'P', 'B', 'M', 'X')),
      ('V', ('A', 'P', 'V', 'B', 'X')),
  )
  def test_order_and_drop_atoms_to_match_real(self, keep_ligand, exp_chains):
    gt_cif = test_utils.read_test_cif_file('6xq5')
    gt_struc = structure.from_mmcif(gt_cif, include_bonds=False)
    gt_struc = gt_struc.generate_bioassembly()

    pred_cif = test_utils.read_test_cif_file('6xq5_pred')
    pred_struc = structure.from_mmcif(pred_cif, include_bonds=False)

    gt_struc = gt_struc.filter(chain_id=('A', 'P', 'B', 'X', keep_ligand))
    pred_struc = pred_struc.filter(chain_id=('A', 'P', 'B', 'X', keep_ligand))
    matched_struc = pred_struc.order_and_drop_atoms_to_match(gt_struc)
    matched_cif = pred_struc.to_mmcif_dict()

    self.assertEqual(gt_struc.num_atoms, 6307)
    self.assertEqual(pred_struc.num_atoms, 6307)
    self.assertEqual(matched_struc.num_atoms, 6307)
    self.assertSequenceEqual(
        matched_struc.chain_id[matched_struc.chain_boundaries].tolist(),
        exp_chains,
    )
    self.assertEqual(matched_struc.chains, exp_chains)
    self.assertEqual(sum(matched_struc.is_ligand_mask), 19)
    self.assertEqual(sum(matched_struc.is_protein_mask), 6307 - 19)
    self.assertSequenceEqual(
        matched_cif['_struct_asym.id'], ('A', 'B', 'P', 'X', keep_ligand)
    )

  def test_copy_and_update(self):
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(['A', 'B'], dtype=object),
        res_id=np.array([1, 1], dtype=np.int32),
    )
    updated = struc.copy_and_update_globals(name='updated')
    self.assertEqual(struc.name, 'test_structure')
    self.assertEqual(updated.name, 'updated')
    np.testing.assert_array_equal(struc.chain_id, updated.chain_id)

  def test_copy_and_update_with_no_changes(self):
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(['A'], dtype=object),
        res_id=np.array([1], dtype=np.int32),
    )
    with self.assertRaisesRegex(
        ValueError, 'Unnecessary call to copy_and_update with no changes..*'
    ):
      struc.copy_and_update()

  def test_copy_and_update_atom_fields(self):
    mmcif_str = test_utils.read_test_cif_file('2pv7')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)
    updated = struc.copy_and_update_atoms(atom_x=struc.atom_x + 42.0)
    self.assertStructuresEqual(struc, updated, ignore_fields=['atom_x'])

  def test_copy_and_update_residue_fields(self):
    mmcif_str = test_utils.read_test_cif_file('2pv7')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)
    updated = struc.copy_and_update(
        residues=struc.residues_table.copy_and_update(
            name=np.array(['ARG'] * len(struc.residues_table), dtype=object),
        )
    )
    self.assertStructuresEqual(
        struc, updated, ignore_fields=['res_name', 'all_residues']
    )

  def test_copy_and_update_chain_fields(self):
    mmcif_str = test_utils.read_test_cif_file('2pv7')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)
    old_chains = struc.chains_table
    updated = struc.copy_and_update(
        chains=old_chains.copy_and_update(
            type=np.array(
                [mmcif_names.PROTEIN_CHAIN] * len(old_chains), dtype=object
            )
        )
    )
    self.assertStructuresEqual(struc, updated, ignore_fields=['chain_type'])

  def test_copy_and_update_coords(self):
    mmcif_str = test_utils.read_test_cif_file('2pv7')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)

    updated_coords = struc.coords + 42.0
    updated = struc.copy_and_update_coords(coords=updated_coords)
    np.testing.assert_array_almost_equal(
        updated.coords - struc.coords, [[42.0] * 3] * struc.num_atoms, decimal=5
    )

  def test_copy_and_update_coords_shape_validation(self):
    mmcif_str = test_utils.read_test_cif_file('2pv7')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)

    with self.subTest('leading_dims_ok'):
      struc.copy_and_update_coords(
          coords=np.zeros((10, struc.num_atoms, 3), dtype=np.float32)
      )

    with self.subTest('bad_num_atoms'):
      with self.assertRaises(ValueError):
        struc.copy_and_update_coords(coords=np.zeros((10, 3)))

    with self.subTest('bad_num_coords'):
      with self.assertRaises(ValueError):
        struc.copy_and_update_coords(coords=np.zeros((struc.num_atoms, 2)))

  def test_copy_and_update_globals(self):
    struc = self.make_test_structure(
        name='test_structure',
        release_date=datetime.date(2024, 2, 2),
        chain_id=np.array(['A', 'B'], dtype=object),
        res_id=np.array([1, 1], dtype=np.int32),
    )
    updated = struc.copy_and_update_globals(
        name='updated',
        release_date=datetime.date(3000, 1, 1),
        resolution=42.0,
        structure_method='test',
    )
    self.assertNotEqual(struc.name, updated.name)
    self.assertNotEqual(struc.release_date, updated.release_date)
    self.assertNotEqual(struc.resolution, updated.resolution)
    self.assertNotEqual(struc.structure_method, updated.structure_method)
    self.assertEqual(updated.name, 'updated')
    self.assertEqual(updated.release_date, datetime.date(3000, 1, 1))
    self.assertEqual(updated.resolution, 42.0)
    self.assertEqual(updated.structure_method, 'test')

  def test_copy_and_update_atoms(self):
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(['A', 'B'], dtype=object),
        res_id=np.array([1, 1], dtype=np.int32),
        atom_name=np.array(['C', 'C'], dtype=object),
    )

    updated = struc.copy_and_update_atoms(
        atom_name=np.array(['CA', 'N'], dtype=object),
        atom_element=np.array(['C', 'N'], dtype=object),
        atom_x=np.array([1.0, 2.0], dtype=np.float32),
        atom_y=np.array([3.0, 4.0], dtype=np.float32),
        atom_z=np.array([5.0, 6.0], dtype=np.float32),
        atom_b_factor=np.array([7.0, 8.0], dtype=np.float32),
        atom_occupancy=np.array([4.0, 7.0], dtype=np.float32),
    )
    np.testing.assert_array_equal(updated.atom_name, ['CA', 'N'])
    np.testing.assert_array_equal(updated.atom_element, ['C', 'N'])
    np.testing.assert_array_equal(updated.atom_x, [1.0, 2.0])
    np.testing.assert_array_equal(updated.atom_y, [3.0, 4.0])
    np.testing.assert_array_equal(updated.atom_z, [5.0, 6.0])
    np.testing.assert_array_equal(updated.atom_b_factor, [7.0, 8.0])
    np.testing.assert_array_equal(updated.atom_occupancy, [4.0, 7.0])

  def test_copy_and_update_residues(self):
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(['A', 'B'], dtype=object),
        res_id=np.array([1, 1], dtype=np.int32),
        res_name=np.array(['CA', 'CA'], dtype=object),
    )

    updated = struc.copy_and_update_residues(
        res_id=np.array([2, 3], dtype=np.int32),
        res_name=np.array(['O', 'N'], dtype=object),
    )
    np.testing.assert_array_equal(updated.res_id, [2, 3])
    np.testing.assert_array_equal(updated.res_name, ['O', 'N'])

  @parameterized.parameters(
      dict(
          chain_id=['A', 'A', 'A'],
          res_id=[1, 1, 1],
          atom_name=['C', 'CA', 'N'],
          atom_b_factor=[5.0],
          expected_atom_b_factor=[5.0, 5.0, 5.0],
      ),
      dict(
          chain_id=['A', 'A', 'A'],
          res_id=[1, 2, 3],
          atom_name=['C', 'C', 'C'],
          atom_b_factor=[5.0, 6.0, 7.0],
          expected_atom_b_factor=[5.0, 6.0, 7.0],
      ),
      dict(
          chain_id=['B', 'B', 'A', 'A'],
          res_id=[1, 1, 1, 1],
          atom_name=['C', 'CA', 'C', 'CA'],
          atom_b_factor=[5.0, 6.0],
          expected_atom_b_factor=[5.0, 5.0, 6.0, 6.0],
      ),
      dict(
          chain_id=['A', 'A', 'A', 'A', 'A', 'B', 'B', 'B'],
          res_id=[1, 1, 1, 2, 2, 2, 3, 3],
          atom_name=['C', 'CA', 'N', 'C', 'CA', 'C', 'CA', 'C'],
          atom_b_factor=[5.0, 6.0, 7.0, 8.0],
          expected_atom_b_factor=[5.0, 5.0, 5.0, 6.0, 6.0, 7.0, 8.0, 8.0],
      ),
  )
  def test_copy_and_update_from_res_arrays(
      self, chain_id, res_id, atom_name, atom_b_factor, expected_atom_b_factor
  ):
    assert len(chain_id) == len(res_id)
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(chain_id, dtype=object),
        res_id=np.array(res_id, dtype=np.int32),
        atom_name=np.array(atom_name, dtype=object),
    )
    updated = struc.copy_and_update_from_res_arrays(
        atom_b_factor=np.array(atom_b_factor, dtype=np.float32)
    )
    self.assertStructuresEqual(struc, updated, ignore_fields=['atom_b_factor'])
    np.testing.assert_array_equal(expected_atom_b_factor, updated.atom_b_factor)

  def test_copy_and_update_from_res_arrays_with_unresolved(self):
    unfiltered_struc = self.make_test_structure(
        res_id=np.array([1, 2, 3, 4, 5, 6, 7, 8], dtype=np.int32),
        chain_id=np.array(['A'] * 4 + ['B'] * 4, dtype=object),
    )
    struc = unfiltered_struc.filter(
        np.array([False, False, True, True, False, False, True, True]),
        cascade_delete=structure.CascadeDelete.NONE,
    )
    with self.subTest('handle_unresolved'):
      updated_struc = struc.copy_and_update_from_res_arrays(
          include_unresolved=True,
          atom_b_factor=np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]),
      )
      np.testing.assert_array_equal(
          updated_struc.atom_b_factor, [3.0, 4.0, 7.0, 8.0]
      )
    with self.subTest('raise_on_incorrect_num_residues'):
      with self.assertRaisesRegex(
          ValueError,
          'atom_b_factor array of length 8 does not match num_residues=4.*',
      ):
        updated_struc.copy_and_update_from_res_arrays(
            include_unresolved=False,
            atom_b_factor=np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]),
        )

  @parameterized.parameters([True, False])
  def test_add_bonds(self, with_preexisting_bonds):
    if with_preexisting_bonds:
      orig_bonds = structure.Bonds(
          key=np.array([1, 2], dtype=np.int64),
          type=np.array(['type1', 'type2'], dtype=object),
          role=np.array(['role1', 'role2'], dtype=object),
          from_atom_key=np.array([1, 0], dtype=np.int64),
          dest_atom_key=np.array([0, 2], dtype=np.int64),
      )
      atom_key = np.array([0, 1, 2, 3], dtype=np.int64)
    else:
      orig_bonds = None
      atom_key = None
    struc = self.make_test_structure(
        name='test',
        chain_id=np.array(['A', 'A', 'A', 'B'], dtype=object),
        res_id=np.array([1, 1, 3, 2], dtype=np.int32),
        atom_name=np.array(['N', 'C', 'C', 'N'], dtype=object),
        atom_key=atom_key,
        bonds=orig_bonds,
        shuffle_atom_keys=False,
    )
    with self.assertRaises(KeyError):
      struc.add_bonds([(('A', 1, 'N'), ('A', 2, 'C'))])
    struc_with_added_bonds = struc.add_bonds(
        [
            (('A', 1, 'N'), ('A', 3, 'C')),
            (('B', 2, 'N'), ('A', 1, 'C')),
        ],
        bond_type='metalc',
    )
    self.assertStructuresEqual(
        struc, struc_with_added_bonds, ignore_fields=['bonds']
    )
    if with_preexisting_bonds:
      np.testing.assert_array_equal(
          struc_with_added_bonds.bonds.from_atom_key, [1, 0, 0, 3]
      )
      np.testing.assert_array_equal(
          struc_with_added_bonds.bonds.dest_atom_key, [0, 2, 2, 1]
      )
      np.testing.assert_array_equal(
          struc_with_added_bonds.bonds.type,
          np.array(['type1', 'type2', 'metalc', 'metalc'], dtype=object),
      )
    else:
      np.testing.assert_array_equal(
          struc_with_added_bonds.bonds.from_atom_key, [0, 3]
      )
      np.testing.assert_array_equal(
          struc_with_added_bonds.bonds.dest_atom_key, [2, 1]
      )
      np.testing.assert_array_equal(
          struc_with_added_bonds.bonds.type,
          np.array(['metalc', 'metalc'], dtype=object),
      )

  def test_coords(self):
    struc = self.make_test_structure(
        res_id=np.array([1, 1, 1], dtype=np.int32),
        atom_x=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float32),
        atom_y=np.array(
            [[11.0, 12.0, 13.0], [14.0, 15.0, 16.0]], dtype=np.float32
        ),
        atom_z=np.array(
            [[21.0, 22.0, 23.0], [24.0, 25.0, 26.0]], dtype=np.float32
        ),
    )
    np.testing.assert_array_equal(
        struc.coords,
        np.array([
            [[1.0, 11.0, 21.0], [2.0, 12.0, 22.0], [3.0, 13.0, 23.0]],
            [[4.0, 14.0, 24.0], [5.0, 15.0, 25.0], [6.0, 16.0, 26.0]],
        ]),
    )

  @parameterized.parameters([True, False])
  def test_chain_single_letter_sequence(self, include_missing_residues):
    res = ['ALA', 'MSE', 'UNK', 'A', 'UNK', 'CCC', 'ALA', 'VWW', 'BR1', 'BR2']
    struc = self.make_test_structure(
        chain_id=np.array(
            ['A', 'A', 'A', 'B', 'B', 'B', 'C', 'D', 'E', 'E'], dtype=object
        ),
        chain_type=np.array(
            [
                *['polypeptide(L)'] * 3,
                *['polyribonucleotide'] * 3,
                *['unknown', 'non-polymer'],
                *['branched'] * 2,
            ],
            dtype=object,
        ),
        res_name=np.array(res, dtype=object),
        res_id=np.array([1, 2, 3, 1, 2, 3, 1, 1, 1, 2], dtype=np.int32),
        all_residues={
            'A': [('ALA', 1), ('MSE', 2), ('UNK', 3), ('GLY', 4)],
            'B': [('A', 1), ('UNK', 2), ('CCC', 3), ('U', 4)],
            'C': [('ALA', 1)],
            'D': [('VWW', 1)],
            'E': [('BR1', 1), ('BR2', 2)],
        },
    )
    if include_missing_residues:
      expected = {'A': 'AMXG', 'B': 'AXCU', 'C': 'X', 'D': 'X', 'E': 'XX'}
    else:
      expected = {'A': 'AMX', 'B': 'AXC', 'C': 'X', 'D': 'X', 'E': 'XX'}
    self.assertDictEqual(
        struc.chain_single_letter_sequence(include_missing_residues), expected
    )

    struc = struc.filter(chain_id='A', res_id=1)
    if include_missing_residues:
      expected = {'A': 'AMXG'}
    else:
      expected = {'A': 'A'}
    self.assertDictEqual(
        struc.chain_single_letter_sequence(include_missing_residues), expected
    )

  @parameterized.parameters(
      ('2pv7', {'A': 'A', 'B': 'B'}, True, True, True, True),
      ('2ppb_trimmed', {'G': 'A', 'H': 'B', 'I': 'C', 'X': 'D', 'Y': 'E',
                        'Z': 'F', 'A': 'G', 'B': 'H', 'C': 'I', 'D': 'J',
                        'E': 'K', 'K': 'L', 'L': 'M', 'M': 'N', 'N': 'O',
                        'O': 'P'},
       True, True, True, True),
      ('2ppb_trimmed', {'A': 'G', 'B': 'H', 'C': 'I', 'D': 'J', 'E': 'K',
                        'K': 'L', 'L': 'M', 'M': 'N', 'N': 'O', 'O': 'P'},
       True, False, False, True),
      ('2ppb_trimmed', {'H': 'B', 'Y': 'E'}, False, True, False, True),
      ('2ppb_trimmed', {'G': 'A', 'I': 'C', 'X': 'D', 'Z': 'F'},
       False, False, True, True),
      ('2ppb_trimmed', {}, False, False, False, True),
      ('3b9f', {'L': 'A', 'H': 'B', 'I': 'C'}, True, True, True, True),
      # 2kql has only 1 polypeptide(D).
      ('2kql_trimmed', {}, True, True, True, False),
      ('2kql_trimmed', {'A': 'A'}, True, True, True, True),
  )  # pyformat: disable
  def test_polymer_auth_asym_id_to_label_asym_id(
      self,
      pdb_id,
      exp_auth_to_label,
      protein,
      rna,
      dna,
      other,
  ):
    mmcif_str = test_utils.read_test_cif_file(pdb_id)
    struc = structure.from_mmcif(
        mmcif_str, include_other=True, include_bonds=False
    )
    auth_to_label = struc.polymer_auth_asym_id_to_label_asym_id(
        protein=protein, rna=rna, dna=dna, other=other
    )
    self.assertDictEqual(auth_to_label, exp_auth_to_label)

  @parameterized.parameters(
      ('2pv7', True, ('A', 'B')),
      ('2pv7', False, ('A', 'B')),
      ('2ppb_trimmed', True, set('ABCDEGHIKLMNOXYZ')),
      ('2ppb_trimmed', False, set('ABCDEGHIKLMNOXYZ')),
      ('3b9f', True, ('L', 'H', 'I')),
      ('3b9f', False, ('L', 'H', 'I')),
  )
  def test_polymer_author_chain_single_letter_sequence(
      self, pdb_id, include_missing_residues, exp_auth_label_ids
  ):
    mmcif_str = test_utils.read_test_cif_file(pdb_id)
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)

    auth_to_seq = struc.polymer_author_chain_single_letter_sequence(
        include_missing_residues=include_missing_residues
    )
    self.assertCountEqual(auth_to_seq.keys(), exp_auth_label_ids)

    label_to_seq = struc.chain_single_letter_sequence(include_missing_residues)
    self.assertContainsSubset(auth_to_seq.values(), label_to_seq.values())

  @parameterized.product(
      include_missing_residues=[True, False],
      fix_non_standard_polymer_res=[True, False],
  )
  def test_chain_res_name_sequence(
      self, include_missing_residues, fix_non_standard_polymer_res
  ):
    struc = self.make_test_structure(
        chain_id=np.array([
            'A', 'A', 'A', 'A', 'B', 'B', 'B', 'C', 'D', 'E', 'E', 'F', 'F',
            'F', 'G', 'G', 'G'
        ], dtype=object),
        chain_type=np.array((
            *['polypeptide(L)'] * 4,
            *['polyribonucleotide'] * 3,
            *['unknown'] * 1,
            *['non-polymer'] * 1,
            *['branched'] * 2,
            *['polypeptide(L)'] * 3,
            *['polydeoxyribonucleotide'] * 3,
        ), dtype=object),
        res_name=np.array([
            'ALA', 'MET', 'MET', 'UNK', 'A', 'G', 'C', 'ALA', 'VWW', 'BR1',
            'BR2', 'MSE', 'ALA', 'XXXYYY', 'VVVVWWW', 'XTR', 'DT',
        ], dtype=object),
        res_id=np.array(
            [1, 2, 2, 3, 1, 2, 3, 1, 1, 1, 2, 1, 2, 3, 1, 2, 3],
            dtype=np.int32),
        atom_name=np.array([*['C'] * 2, 'CA', *['C'] * 14], dtype=object),
        all_residues={
            'A': [('ALA', 1), ('MET', 2), ('UNK', 3), ('GLY', 4)],
            'B': [('A', 1), ('G', 2), ('C', 3), ('U', 4)],
            'C': [('ALA', 1)],
            'D': [('VWW', 1)],
            'E': [('BR1', 1), ('BR2', 2)],
            'F': [('MSE', 1), ('ALA', 2), ('XXXYYY', 3)],
            'G': [('VVVVWWW', 1), ('XTR', 2), ('DT', 3)],
        },
    )  # pyformat: disable
    if include_missing_residues:
      expected = {
          'A': ('ALA', 'MET', 'UNK', 'GLY'),
          'B': ('A', 'G', 'C', 'U'),
          'C': ('ALA',),
          'D': ('VWW',),
          'E': ('BR1', 'BR2'),
      }
    else:
      expected = {
          'A': ('ALA', 'MET', 'UNK'),
          'B': ('A', 'G', 'C'),
          'C': ('ALA',),
          'D': ('VWW',),
          'E': ('BR1', 'BR2'),
      }
    if fix_non_standard_polymer_res:
      expected.update({'F': ('MET', 'ALA', 'UNK'), 'G': ('N', 'DT', 'DT')})  # pyrefly: ignore[no-matching-overload]
    else:
      expected.update({  # pyrefly: ignore[no-matching-overload]
          'F': ('MSE', 'ALA', 'XXXYYY'),
          'G': ('VVVVWWW', 'XTR', 'DT'),
      })

    self.assertDictEqual(
        struc.chain_res_name_sequence(
            include_missing_residues=include_missing_residues,
            fix_non_standard_polymer_res=fix_non_standard_polymer_res,
        ),
        expected,
    )

    struc = struc.filter(chain_id='A', res_id=1)
    if include_missing_residues:
      expected = {'A': ('ALA', 'MET', 'UNK', 'GLY')}
    else:
      expected = {'A': ('ALA',)}

    self.assertDictEqual(
        struc.chain_res_name_sequence(
            include_missing_residues=include_missing_residues,
            fix_non_standard_polymer_res=False,
        ),
        expected,
    )

  @parameterized.parameters([True, False])
  def test_fix_non_standard_polymer_res(self, without_all_residues):
    struc = self.make_test_structure(
        name='test_structure',
        all_residues={
            'A': [('ALA', 1), ('MSE', 2)],
            'B': [('MET', 1), ('ACE', 2), ('UNK', 3)],
            'C': [('S04', 1)],
        },
        chain_id=np.array(['A', 'A', 'B', 'B', 'B', 'C'], dtype=object),
        chain_type=np.array(
            ['polypeptide(L)'] * 5 + ['non-polymer'], dtype=object
        ),
        res_id=np.array([1, 2, 1, 2, 3, 1], dtype=np.int32),
        res_name=np.array(
            ['ALA', 'MSE', 'MET', 'ACE', 'UNK', 'S04'], dtype=object
        ),
    )
    if without_all_residues:
      struc = struc.copy_and_update(residues=struc.present_residues)
    fixed = struc.fix_non_standard_polymer_res()
    self.assertEqual(fixed.name, struc.name)
    np.testing.assert_array_equal(fixed.chain_id, struc.chain_id)
    np.testing.assert_array_equal(fixed.chain_type, struc.chain_type)
    np.testing.assert_array_equal(fixed.res_id, struc.res_id)
    np.testing.assert_array_equal(
        fixed.res_name, ['ALA', 'MET', 'MET', 'UNK', 'UNK', 'S04']
    )
    self.assertDictEqual(
        fixed.all_residues,
        {
            'A': [('ALA', 1), ('MET', 2)],
            'B': [('MET', 1), ('UNK', 2), ('UNK', 3)],
            'C': [('S04', 1)],
        },
    )

  def test_fix_and_drop_non_standard_polymer_res_real(self):
    mmcif_str = test_utils.read_test_cif_file('1al4')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)
    # Make the last/first chain proteins to catch off-by-one bugs at the edge.
    struc = struc.filter_to_entity_type(protein=True)
    struc = struc.fix_non_standard_polymer_res()
    struc = struc.drop_non_standard_atoms(
        ccd=chemical_components.Ccd(), drop_unk=True, drop_non_ccd=False
    )

    self.assertEqual(struc.num_atoms, 526)
    # fmt: off
    self.assertEqual(
        struc.group_by_residue.res_name.tolist(),
        [
            'VAL', 'GLY', 'ALA', 'LEU', 'ALA', 'VAL', 'VAL', 'VAL', 'TRP',
            'LEU', 'TRP', 'LEU', 'TRP', 'LEU', 'TRP', 'VAL', 'GLY', 'ALA',
            'LEU', 'ALA', 'VAL', 'VAL', 'VAL', 'TRP', 'LEU', 'TRP', 'LEU',
            'TRP', 'LEU', 'TRP',
        ],
    )
    # fmt: on

  def test_slice_leading_dims(self):
    struc = self.make_test_structure(
        name='test_structure',
        res_id=np.array([0, 1, 2], dtype=np.int32),
        atom_x=np.array(
            [[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]], dtype=np.float32
        ),  # shape: (1, 2, 3)
        atom_y=np.array(
            [[[11.0, 12.0, 13.0], [14.0, 15.0, 16.0]]], dtype=np.float32
        ),
        atom_z=np.array(
            [[[21.0, 22.0, 23.0], [24.0, 25.0, 26.0]]], dtype=np.float32
        ),
        atom_b_factor=np.array([0.5, 1.0, 1.5], dtype=np.float32),
    )
    sliced = struc.slice_leading_dims[0, 1]
    np.testing.assert_array_equal(
        sliced.coords,
        np.array([
            [4.0, 14.0, 24.0],
            [5.0, 15.0, 25.0],
            [6.0, 16.0, 26.0],
        ]),
    )
    np.testing.assert_array_equal(
        sliced.atom_b_factor, np.array([0.5, 1.0, 1.5])
    )
    with self.assertRaises(ValueError):
      _ = struc.slice_leading_dims[0, 0, 0]

  def test_multidim_b_factor(self):
    struc = self.make_test_structure(
        name='test_structure',
        res_id=np.array([0, 1, 2], dtype=np.int32),
        # shape: (1, 2, 3)
        atom_x=np.array([[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]], dtype=np.float32),
        atom_y=np.array(
            [[[11.0, 12.0, 13.0], [14.0, 15.0, 16.0]]], dtype=np.float32
        ),
        atom_z=np.array(
            [[[21.0, 22.0, 23.0], [24.0, 25.0, 26.0]]], dtype=np.float32
        ),
        atom_b_factor=np.array(
            [[[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]], dtype=np.float32
        ),
        atom_occupancy=np.array(
            [[[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]], dtype=np.float32
        ),
    )
    sliced = struc.slice_leading_dims[0, 1]
    np.testing.assert_array_equal(
        sliced.atom_b_factor, np.array([0.4, 0.5, 0.6], dtype=np.float32)
    )
    atom = list(struc.iter_atoms())[0]
    expected_atom = {
        'atom_x': np.array([[1.0, 4.0]], dtype=np.float32),
        'atom_y': np.array([[11.0, 14.0]], dtype=np.float32),
        'atom_z': np.array([[21.0, 24.0]], dtype=np.float32),
        'atom_b_factor': np.array([[0.1, 0.4]], dtype=np.float32),
        'atom_occupancy': np.array([[0.1, 0.4]], dtype=np.float32),
        'res_id': 0,
    }
    for field in expected_atom:
      if isinstance(expected_atom[field], np.ndarray):
        np.testing.assert_array_equal(atom[field], expected_atom[field])
      else:
        self.assertEqual(atom[field], expected_atom[field])
    with self.assertRaises(ValueError):
      _ = struc.filter(atom_b_factor=lambda b: np.any(b < 0.5))  # pyrefly: ignore[bad-argument-type]

  @parameterized.parameters([-2, 0, 1])
  def test_unstack(self, axis):
    struc = self.make_test_structure(
        name='test_structure',
        res_id=np.array([0, 1, 2], dtype=np.int32),
        atom_x=np.array(
            [[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]], dtype=np.float32
        ),  # shape: (1, 2, 3)
        atom_y=np.array(
            [[[11.0, 12.0, 13.0], [14.0, 15.0, 16.0]]], dtype=np.float32
        ),
        atom_z=np.array(
            [[[21.0, 22.0, 23.0], [24.0, 25.0, 26.0]]], dtype=np.float32
        ),
        atom_b_factor=np.array([0.5, 1.0, 1.5], dtype=np.float32),
    )
    unstacked = struc.unstack(axis=axis)
    self.assertLen(unstacked, struc.atom_x.shape[axis])
    expected_shape = list(struc.atom_x.shape)
    del expected_shape[axis]
    for i, unstacked_struc in enumerate(unstacked):
      with self.subTest(f'struc_{i}'):
        self.assertSequenceEqual(unstacked_struc.atom_x.shape, expected_shape)
        np.testing.assert_array_equal(
            unstacked_struc.atom_x, np.take(struc.atom_x, i, axis=axis)
        )
        np.testing.assert_array_equal(
            unstacked_struc.atom_y, np.take(struc.atom_y, i, axis=axis)
        )
        np.testing.assert_array_equal(
            unstacked_struc.atom_z, np.take(struc.atom_z, i, axis=axis)
        )

  @parameterized.parameters([-5, -1, 2, 3])
  def test_unstack_invalid_axis(self, axis):
    struc = self.make_test_structure(
        name='test_structure',
        res_id=np.array([0, 1, 2], dtype=np.int32),
        atom_x=np.array(
            [[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]], dtype=np.float32
        ),  # shape: (1, 2, 3)
        atom_y=np.array(
            [[[11.0, 12.0, 13.0], [14.0, 15.0, 16.0]]], dtype=np.float32
        ),
        atom_z=np.array(
            [[[21.0, 22.0, 23.0], [24.0, 25.0, 26.0]]], dtype=np.float32
        ),
        atom_b_factor=np.array([0.5, 1.0, 1.5], dtype=np.float32),
    )
    with self.assertRaises(IndexError):
      struc.unstack(axis=axis)

  def test_reorder_chains(self):
    all_residues = {
        'A': [('ARG', 1), ('ARG', 2), ('ALA', 3), ('UNK', 4)],
        'B': [('GLU', 1), ('THR', 2), ('GLY', 3)],
        'C': [('CYS', 1)],
    }
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(['A', 'A', 'A', 'B', 'B', 'C'], dtype=object),
        chain_type=np.array([mmcif_names.PROTEIN_CHAIN] * 6, dtype=object),
        res_id=np.array([1, 2, 3, 1, 3, 1], dtype=np.int32),
        res_name=np.array(
            ['ARG', 'ARG', 'ALA', 'GLU', 'GLY', 'CYS'], dtype=object
        ),
        all_residues=all_residues,
        atom_x=np.arange(6, dtype=np.float32),
    )
    new_chain_order = ('B', 'A', 'C')
    new_res_name_order = ('GLU', 'GLY', 'ARG', 'ARG', 'ALA', 'CYS')
    new_res_id_order = (1, 3, 1, 2, 3, 1)
    reordered = struc.reorder_chains(new_order=new_chain_order)

    with self.subTest('chains_reordered'):
      self.assertSequenceEqual(reordered.chains, new_chain_order)
      np.testing.assert_array_equal(
          reordered.chain_id,
          np.array(['B', 'B', 'A', 'A', 'A', 'C'], dtype=object),
      )

    with self.subTest('residues_reordered'):
      np.testing.assert_array_equal(
          reordered.res_id, np.array(new_res_id_order, dtype=np.int64)
      )
      np.testing.assert_array_equal(
          reordered.res_name,
          np.array(new_res_name_order, dtype=object),
      )
      self.assertSequenceEqual(
          [r['res_id'] for r in reordered.iter_residues()], new_res_id_order
      )
      self.assertSequenceEqual(
          [r['res_name'] for r in reordered.iter_residues()], new_res_name_order
      )
      self.assertDictEqual(
          reordered.chain_single_letter_sequence(
              include_missing_residues=False
          ),
          {'A': 'RRA', 'B': 'EG', 'C': 'C'},
      )
      self.assertDictEqual(
          reordered.chain_single_letter_sequence(include_missing_residues=True),
          {'A': 'RRAX', 'B': 'ETG', 'C': 'C'},
      )

    with self.subTest('all_residues_unchanged'):
      self.assertDictEqual(reordered.all_residues, all_residues)

    with self.subTest('atoms_reordered'):
      np.testing.assert_array_equal(
          reordered.atom_x,
          np.array([3.0, 4.0, 0.0, 1.0, 2.0, 5.0], dtype=np.float32),
      )

  @parameterized.named_parameters(
      ('different_count', ('A', 'B', 'C')),
      ('not_unique', ('A', 'B', 'C', 'D', 'F', 'F')),
      ('different_set', ('A', 'B', 'C', 'D', 'E', 'G')),
  )
  def test_reorder_chains_errors(self, new_chain_order):
    struc = structure.from_mmcif(
        test_utils.read_test_cif_file('2pv7'), include_bonds=False
    )

    with self.assertRaises(ValueError):
      struc.reorder_chains(new_order=new_chain_order)

  def test_auth_cached_properties(self):
    struc = structure.Structure(
        chains=structure.Chains.from_defaults(
            key=np.arange(3, dtype=np.int64),
            id=np.array(['A', 'B', 'C'], dtype=object),
            auth_asym_id=np.array(['X', 'Y', 'A'], dtype=object),
            entity_id=np.array(['3', '1', '2'], dtype=object),
            entity_desc=np.array(
                ['THIRD', 'FIRST', 'SECOND'],
                dtype=object,
            ),
        ),
        residues=structure.Residues.from_defaults(
            key=np.array([0, 1, 2, 3, 4, 5], dtype=np.int64),
            chain_key=np.array([0, 0, 0, 1, 1, 2], dtype=np.int64),
            id=np.array([1, 2, 3, 1, 3, 1], dtype=np.int32),
            auth_seq_id=np.array(
                ['5', '4', '10', '-2', '3', '4'], dtype=object
            ),
            insertion_code=np.array(
                ['A', 'B', '?', 'B', 'A', 'A'], dtype=object
            ),
        ),
        atoms=structure.Atoms.from_defaults(
            chain_key=np.array([0, 0, 0, 1, 1, 2], np.int64),
            res_key=np.array([0, 1, 2, 3, 4, 5], dtype=np.int64),
        ),
        bonds=structure.Bonds.make_empty(),
    )
    np.testing.assert_array_equal(
        struc.chain_auth_asym_id, ('X', 'X', 'X', 'Y', 'Y', 'A')
    )
    np.testing.assert_array_equal(
        struc.chain_entity_id, ('3', '3', '3', '1', '1', '2')
    )
    np.testing.assert_array_equal(
        struc.chain_entity_desc,
        ('THIRD', 'THIRD', 'THIRD', 'FIRST', 'FIRST', 'SECOND'),
    )
    np.testing.assert_array_equal(
        struc.res_auth_seq_id, ('5', '4', '10', '-2', '3', '4')
    )
    np.testing.assert_array_equal(
        struc.res_insertion_code, ('A', 'B', '?', 'B', 'A', 'A')
    )

  def test_rename_auth_asym_ids(self):
    struc = structure.from_mmcif(
        test_utils.read_test_cif_file('2pv7'), include_bonds=False
    )
    renamed = struc.rename_auth_asym_ids({'A': 'B', 'B': 'A'})
    with self.subTest('renames_correctly'):
      self.assertEqual(renamed.author_naming_scheme.auth_asym_id['A'], 'B')
      self.assertEqual(renamed.author_naming_scheme.auth_asym_id['B'], 'A')

    new_struc = structure.from_mmcif(renamed.to_mmcif(), include_bonds=False)
    with self.subTest('round_trip_to_mmcif'):
      self.assertEqual(new_struc.author_naming_scheme.auth_asym_id['A'], 'B')
      self.assertEqual(new_struc.author_naming_scheme.auth_asym_id['B'], 'A')

  def test_rename_chain_ids_all(self):
    mmcif_str = test_utils.read_test_cif_file('7a9y')
    # Before the rename - note that F and G are water molecules, which are
    # filtered out when loading the structure, but kept in the metadata fields.
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)
    self.assertSequenceEqual(struc.chains, ('A', 'B', 'C', 'D', 'E'))
    self.assertEqual(
        struc.author_naming_scheme.auth_asym_id,
        {
            **struc.author_naming_scheme.auth_asym_id,
            'A': 'BBB',
            'B': 'AAA',
            'C': 'BBB',
            'D': 'BBB',
            'E': 'AAA',
        },
    )
    self.assertCountEqual(
        struc.author_naming_scheme.auth_seq_id.keys(), {'A', 'B', 'C', 'D', 'E'}
    )
    self.assertCountEqual(
        struc.author_naming_scheme.insertion_code.keys(),
        {'A', 'B', 'C', 'D', 'E'},
    )
    self.assertEqual(
        struc.author_naming_scheme.entity_id,
        {
            **struc.author_naming_scheme.entity_id,
            'A': '1',
            'B': '1',
            'C': '2',
            'D': '3',
            'E': '4',
        },
    )
    self.assertEqual(
        struc.author_naming_scheme.entity_desc,
        {
            **struc.author_naming_scheme.entity_desc,
            '1': 'Cellular retinoic acid-binding protein 1',
            '2': 'MYRISTIC ACID',
            '3': 'GLYCEROL',
            '4': 'N-TRIDECANOIC ACID',
        },
    )
    old_bd = struc.bioassembly_data
    self.assertIsNotNone(old_bd)
    self.assertSequenceEqual(
        old_bd.to_mmcif_dict()['_pdbx_struct_assembly_gen.asym_id_list'],
        ['B,E,G', 'A,C,D,F', 'A,C,D,F', 'B,E,G'],
    )

    # After the rename.
    new_struc = struc.rename_chain_ids({'A': 'X', 'B': 'ZZ'})
    self.assertSequenceEqual(new_struc.chains, ('X', 'ZZ', 'C', 'D', 'E'))
    new_bd = new_struc.bioassembly_data
    self.assertIsNotNone(new_bd)
    self.assertSequenceEqual(
        new_bd.to_mmcif_dict()['_pdbx_struct_assembly_gen.asym_id_list'],
        ['ZZ,E', 'X,C,D', 'X,C,D', 'ZZ,E'],
    )
    # Check that the renaming has no impact on bioassembly generation.
    np.testing.assert_array_equal(
        struc.generate_bioassembly().coords,
        new_struc.generate_bioassembly().coords,
    )

  def test_rename_chain_ids_chains(self):
    mmcif_str = test_utils.read_test_cif_file('6gnx')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)

    with self.subTest(msg='Raises for many to one maps'):
      with self.assertRaises(ValueError):
        struc.rename_chain_ids({'A': 'A', 'B': 'A'})

    with self.subTest(msg='Renames correctly'):
      renamed = struc.rename_chain_ids(
          {'A': 'A2', 'C': 'C2', 'B': 'B2', 'D': 'D2'}
      )
      self.assertSequenceEqual(renamed.chains, ('A2', 'B2', 'C2', 'D2'))

    with self.subTest(msg='Noop'):
      renamed = struc.rename_chain_ids(new_id_by_old_id={})
      self.assertSequenceEqual(renamed.chains, ('A', 'B', 'C', 'D'))

  def test_rename_chain_ids_after_filter(self):
    mmcif_str = test_utils.read_test_cif_file('6gnx')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)
    self.assertSequenceEqual(struc.chains, ('A', 'B', 'C', 'D'))
    struc = struc.filter(chain_id='A')
    self.assertSequenceEqual(struc.chains, 'A')
    # Rename A to a chain ID B that existed in the structure before filtering.
    struc = struc.rename_chain_ids({'A': 'B'})
    self.assertSequenceEqual(struc.chains, 'B')

  def test_rename_chain_ids_author_names(self):
    mmcif_str = test_utils.read_test_cif_file('2pv7')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False).filter(
        chain_id='B'
    )
    renamed = struc.rename_chain_ids({'B': 'Z'})

    original_cif_dict = struc.to_mmcif_dict()
    renamed_cif_dict = renamed.to_mmcif_dict()
    # The first present residue in chain B has author ID 92 and internal ID 13.
    # Since we remove author names, the author ID should equal to 13.
    with self.subTest('author_residue_name'):
      self.assertEqual(original_cif_dict['_atom_site.auth_seq_id'][0], '92')
      self.assertEqual(original_cif_dict['_atom_site.label_seq_id'][0], '13')
      self.assertEqual(renamed_cif_dict['_atom_site.auth_seq_id'][0], '13')
      self.assertEqual(renamed_cif_dict['_atom_site.label_seq_id'][0], '13')

    with self.subTest('author_chain_name'):
      self.assertEqual(original_cif_dict['_atom_site.auth_asym_id'][0], 'B')
      self.assertEqual(original_cif_dict['_atom_site.label_asym_id'][0], 'B')
      self.assertEqual(renamed_cif_dict['_atom_site.auth_asym_id'][0], 'Z')
      self.assertEqual(renamed_cif_dict['_atom_site.label_asym_id'][0], 'Z')

  def test_rename_chain_ids_error(self):
    mmcif_str = test_utils.read_test_cif_file('6tug_trimmed')
    struc = structure.from_mmcif(
        mmcif_str, include_bonds=False
    ).filter_to_entity_type(protein=True)
    self.assertSequenceEqual(
        struc.chains, ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')
    )
    with self.assertRaises(ValueError):
      struc.rename_chain_ids({'A': 'G'})

  def test_reset_author_naming_scheme(self):
    # Test that reset_author_naming_scheme undoes author_naming_scheme.
    example_author_naming_scheme = structure.AuthorNamingScheme(
        auth_asym_id={'A': 'X', 'B': 'Y'},
        auth_seq_id={
            'A': {1: '5', 2: '4', 3: '10'},
            'B': {1: '-2', 2: '3', 3: '3'},
        },
        insertion_code={
            'A': {1: None, 2: None, 3: None},
            'B': {1: None, 2: 'A', 3: 'B'},
        },
        entity_id={'A': '2', 'B': '1'},
        entity_desc={'2': 'ACID A', '1': 'ALKALI B'},
    )

    no_scheme = structure.Structure(
        chains=structure.Chains.from_defaults(
            key=np.arange(2, dtype=np.int64),
            id=np.array(['A', 'B'], dtype=object),
        ),
        residues=structure.Residues.from_defaults(
            key=np.arange(6, dtype=np.int64),
            chain_key=np.array([0, 0, 0, 1, 1, 1], dtype=np.int64),
            id=np.array([1, 2, 3, 1, 2, 3], dtype=np.int32),
            name=np.array(
                ['ALA', 'ILE', 'GLY', 'ASN', 'PRO', 'MET'], dtype=object
            ),
        ),
        atoms=structure.Atoms.from_defaults(
            chain_key=np.array([0, 0, 0, 1, 1], dtype=np.int64),
            res_key=np.array([0, 1, 2, 3, 4], dtype=np.int64),
            name=np.array(['CA', 'CA', 'CA', 'CA', 'CA'], dtype=object),
        ),
        bonds=structure.Bonds.make_empty(),
    )
    with_scheme = no_scheme.copy_and_update(
        chains=no_scheme.chains_table.copy_and_update(
            auth_asym_id=np.array(['X', 'Y'], dtype=object),
            entity_id=np.array(['2', '1'], dtype=object),
            entity_desc=np.array(['ACID A', 'ALKALI B'], dtype=object),
        ),
        residues=no_scheme.residues_table.copy_and_update(
            auth_seq_id=np.array(
                ['5', '4', '10', '-2', '3', '3'], dtype=object
            ),
            insertion_code=np.array(
                [None, None, None, None, 'A', 'B'], dtype=object
            ),
        ),
    )

    print(with_scheme.author_naming_scheme)
    print(example_author_naming_scheme)

    self.assertEqual(
        with_scheme.author_naming_scheme, example_author_naming_scheme
    )
    scheme_removed = with_scheme.reset_author_naming_scheme()

    self.assertEqual(
        no_scheme.author_naming_scheme, scheme_removed.author_naming_scheme
    )

    self.assertStructuresEqual(no_scheme, scheme_removed)

  @parameterized.parameters(
      '2pv7', '2ppb_trimmed', '1wib', '3b9f', '5pti', '402d'
  )
  def test_reset_author_naming_scheme_real(self, pdb_id):
    mmcif_str = test_utils.read_test_cif_file(pdb_id)
    raw_struc = structure.from_mmcif(mmcif_str, include_bonds=False)
    reset_struc = raw_struc.reset_author_naming_scheme()
    np.testing.assert_array_equal(
        reset_struc.chains_table.auth_asym_id,
        reset_struc.chains_table.id,
    )
    np.testing.assert_array_equal(
        reset_struc.chains_table.entity_id,
        (reset_struc.chains_table.key + 1).astype(str).astype(object),
    )
    np.testing.assert_array_equal(
        reset_struc.chains_table.entity_desc,
        np.full(len(reset_struc.chains_table), fill_value='.', dtype=object),
    )
    np.testing.assert_array_equal(
        reset_struc.residues_table.auth_seq_id,
        reset_struc.residues_table.id.astype(str).astype(object),
    )
    np.testing.assert_array_equal(
        reset_struc.residues_table.insertion_code,
        np.full(len(reset_struc.residues_table), fill_value='?', dtype=object),
    )

  def test_rename_chains_to_match(self):
    mmcif_str = test_utils.read_test_cif_file('2pv7')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)

    with self.subTest(name='bad'):
      # Raises when chains don't match.
      with self.assertRaises(ValueError):
        struc_without_chain_a = struc.filter_out(chain_id='A')
        struc_without_chain_a.rename_chains_to_match(struc)

    with self.subTest(name='good'):
      # Standard use case.
      renamed_struc = struc.rename_chain_ids({'A': 'X'})

      # 2pv7 contains MSE residues which can be mapped to standard MET.
      for fuzzy_match in (True, False):
        restored_struc = renamed_struc.rename_chains_to_match(
            struc, fuzzy_match_non_standard_res=fuzzy_match
        )
        self.assertNotEqual(renamed_struc.chains, struc.chains)
        self.assertSameElements(restored_struc.chains, struc.chains)
        self.assertDictEqual(
            restored_struc.chain_single_letter_sequence(),
            struc.chain_single_letter_sequence(),
        )

        self.assertDictEqual(
            restored_struc.chain_res_name_sequence(),
            struc.chain_res_name_sequence(),
        )

  def test_rename_chains_to_match_fuzzy(self):
    # Check fuzzy matching is required for chains that differ only by one
    # standard/non-standard pair.
    mmcif_str = test_utils.read_test_cif_file('2pv7')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)

    standardized_struc = struc.rename_res_name(
        {'MSE': 'MET'}, fail_if_not_found=False
    )

    with self.subTest("Doesn't match with fuzzy flag off."):
      with self.assertRaises(ValueError):
        standardized_struc.rename_chains_to_match(
            struc, fuzzy_match_non_standard_res=False
        )

      with self.assertRaises(ValueError):
        struc.rename_chains_to_match(
            standardized_struc, fuzzy_match_non_standard_res=False
        )

    with self.subTest('Does match with fuzzy flag on.'):
      struc.rename_chains_to_match(
          standardized_struc, fuzzy_match_non_standard_res=True
      )
      standardized_struc.rename_chains_to_match(
          struc, fuzzy_match_non_standard_res=True
      )

  @parameterized.parameters(True, False)
  def test_rename_res_name(self, use_all_residues):
    struc = self.make_test_structure(
        all_residues=(
            {  # pylint: disable=g-long-ternary
                'A': [('ALA', 1), ('GLY', 2), ('MSE', 3)],
                'B': [('GLY', 1), ('MSE', 2)],
            }
            if use_all_residues
            else None
        ),
        chain_id=np.array(['A', 'A', 'A', 'B'], dtype=object),
        res_name=np.array(['ALA', 'GLY', 'MSE', 'GLY'], dtype=object),
        res_id=np.array([1, 2, 3, 1], dtype=np.int32),
        chemical_components_data=(
            structure.get_data_for_ccd_components(
                ccd=chemical_components.Ccd(),
                chemical_component_ids=['ALA', 'GLY', 'MSE'],
            )
        ),
    )
    with self.assertRaises(ValueError):
      struc.rename_res_name(
          {'GLY': 'MSE', 'MSE': 'HIS', 'CYS': 'PRO'}, fail_if_not_found=True
      )
    renamed = struc.rename_res_name(
        {'GLY': 'MSE', 'MSE': 'HIS', 'CYS': 'PRO'}, fail_if_not_found=False
    )
    np.testing.assert_array_equal(
        renamed.res_name, ['ALA', 'MSE', 'HIS', 'MSE']
    )
    renamed_data = renamed.chemical_components_data
    self.assertIsNotNone(renamed_data)
    np.testing.assert_array_equal(
        sorted(renamed_data.chem_comp.keys()),
        ['ALA', 'HIS', 'MSE'],
    )

    if use_all_residues:
      self.assertDictEqual(
          renamed.all_residues,
          {
              'A': [('ALA', 1), ('MSE', 2), ('HIS', 3)],
              'B': [('MSE', 1), ('HIS', 2)],
          },
      )
    else:
      self.assertDictEqual(
          renamed.all_residues,
          {
              'A': [('ALA', 1), ('MSE', 2), ('HIS', 3)],
              'B': [('MSE', 1)],
          },
      )

  def test_remap_res_id(self):
    struc = self.make_test_structure(
        all_residues=({
            'A': [('ALA', 1), ('GLY', 2), ('MSE', 3)],
            'B': [('GLY', 1), ('MSE', 2)],
        }),
        chain_id=np.array(['A', 'A', 'A', 'B'], dtype=object),
        res_name=np.array(['ALA', 'GLY', 'MSE', 'GLY'], dtype=object),
        res_id=np.array([1, 2, 3, 1], dtype=np.int32),
    )
    new_struc = struc.remap_res_id(
        {'A': {1: 1, 2: 3, 3: 4}, 'B': {1: 10, 2: 11}}
    )
    np.testing.assert_array_equal(
        new_struc.res_id, np.array([1, 3, 4, 10], dtype=np.int32), strict=True
    )
    np.testing.assert_array_equal(
        new_struc.residues_table.id,
        np.array([1, 3, 4, 10, 11], dtype=np.int32),
        strict=True,
    )

  def test_remap_res_id_missing_key(self):
    struc = self.make_test_structure(
        chain_id=np.array(['A'], dtype=object),
        res_name=np.array(['ALA'], dtype=object),
        res_id=np.array([1], dtype=np.int32),
    )
    with self.assertRaises(KeyError):
      struc.remap_res_id({'A': {2: 1}})

  def test_remap_res_id_not_uniquely_mapped(self):
    struc = self.make_test_structure(
        chain_id=np.array(['A', 'A'], dtype=object),
        res_name=np.array(['ALA', 'GLY'], dtype=object),
        res_id=np.array([1, 2], dtype=np.int32),
    )
    with self.assertRaises(ValueError):
      struc.remap_res_id({'A': {1: 2, 2: 2}})

  def test_generate_bioassemblies_atom_key_numbering(self):
    mmcif_str = test_utils.read_test_cif_file('4bxt')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)
    assembly_1 = struc.generate_bioassembly(assembly_id='1')
    assembly_2 = struc.generate_bioassembly(assembly_id='2')

    with self.subTest('assembly 1 chains E, F, G, H'):
      offset = struc.filter(chain_id=('A', 'B', 'C', 'D')).num_atoms
      np.testing.assert_array_equal(
          assembly_1.atom_key,
          np.arange(offset, offset + assembly_1.num_atoms, dtype=np.int64),
      )

    with self.subTest('assembly 2 chains A, B, C, D'):
      np.testing.assert_array_equal(
          assembly_2.atom_key, np.arange(assembly_2.num_atoms, dtype=np.int64)
      )

  # Info on bioassembly test cases:
  # 2bfu: A capsid structure with (mainly) one-to-many bioassemblies.
  # 7a9y: A structure with one-to-one bioassemblies where one transform is
  #    applied to one set of chains and another transform is applied to a
  #    different set of chains.
  # 6tug: A structure with bioassemblies that simply crop out chains from
  #   the asymmetric unit.
  @parameterized.parameters('2bfu', '7a9y', '6tug_trimmed')
  def test_generate_bioassemblies(self, pdb_id):
    # This tests that we can generate structures for all the available
    # bioassemblies and that they contain the expected oligomeric count.
    mmcif_str = test_utils.read_test_cif_file(pdb_id)
    struc = structure.from_mmcif(mmcif_str, include_bonds=True)
    struc_bioassembly_data = struc.bioassembly_data
    self.assertIsNotNone(struc_bioassembly_data)
    for assembly_id in struc_bioassembly_data.assembly_ids:
      with self.subTest(assembly_id=assembly_id):
        assembly_info = struc_bioassembly_data.get_assembly_info(assembly_id)
        expected_num_protein_chains = int(assembly_info['oligomeric_count'])
        assembly_struc = struc.generate_bioassembly(assembly_id)
        num_protein_chains = assembly_struc.filter_to_entity_type(
            protein=True
        ).num_chains
        self.assertIsNone(assembly_struc.bioassembly_data)
        self.assertEqual(assembly_struc.name, struc.name)
        self.assertEqual(assembly_struc.release_date, struc.release_date)
        self.assertEqual(assembly_struc.resolution, struc.resolution)
        self.assertEqual(
            assembly_struc.structure_method, struc.structure_method
        )
        self.assertEqual(
            assembly_struc.chemical_components_data,
            struc.chemical_components_data,
        )
        self.assertEqual(num_protein_chains, expected_num_protein_chains)

        orig_seqs = struc.chain_res_name_sequence()
        assembly_seqs = assembly_struc.chain_res_name_sequence()
        orig_ans = struc.author_naming_scheme
        assembly_ans = assembly_struc.author_naming_scheme
        assert orig_ans is not None
        assert assembly_ans is not None
        for orig_chain_id, orig_entity_id in orig_ans.entity_id.items():
          for chain_id, entity_id in assembly_ans.entity_id.items():
            if entity_id == orig_entity_id:
              self.assertEqual(
                  assembly_seqs[chain_id], orig_seqs[orig_chain_id]
              )
          for chain_id, auth_chain_id in assembly_ans.auth_asym_id.items():
            self.assertEqual(chain_id, auth_chain_id)

        for (
            assembly_chain_id,
            chain_id,
        ) in struc_bioassembly_data.asym_id_by_assembly_chain_id(
            assembly_id
        ).items():
          struc_chain = struc.filter(chain_id=chain_id)
          assembly_chain = assembly_struc.filter(chain_id=assembly_chain_id)
          np.testing.assert_array_equal(
              struc_chain.chain_type, assembly_chain.chain_type
          )
          np.testing.assert_array_equal(
              struc_chain.res_id, assembly_chain.res_id
          )
          self.assertBondsEqual(
              struc_chain.bonds,
              assembly_chain.bonds,
              struc_chain.atom_key,
              assembly_chain.atom_key,
          )

  @parameterized.parameters('2bfu', '7a9y', '6tug_trimmed', '1wib')
  def test_generate_default_bioassembly(self, pdb_id):
    mmcif_str = test_utils.read_test_cif_file(pdb_id)
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)

    if struc.bioassembly_data:
      assembly_id = struc.bioassembly_data.get_default_assembly_id()
      assembly_info = struc.bioassembly_data.get_assembly_info(assembly_id)
      expected_num_protein_chains = int(assembly_info['oligomeric_count'])
    else:
      expected_num_protein_chains = struc.filter_to_entity_type(
          protein=True
      ).num_chains

    assembly_struc = struc.generate_bioassembly()
    num_protein_chains = assembly_struc.filter_to_entity_type(
        protein=True
    ).num_chains
    self.assertIsNone(assembly_struc.bioassembly_data)
    self.assertEqual(num_protein_chains, expected_num_protein_chains)

  def test_generate_bioassemblies_2bfu(self):
    # For 2bfu we check that the bioassembly structure contains a chain who's
    # first CA atom has coordinates equal to a known value that isn't present
    # in the raw structure.
    mmcif_str = test_utils.read_test_cif_file('2bfu')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)
    target = np.array([7.839, -109.065, 23.082], dtype=np.float32)
    assembly_struc = struc.generate_bioassembly('1')
    found = False
    for chain_id in assembly_struc.chains:
      chain_ca_struc = assembly_struc.filter(chain_id=chain_id, atom_name='CA')
      first_ca_coords = chain_ca_struc.coords[0]
      if np.linalg.norm(first_ca_coords - target) < 0.01:
        found = True
        break
    self.assertTrue(found)

  def test_generate_bioassemblies_6tug(self):
    # For 6tug we check that the bioassembly only contains the chains in the
    # bioassembly as the bioassembly is just a crop of the asymmetric unit.
    mmcif_str = test_utils.read_test_cif_file('6tug_trimmed')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)
    assembly_struc = struc.generate_bioassembly('1')
    self.assertSequenceEqual(
        assembly_struc.chains, ('A', 'B', 'I', 'J', 'K', 'L', 'M', 'N', 'O')
    )

  def test_generate_bioassemblies_consistent_mmcif_ordering(self):
    # For 2fid we check that the chains are serialised in an order consistent
    # with the atom site ordering after bioassembly expansion for both
    # _pdbx_nonpoly_scheme.pdb_strand_id and _pdbx_poly_seq_scheme.seq_id.
    mmcif_str = test_utils.read_test_cif_file('2fid')
    struc = structure.from_mmcif(
        mmcif_str, include_water=True, include_other=True, include_bonds=False
    )
    assembly_struc = struc.generate_bioassembly()
    mmcif_dict = assembly_struc.to_mmcif_dict()

    with self.subTest('pdbx_poly_seq_scheme'):
      np.testing.assert_array_equal(
          mmcif_dict['_pdbx_poly_seq_scheme.asym_id'],
          np.array(['A'] * 76 + ['F'] * 76 + ['B'] * 70, dtype=object),
      )
    with self.subTest('pdbx_nonpoly_scheme'):
      np.testing.assert_array_equal(
          mmcif_dict['_pdbx_nonpoly_scheme.asym_id'],
          np.array(['D', 'C', 'G', 'E'], dtype=object),
      )

    with self.subTest('chains'):
      self.assertSequenceEqual(
          assembly_struc.chains, ('A', 'D', 'F', 'B', 'C', 'G', 'E')
      )

  def test_generate_bioassemblies_5o61(self):
    # 5o61 has a duplicate author residue IDs in chain QA.
    mmcif_str = test_utils.read_test_cif_file('5o61_bad_auth_asym_ids')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)
    self.assertEqual(struc.num_atoms, 4)
    self.assertEqual(struc.num_residues(count_unresolved=True), 4)
    self.assertEqual(struc.num_residues(count_unresolved=False), 4)

  def test_generate_bioassemblies_partial(self):
    # 6tug assembly 1 applies only to chains A,B,I,J,K,L,M,N,O,KA,LA.
    # We should produce an empty structure if none of these chains are present.
    mmcif_str = test_utils.read_test_cif_file('6tug_trimmed')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)
    partial_struc = struc.filter(chain_id=('C', 'D'))
    assembly_struc = partial_struc.generate_bioassembly('1')
    self.assertSequenceEqual(assembly_struc.chains, [])
    self.assertEqual(assembly_struc.num_atoms, 0)

  @parameterized.parameters(
      ('100d', True, True, True, False, False, []),
      ('100d', True, True, True, False, True, [2]),
      ('2pv7', True, True, True, False, False, [2]),
      ('2pv7', True, True, True, True, False, [2, 2, 2]),
      ('2ppb_trimmed', True, False, False, False, False, [4, 2, 2, 2]),
      ('2ppb_trimmed', True, True, True, False, False, [4, 2, 2, 2, 2, 2, 2]),
      (
          '2ppb_trimmed',
          True,
          True,
          True,
          True,
          False,
          [4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2],
      ),
      ('3b9f', True, False, False, False, False, [1, 1, 1]),
      ('3b9f', True, True, True, False, False, [1, 1, 1]),
      ('3b9f', True, True, True, True, False, [5, 3, 1, 1, 1, 1, 1]),
  )
  def test_get_stoichiometry(
      self,
      pdb_id,
      keep_protein,
      keep_dna,
      keep_rna,
      keep_ligands,
      keep_dna_rna_hybrid,
      expected_stoichiometry,
  ):
    mmcif_str = test_utils.read_test_cif_file(pdb_id)
    struc = structure.from_mmcif(
        mmcif_str, include_bonds=False
    ).filter_to_entity_type(
        protein=keep_protein,
        dna=keep_dna,
        rna=keep_rna,
        dna_rna_hybrid=keep_dna_rna_hybrid,
        ligand=keep_ligands,
    )
    self.assertSequenceEqual(struc.get_stoichiometry(), expected_stoichiometry)

  @parameterized.parameters(
      # fmt: off
      ('2pv7', ('A', 'B'), tuple(), tuple(), ('C', 'D', 'E', 'F'), ('G', 'H')),
      ('2ppb_trimmed', ('G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P'),
       ('B', 'E'), ('A', 'C', 'D', 'F'),
       ('Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'AA', 'BA'),
       ('CA', 'DA', 'EA', 'FA', 'GA', 'HA', 'IA', 'JA', 'KA', 'LA', 'MA', 'NA',
        'OA', 'PA', 'QA', 'RA')),
      ('3b9f', ('A', 'B', 'C'), tuple(), tuple(),
       ('D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M'), ('N', 'O', 'P')),
      # fmt: on
  )
  def test_filter_to_entity_type(
      self,
      pdb_id,
      protein_chains,
      rna_chains,
      dna_chains,
      ligand_chains,
      water_chains,
  ):
    mmcif_str = test_utils.read_test_cif_file(pdb_id)
    struc = structure.from_mmcif(
        mmcif_str, include_water=True, include_bonds=False
    )
    filtered = struc.filter_to_entity_type(protein=True)
    self.assertSequenceEqual(filtered.chains, protein_chains)
    filtered = struc.filter_to_entity_type(rna=True)
    self.assertSequenceEqual(filtered.chains, rna_chains)
    filtered = struc.filter_to_entity_type(dna=True)
    self.assertSequenceEqual(filtered.chains, dna_chains)
    filtered = struc.filter_to_entity_type(ligand=True)
    self.assertSequenceEqual(filtered.chains, ligand_chains)
    filtered = struc.filter_to_entity_type(water=True)
    self.assertSequenceEqual(filtered.chains, water_chains)

    # Try some combinations:
    filtered = struc.filter_to_entity_type(protein=True, rna=True, dna=True)
    self.assertSameElements(
        filtered.chains, protein_chains + rna_chains + dna_chains
    )
    filtered = struc.filter_to_entity_type(ligand=True, dna=True)
    self.assertSameElements(filtered.chains, ligand_chains + dna_chains)

  def test_filter_out(self):
    struc = self.make_test_structure(
        chain_id=np.array(
            ['A', 'A', 'A', 'A', 'B', 'B', 'B', 'C', 'C'], dtype=object
        ),
        chain_type=np.array(
            ['polypeptide(L)'] * 7 + ['other'] * 2, dtype=object
        ),
        res_id=np.array([1, 1, 2, 2, 1, 2, 3, 1, 2], dtype=np.int32),
        atom_name=np.array(
            ['C', 'N', 'C', 'N', 'C', 'N', 'C', 'C', 'C'], dtype=object
        ),
    )
    with self.subTest('single_chain_id'):
      filtered = struc.filter_out(chain_id='A')
      self.assertSequenceEqual(filtered.chains, ('B', 'C'))
      self.assertEqual(filtered.num_atoms, 5)
    with self.subTest('multiple_chain_ids'):
      filtered = struc.filter_out(chain_id=('A', 'C'))
      self.assertSequenceEqual(filtered.chains, ('B',))
      self.assertEqual(filtered.num_atoms, 3)
    with self.subTest('chain_id_and_type_nonempty_selection'):
      filtered = struc.filter_out(chain_id='A', chain_type='polypeptide(L)')
      self.assertSequenceEqual(filtered.chains, ('B', 'C'))
      self.assertEqual(filtered.num_atoms, 5)
    with self.subTest('chain_id_and_type_empty_selection'):
      filtered = struc.filter_out(chain_id='A', chain_type='other')
      self.assertSequenceEqual(filtered.chains, ('A', 'B', 'C'))
      self.assertEqual(filtered.num_atoms, 9)
    with self.subTest('res_id'):
      filtered = struc.filter_out(res_id=1)
      self.assertSequenceEqual(filtered.chains, ('A', 'B', 'C'))
      self.assertEqual(filtered.num_atoms, 5)
      self.assertTrue(np.all(filtered.res_id != 1))
    with self.subTest('chain_id_and_res_id'):
      filtered = struc.filter_out(chain_id='A', res_id=1)
      self.assertSequenceEqual(filtered.chains, ('A', 'B', 'C'))
      np.testing.assert_array_equal(filtered.res_id, [2, 2, 1, 2, 3, 1, 2])
    with self.subTest('multiple_chain_ids_and_res_ids'):
      filtered = struc.filter_out(chain_id=('A', 'B'), res_id=(1, 3))
      self.assertSequenceEqual(filtered.chains, ('A', 'B', 'C'))
      np.testing.assert_array_equal(filtered.res_id, [2, 2, 2, 1, 2])
    with self.subTest('atom_name'):
      filtered = struc.filter_out(atom_name='N')
      self.assertSequenceEqual(filtered.chains, ('A', 'B', 'C'))
      self.assertEqual(filtered.num_atoms, 6)
      self.assertTrue(np.all(filtered.atom_name != 'N'))
    with self.subTest('chain_id_and_atom_name'):
      filtered = struc.filter_out(chain_id='A', atom_name='C')
      self.assertSequenceEqual(filtered.chains, ('A', 'B', 'C'))
      self.assertEqual(filtered.num_atoms, 7)
      np.testing.assert_array_equal(
          filtered.atom_name, ['N', 'N', 'C', 'N', 'C', 'C', 'C']
      )
    with self.subTest('chain_id_and_res_id_and_atom_name'):
      filtered = struc.filter_out(chain_id='A', res_id=1, atom_name='C')
      self.assertSequenceEqual(filtered.chains, ('A', 'B', 'C'))
      self.assertEqual(filtered.num_atoms, 8)
      np.testing.assert_array_equal(filtered.res_id, [1, 2, 2, 1, 2, 3, 1, 2])
      np.testing.assert_array_equal(
          filtered.atom_name, ['N', 'C', 'N', 'C', 'N', 'C', 'C', 'C']
      )
    with self.subTest('mask'):
      mask = np.array([1, 0, 0, 1, 1, 0, 1, 0, 1], dtype=bool)
      filtered = struc.filter_out(mask)
      self.assertSequenceEqual(filtered.chains, ('A', 'B', 'C'))
      self.assertEqual(filtered.num_atoms, 4)
      np.testing.assert_array_equal(filtered.res_id, [1, 2, 2, 1])
      np.testing.assert_array_equal(filtered.atom_name, ['N', 'C', 'N', 'C'])
    with self.subTest('mask_and_chain_id'):
      mask = np.array([1, 0, 0, 1, 1, 0, 1, 0, 1], dtype=bool)
      filtered = struc.filter_out(mask, chain_id='B')
      self.assertSequenceEqual(filtered.chains, ('A', 'B', 'C'))
      self.assertEqual(filtered.num_atoms, 7)
      np.testing.assert_array_equal(filtered.res_id, [1, 1, 2, 2, 2, 1, 2])
      np.testing.assert_array_equal(
          filtered.atom_name, ['C', 'N', 'C', 'N', 'N', 'C', 'C']
      )

  def test_chains(self):
    struc = structure.from_atom_arrays(
        chain_id=np.array(['A', 'A', 'B', 'D'], dtype=object),
        res_id=np.array([1, 2, 1, 1], dtype=np.int32),
        res_name=np.array(['ALA', 'ARG', 'CYS', 'PRO'], dtype=object),
        all_residues={
            'A': [('ALA', 1), ('ARG', 2), ('TYR', 3)],  # (TYR, 3) is missing.
            'B': [('CYS', 1)],
            'D': [('PRO', 1)],
        },
    )
    self.assertSequenceEqual(struc.chains, ('A', 'B', 'D'))

  def test_boundaries(self):
    mmcif_str = test_utils.read_test_cif_file('2pv7')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)

    self.assertLen(struc.chain_boundaries, struc.num_chains)
    self.assertLen(
        struc.res_boundaries, struc.num_residues(count_unresolved=False)
    )

    with self.subTest('chain_boundaries'):
      np.testing.assert_array_equal(
          struc.chain_boundaries, (0, 2205, 4427, 4440, 4484, 4497)
      )
    # Test the start and the end for residue boundaries.
    with self.subTest('res_boundaries'):
      np.testing.assert_array_equal(
          struc.res_boundaries[:5], (0, 4, 15, 24, 31)
      )
      np.testing.assert_array_equal(
          struc.res_boundaries[-5:], (4425, 4427, 4440, 4484, 4497)
      )
      self.assertContainsSubset(
          set(struc.chain_boundaries), set(struc.res_boundaries)
      )

  @parameterized.parameters([1, 2, 10])
  def test_transform_states_to_chains(self, num_states):
    multi_state_struc = self.make_test_structure(
        name='test',
        all_residues={'A': [('ALA', 1), ('MSE', 2)], 'B': [('VWW', 1)]},
        chain_id=np.array(['A', 'A', 'B'], dtype=object),
        chain_type=np.array(['polypeptide(L)'] * 3, dtype=object),
        res_id=np.array([1, 2, 1], dtype=np.int32),
        res_name=np.array(['ALA', 'MSE', 'VWW'], dtype=object),
        atom_name=np.array(['C', 'C', 'C'], dtype=object),
        atom_element=np.array(['C', 'C', 'C'], dtype=object),
        atom_x=np.zeros(shape=(num_states, 3), dtype=np.float32),
        atom_y=np.zeros(shape=(num_states, 3), dtype=np.float32),
        atom_z=np.zeros(shape=(num_states, 3), dtype=np.float32),
        atom_b_factor=np.zeros(shape=(num_states, 3), dtype=np.float32),
        atom_occupancy=np.zeros(shape=(num_states, 3), dtype=np.float32),
    )
    multi_chain_struc = multi_state_struc.transform_states_to_chains()
    self.assertEqual(
        multi_chain_struc.num_chains, num_states * multi_state_struc.num_chains
    )
    self.assertEqual(multi_chain_struc.num_models, 1)
    all_chain_ids = multi_chain_struc.chains
    with self.subTest('compare_each_slice'):
      for i in range(num_states):
        start_chain_idx = i * multi_state_struc.num_chains
        end_chain_idx = start_chain_idx + multi_state_struc.num_chains
        chain_id = all_chain_ids[start_chain_idx:end_chain_idx]
        sliced_struc = multi_state_struc.slice_leading_dims[i]
        filtered_struc = multi_chain_struc.filter(chain_id=chain_id)
        renaming = {
            c: new_c
            for c, new_c in zip(chain_id, multi_state_struc.chains, strict=True)
        }
        renamed_struc = filtered_struc.rename_chain_ids(renaming)
        self.assertStructuresEqual(
            renamed_struc,
            sliced_struc,
            ignore_fields=['name', 'author_naming_scheme'],
        )
    with self.subTest('to_mmcif'):
      _ = multi_chain_struc.to_mmcif()

  @parameterized.named_parameters(
      ('default_chain_ids', None, ('A', 'B', 'C', 'D', 'E', 'F')),
      (
          'with_chain_ids',
          ('Z', 'B', 'X', 'A', 'D', 'F'),
          ('Z', 'B', 'X', 'A', 'D', 'F'),
      ),
  )
  def test_merge_chains(self, group_chain_ids, exp_chain_ids):
    mmcif_str = test_utils.read_test_cif_file('2ppb_trimmed')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)

    merged = struc.merge_chains(
        chain_groups=(
            ('A', 'C', 'D', 'F'),  # DNA.
            ('B', 'E'),  # RNA.
            ('G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'),  # Protein.
            ('P',),  # Single protein.
            ('Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'),  # Ligands.
            ('AA', 'BA'),  # More ligands.
        ),
        chain_group_ids=group_chain_ids,
    )
    with self.subTest('chains'):
      self.assertEqual(merged.num_chains, 6)
      self.assertEqual(merged.chains, exp_chain_ids)
      self.assertEqual(
          merged.chains_table.type.tolist(),
          [mmcif_names.DNA_CHAIN, mmcif_names.RNA_CHAIN]
          + [mmcif_names.PROTEIN_CHAIN] * 2
          + [mmcif_names.NON_POLYMER_CHAIN] * 2,
      )
    with self.subTest('residues'):
      self.assertEqual(
          merged.num_residues(count_unresolved=True),
          struc.num_residues(count_unresolved=True),
      )
      self.assertEqual(
          merged.num_residues(count_unresolved=False),
          struc.num_residues(count_unresolved=False),
      )
    with self.subTest('atoms'):
      self.assertEqual(merged.num_atoms, struc.num_atoms)
      self.assertEqual(
          struc.filter(chain_id=('B', 'E')).num_atoms,
          merged.filter(chain_id='B').num_atoms,
      )
    with self.subTest('residue_ids_unique_in_chains'):
      original_res_count = len(
          set(zip(struc.residues_table.chain_key, struc.residues_table.id))
      )
      merged_res_count = len(
          set(zip(merged.residues_table.chain_key, merged.residues_table.id))
      )
      self.assertEqual(original_res_count, merged_res_count)
      for chain_id in merged.chains:
        merged_chain_struc = merged.filter(chain_id=chain_id)
        exp_unique_res_ids = list(
            range(1, len(merged_chain_struc.residues_table) + 1)
        )
        self.assertEqual(
            merged_chain_struc.residues_table.id.tolist(),
            exp_unique_res_ids,
            msg=f'merged chain {chain_id}',
        )
        self.assertEqual(
            merged_chain_struc.residues_table.auth_seq_id.tolist(),
            [str(res_id) for res_id in exp_unique_res_ids],
            msg=f'merged chain {chain_id}',
        )

  def test_merge_chains_set_group_properties(self):
    mmcif_str = test_utils.read_test_cif_file('2ppb_trimmed')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)
    merged = struc.merge_chains(
        chain_groups=[struc.chains],
        chain_group_types=[mmcif_names.DNA_CHAIN],
        chain_group_entity_ids=['42'],
    )
    self.assertEqual(merged.chains, ('A',))
    self.assertEqual(merged.chains_table.type.tolist(), [mmcif_names.DNA_CHAIN])
    self.assertEqual(merged.chains_table.entity_id.tolist(), ['42'])
    self.assertEqual(
        merged.residues_table.id.tolist(),
        list(range(1, len(merged.residues_table) + 1)),
    )

  def test_merge_chains_no_op(self):
    mmcif_str = test_utils.read_test_cif_file('2ppb_trimmed')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)
    merged = struc.merge_chains(chain_groups=[[cid] for cid in struc.chains])
    self.assertStructuresEqual(merged, struc)

  def test_merge_chains_bad_ids(self):
    mmcif_str = test_utils.read_test_cif_file('2ppb_trimmed')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)
    with self.assertRaises(ValueError):
      struc.merge_chains(chain_groups=[('A', 'C', 'D', 'F')])

  def test_merge_chains_bad_groups(self):
    mmcif_str = test_utils.read_test_cif_file('2ppb_trimmed')
    struc = structure.from_mmcif(mmcif_str, include_bonds=False)
    with self.assertRaises(ValueError):
      struc.merge_chains(chain_groups=[struc.chains])

  @parameterized.parameters(
      ('2pv7', atom_types.ATOM37, atom_types.ATOM37_ORDER, 'A'),
      ('402d', atom_types.ATOM29, atom_types.ATOM29_ORDER, 'B'),
  )
  def test_to_res_arrays(self, pdb_id, atom_name, atom_order, chain_id):
    mmcif_str = test_utils.read_test_cif_file(pdb_id)
    full_struc = structure.from_mmcif(mmcif_str, include_bonds=False)
    # Filter to chain A and only standard protein atom types.
    orig_chain_struc = full_struc.filter(chain_id=chain_id, atom_name=atom_name)
    res_arrays = orig_chain_struc.to_res_arrays(
        include_missing_residues=True, atom_order=atom_order
    )
    atom_positions = res_arrays.atom_positions
    atom_mask = res_arrays.atom_mask
    res_name, res_id = zip(*full_struc.all_residues[chain_id], strict=True)
    num_res = len(res_name)
    atom_name = np.tile(np.array(atom_name, dtype=object), (num_res, 1))
    struc = structure.from_res_arrays(
        name=pdb_id.upper(),
        atom_mask=atom_mask,
        atom_x=atom_positions[..., 0],
        atom_y=atom_positions[..., 1],
        atom_z=atom_positions[..., 2],
        atom_name=atom_name,
        res_name=np.array(res_name, dtype=object),
        res_id=np.array(res_id, dtype=np.int32),
        chain_id=np.full((num_res,), chain_id, dtype=object),
    )
    self.assertStructuresEqual(
        struc,
        orig_chain_struc,
        ignore_fields=[
            'resolution',
            'release_date',
            'structure_method',
            'author_naming_scheme',
            'all_residues',
            'bioassembly_data',
            'chemical_components_data',
            'bonds',
            'atom_element',
            'atom_b_factor',
            'atom_occupancy',
            'chain_type',
        ],
    )
    self.assertEqual(atom_positions.shape, (num_res, len(atom_order), 3))
    self.assertEqual(atom_mask.shape, (num_res, len(atom_order)))

  @parameterized.parameters(True, False)
  def test_to_res_arrays_multi_chain(self, include_missing_residues):
    struc = structure.from_mmcif(
        test_utils.read_test_cif_file('2pv7'), include_bonds=False
    )
    self.assertGreater(struc.num_chains, 1)
    expected = [
        chain.to_res_arrays(include_missing_residues=include_missing_residues)
        for chain in struc.split_by_chain()
    ]
    actual = struc.to_res_arrays(
        include_missing_residues=include_missing_residues
    )
    np.testing.assert_array_equal(
        actual.atom_positions,
        np.concatenate([c.atom_positions for c in expected], axis=0),
        strict=True,
    )
    np.testing.assert_array_equal(
        actual.atom_mask,
        np.concatenate([c.atom_mask for c in expected], axis=0),
        strict=True,
    )
    np.testing.assert_array_equal(
        actual.atom_b_factor,
        np.concatenate([c.atom_b_factor for c in expected], axis=0),
        strict=True,
    )
    np.testing.assert_array_equal(
        actual.atom_occupancy,
        np.concatenate([c.atom_occupancy for c in expected], axis=0),
        strict=True,
    )

  @parameterized.parameters(True, False)
  def test_to_res_arrays_unresolved_leading_chain(
      self, include_missing_residues
  ):
    struc = structure.from_mmcif(
        test_utils.read_test_cif_file('2pv7'), include_bonds=False
    )
    zero_first = struc.filter(
        struc.chain_id != struc.chains[0],
        cascade_delete=structure.CascadeDelete.NONE,
    )
    actual = zero_first.to_res_arrays(
        include_missing_residues=include_missing_residues
    )
    remaining = struc.filter(struc.chain_id != struc.chains[0]).to_res_arrays(
        include_missing_residues=include_missing_residues
    )
    offset = (
        len(struc.all_residues[struc.chains[0]])
        if include_missing_residues
        else 0
    )
    np.testing.assert_array_equal(
        actual.atom_mask[:offset],
        np.zeros((offset, atom_types.ATOM37_NUM), dtype=np.float32),
        strict=True,
    )
    np.testing.assert_array_equal(
        actual.atom_positions[offset:],
        remaining.atom_positions,
        strict=True,
    )
    np.testing.assert_array_equal(
        actual.atom_mask[offset:],
        remaining.atom_mask,
        strict=True,
    )
    np.testing.assert_array_equal(
        actual.atom_b_factor[offset:],
        remaining.atom_b_factor,
        strict=True,
    )
    np.testing.assert_array_equal(
        actual.atom_occupancy[offset:],
        remaining.atom_occupancy,
        strict=True,
    )

  @parameterized.parameters(True, False)
  def test_to_res_arrays_empty_structure(self, include_missing_residues):
    res_arrays = structure.make_empty_structure().to_res_arrays(
        include_missing_residues=include_missing_residues
    )
    self.assertEqual(
        res_arrays.atom_positions.shape, (0, atom_types.ATOM37_NUM, 3)
    )
    self.assertEqual(res_arrays.atom_mask.shape, (0, atom_types.ATOM37_NUM))
    self.assertEqual(res_arrays.atom_b_factor.shape, (0, atom_types.ATOM37_NUM))
    self.assertEqual(
        res_arrays.atom_occupancy.shape, (0, atom_types.ATOM37_NUM)
    )

  @parameterized.parameters(True, False)
  def test_to_res_arrays_sparse_res_keys(self, include_missing_residues):
    struc = structure.from_mmcif(
        test_utils.read_test_cif_file('2pv7'),
        model_id=structure.ModelID.FIRST,
        include_bonds=False,
    )
    sparse_keys = struc.residues_table.key * 100_000 + 10_000_000
    perm = np.arange(struc.residues_table.size)[::-1]
    sparse = struc.copy_and_update(
        residues=struc.residues_table.apply_array(perm).copy_and_update(
            key=sparse_keys[perm]
        ),
        atoms=struc.atoms_table.copy_and_update(
            res_key=sparse_keys[struc.atoms_table.res_key]
        ),
        skip_validation=True,
    )
    expected = struc.to_res_arrays(
        include_missing_residues=include_missing_residues
    )
    actual = sparse.to_res_arrays(
        include_missing_residues=include_missing_residues
    )
    expected_perm = np.arange(expected.atom_positions.shape[0])[::-1]
    np.testing.assert_array_equal(
        actual.atom_positions,
        expected.atom_positions[expected_perm],
        strict=True,
    )
    np.testing.assert_array_equal(
        actual.atom_mask, expected.atom_mask[expected_perm], strict=True
    )
    np.testing.assert_array_equal(
        actual.atom_b_factor, expected.atom_b_factor[expected_perm], strict=True
    )
    np.testing.assert_array_equal(
        actual.atom_occupancy,
        expected.atom_occupancy[expected_perm],
        strict=True,
    )
    self.assertNotIn('index_by_key', sparse.residues_table.__dict__)
    self.assertNotIn('index_by_key', sparse.present_residues.__dict__)

  def test_to_res_arrays_single_model_leading_dim_raises(self):
    # 2pv7 is a single-model mmCIF, so ModelID.ALL adds a size-1 leading dim.
    batched = structure.from_mmcif(
        test_utils.read_test_cif_file('2pv7'),
        model_id=structure.ModelID.ALL,
        include_bonds=False,
    )
    self.assertEqual(batched.atoms_table.ndim, 2)
    self.assertEqual(batched.num_models, 1)
    with self.assertRaisesRegex(ValueError, 'Only 1D atom arrays allowed'):
      batched.to_res_arrays(include_missing_residues=False)

  def test_to_res_arrays_multi_model_raises(self):
    multi_model = structure.from_mmcif(
        test_utils.read_test_cif_file('1wib'),
        model_id=structure.ModelID.ALL,
        include_bonds=False,
    )
    with self.assertRaisesRegex(
        ValueError, 'Only single-model Structure allowed'
    ):
      multi_model.to_res_arrays(include_missing_residues=False)

  @parameterized.parameters(
      ('2pv7', structure.ModelID.FIRST, True, 600, 40, ()),
      ('2pv7', structure.ModelID.FIRST, False, 560, 0, ()),
      ('1wib', structure.ModelID.ALL, True, 92, 0, (20,)),
      ('1wib', structure.ModelID.ALL, False, 92, 0, (20,)),
  )
  def test_to_res_atom_lists(
      self,
      pdb_id,
      model_id,
      include_missing_residues,
      exp_len,
      exp_num_missing,
      exp_shape,
  ):
    mmcif_str = test_utils.read_test_cif_file(pdb_id)
    struc = structure.from_mmcif(
        mmcif_str, model_id=model_id, include_bonds=False
    )
    res_atom_lists = struc.to_res_atom_lists(
        include_missing_residues=include_missing_residues
    )
    self.assertLen(res_atom_lists, exp_len)
    self.assertLen([r for r in res_atom_lists if not r], exp_num_missing)
    self.assertEqual(sum([len(r) for r in res_atom_lists]), struc.num_atoms)

    for residue_atoms in res_atom_lists:
      if residue_atoms:
        self.assertEqual(
            residue_atoms[0]['atom_x'].shape, exp_shape, msg=residue_atoms[0]
        )

  @parameterized.product(
      pdb_id_model_id=(
          ('2pv7', structure.ModelID.FIRST),
          ('2ppb_trimmed', structure.ModelID.FIRST),
          ('1wib', structure.ModelID.ALL),
          ('3b9f', structure.ModelID.FIRST),
      ),
      use_auth_naming_scheme=(True, False),
  )
  def test_to_mmcif(self, pdb_id_model_id, use_auth_naming_scheme):
    pdb_id, model_id = pdb_id_model_id

    original_mmcif_str = test_utils.read_test_cif_file(pdb_id)
    original_struc = structure.from_mmcif(
        original_mmcif_str, model_id=model_id, include_bonds=False
    )
    if not use_auth_naming_scheme:
      original_struc = original_struc.reset_author_naming_scheme()

    generated_mmcif_str = original_struc.to_mmcif()
    struc = structure.from_mmcif(generated_mmcif_str, model_id=model_id)

    if not use_auth_naming_scheme:
      struc = struc.reset_author_naming_scheme()
    self.assertStructuresEqual(
        struc,
        original_struc,
        allow_chem_comp_data_extension=True,
    )

  def test_to_mmcif_water_in_nonpoly(self):
    struc = structure.from_mmcif(
        test_utils.read_test_cif_file('3b9f'),
        include_water=True,
        include_bonds=False,
    )
    cif = struc.to_mmcif_dict()

    self.assertIn('HOH', cif['_pdbx_nonpoly_scheme.mon_id'])
    self.assertNotIn('HOH', cif['_pdbx_branch_scheme.mon_id'])

  def test_atom_record_name(self):
    protein = residue_names.PROTEIN_TYPES_WITH_UNKNOWN + (residue_names.MSE,)
    rna = residue_names.RNA_TYPES + (residue_names.UNK_RNA,)
    dna = residue_names.DNA_TYPES + (residue_names.UNK_DNA,)
    ligand = ('NAD', 'ATP')
    water = residue_names.WATER_TYPES
    res_name = protein + rna + dna + ligand + water
    chain_id = (
        ['A'] * len(protein)
        + ['B'] * len(rna)
        + ['C'] * len(dna)
        + ['D'] * len(ligand)
        + ['E'] * len(water)
    )
    chain_type = (
        [mmcif_names.PROTEIN_CHAIN] * len(protein)
        + [mmcif_names.RNA_CHAIN] * len(rna)
        + [mmcif_names.DNA_CHAIN] * len(dna)
        + [mmcif_names.NON_POLYMER_CHAIN] * len(ligand)
        + [mmcif_names.WATER] * len(water)
    )

    num_res = len(res_name)
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(chain_id, dtype=object),
        chain_type=np.array(chain_type, dtype=object),
        res_id=np.arange(num_res, dtype=np.int32),
        res_name=np.array(res_name, dtype=object),
        atom_name=np.array(['C'] * num_res, dtype=object),
        atom_element=np.array(['C'] * num_res, dtype=object),
        atom_x=np.array([0.0] * num_res, dtype=np.float32),
        atom_y=np.array([0.0] * num_res, dtype=np.float32),
        atom_z=np.array([0.0] * num_res, dtype=np.float32),
        atom_b_factor=np.array([0.0] * num_res, dtype=np.float32),
    )
    exp_record_name = (
        (['ATOM'] * (len(protein) - 2) + ['HETATM'] * 2)
        + (['ATOM'] * (len(rna) - 1) + ['HETATM'] * 1)
        + (['ATOM'] * (len(dna) - 1) + ['HETATM'] * 1)
        + (['HETATM'] * len(ligand))
        + (['HETATM'] * len(water))
    )

    cif = struc.to_mmcif_dict()
    self.assertEqual(cif['_atom_site.group_PDB'], exp_record_name)

  @parameterized.parameters((1, ['1.0']), (2, ['1.03']), (4, ['1.0258']))
  def test_to_mmcif_coords_decimal_places(
      self, coords_decimal_places, exp_atom_y
  ):
    struc = self.make_test_structure(
        name='test_structure',
        chain_id=np.array(['A'], dtype=object),
        chain_type=np.array(['polypeptide(L)'], dtype=object),
        res_id=np.array([1], dtype=np.int32),
        res_name=np.array(['ALA'], dtype=object),
        atom_name=np.array(['C'], dtype=object),
        atom_element=np.array(['C'], dtype=object),
        atom_x=np.array([1.0], dtype=np.float32),
        atom_y=np.array([1.02583], dtype=np.float32),
        atom_z=np.array([1.0], dtype=np.float32),
        atom_b_factor=np.array([1.0], dtype=np.float32),
    )
    cif = struc.to_mmcif_dict(coords_decimal_places=coords_decimal_places)
    exp_atom_x_z = ['1.' + '0' * coords_decimal_places]
    self.assertEqual(cif['_atom_site.Cartn_x'], exp_atom_x_z)
    self.assertEqual(cif['_atom_site.Cartn_z'], exp_atom_x_z)
    self.assertEqual(cif['_atom_site.Cartn_y'], exp_atom_y)

  def test_to_mmcif_with_flat_b_factor(self):
    mmcif_str = test_utils.read_test_cif_file('1wib')
    struc = structure.from_mmcif(
        mmcif_str, model_id=structure.ModelID.ALL, include_bonds=False
    )
    struc = struc.copy_and_update_atoms(atom_b_factor=struc.atom_b_factor[0])
    new_mmcif = struc.to_mmcif()
    new_struc = structure.from_mmcif(
        new_mmcif, model_id=structure.ModelID.ALL, include_bonds=False
    )
    new_struc = new_struc.copy_and_update_atoms(
        atom_b_factor=new_struc.atom_b_factor[0]
    )
    self.assertStructuresEqual(new_struc, struc)

  def test_to_mmcif_no_all_residues(self):
    original_mmcif_str = test_utils.read_test_cif_file('2pv7')
    original_struc = structure.from_mmcif(
        original_mmcif_str, include_bonds=False
    )
    original_struc_no_all_residues = original_struc.copy_and_update(
        residues=original_struc.present_residues
    )
    generated_mmcif_str = original_struc_no_all_residues.to_mmcif()
    struc = structure.from_mmcif(generated_mmcif_str, include_bonds=False)

    # Residue 1 in chain A is unresolved so all_residues should start from
    # 13 instead when all_residues has been removed.
    self.assertTupleEqual(struc.all_residues['A'][0], ('GLY', 13))
    self.assertTupleEqual(original_struc.all_residues['A'][0], ('GLY', 1))

    self.assertTupleEqual(struc.all_residues['A'][-1], ('ALA', 292))
    self.assertTupleEqual(original_struc.all_residues['A'][-1], ('GLY', 298))

  def test_to_mmcif_missing_pdbx_poly_seq_scheme(self):
    original_mmcif_str = test_utils.read_test_cif_file(
        '2pv7_missing_pdbx_poly_seq_scheme'
    )
    original_struc = structure.from_mmcif(
        original_mmcif_str, include_bonds=False
    )
    generated_mmcif_str = original_struc.to_mmcif()

    struc = structure.from_mmcif(generated_mmcif_str, include_bonds=False)
    with self.subTest('round_trips_to_mmcif'):
      self.assertStructuresEqual(struc, original_struc)

    # The first residue in chain A is unresolved. Since this mmCIF is missing
    # a _pdbx_poly_seq_scheme table we should still capture this residue in
    # all_residues via _entity_poly_seq but we will lack auth_seq_ids and
    # insertion codes.
    with self.subTest('present_in_all_residues'):
      self.assertTupleEqual(original_struc.all_residues['A'][0], ('GLY', 1))
      self.assertTupleEqual(struc.all_residues['A'][0], ('GLY', 1))
    with self.subTest('missing_from_author_naming_scheme'):
      self.assertNotIn(1, original_struc.author_naming_scheme.auth_seq_id['A'])
      self.assertNotIn(1, struc.author_naming_scheme.auth_seq_id['A'])
      self.assertNotIn(
          1, original_struc.author_naming_scheme.insertion_code['A']
      )
      self.assertNotIn(1, struc.author_naming_scheme.insertion_code['A'])

  def test_multichain_residue_index(self):
    struc = self.make_test_structure(
        chain_id=np.array(
            ['B', 'B', 'B', 'A', 'A', 'G', 'E', 'E'], dtype=object
        ),
        res_id=np.array([0, 1, 1, 0, 1, 5, 2, 3], dtype=np.int32),
        atom_name=np.array(
            ['C', 'C', 'CA', 'C', 'C', 'C', 'C', 'C'], dtype=object
        ),
    )
    multichain_residue_index = structure.multichain_residue_index(struc)
    np.testing.assert_array_equal(
        multichain_residue_index, [0, 1, 1, 10000, 10001, 20005, 30002, 30003]
    )
    self.assertTrue(
        np.all(multichain_residue_index[1:] >= multichain_residue_index[:-1]),
        msg=multichain_residue_index,
    )

  def test_stack_empty(self):
    with self.assertRaises(ValueError):
      _ = structure.stack([])

  def test_stack_incompatible_strucs(self):
    struc_0 = self.make_test_structure(
        res_id=np.arange(4, dtype=np.int32),
        atom_x=np.zeros([4], dtype=np.float32),
        atom_y=np.ones([4], dtype=np.float32),
        atom_z=np.ones([4], dtype=np.float32),
        atom_name=np.array(['CA'] * 4, dtype=object),
    )
    struc_1 = self.make_test_structure(
        res_id=np.arange(5, dtype=np.int32),
        atom_x=np.ones([5], dtype=np.float32),
        atom_y=np.ones([5], dtype=np.float32),
        atom_z=np.ones([5], dtype=np.float32),
        atom_name=np.array(['CA'] * 5, dtype=object),
    )
    with self.subTest('atom_shape_mismatch'):
      with self.assertRaises(ValueError):
        _ = structure.stack([struc_0, struc_1])
    struc_0_atoms_renamed = struc_0.copy_and_update_atoms(
        atom_name=np.array(['N'] * 4, dtype=object)
    )
    with self.subTest('atom_name_mismatch'):
      with self.assertRaises(ValueError):
        _ = structure.stack([struc_0, struc_0_atoms_renamed])
    struc_0_with_bfac = struc_0.copy_and_update_atoms(
        atom_b_factor=np.ones([4]), atom_occupancy=np.ones([4])
    )
    stacked = structure.stack([struc_0_with_bfac, struc_0_with_bfac])
    stacked_shared_bfac = stacked.copy_and_update_atoms(
        atom_b_factor=np.ones([4]), atom_occupancy=np.ones([4])
    )
    with self.subTest('b_fac_shape_mismatch'):
      with self.assertRaises(ValueError):
        _ = structure.stack([stacked, stacked_shared_bfac])

  def test_stack(self):
    coords = np.zeros((5, 4, 3), dtype=np.float32)
    strucs = []
    for i in range(6):
      this_coords = coords + i
      struc = self.make_test_structure(
          res_id=np.arange(4, dtype=np.int32),
          atom_name=np.array(['CA'] * 4, dtype=object),
          atom_x=this_coords[..., 0],
          atom_y=this_coords[..., 1],
          atom_z=this_coords[..., 2],
          atom_b_factor=this_coords[..., 0],
          atom_occupancy=np.ones_like(this_coords[..., 0], dtype=np.float32),
      )
      strucs.append(struc)
    stacked = structure.stack(strucs, axis=1)
    self.assertTupleEqual(stacked.atom_x.shape, (5, 6, 4))
    for i in range(6):
      struc = stacked.slice_leading_dims[:, i]
      expected = np.zeros((5, 4)) + i
      for field in ['atom_x', 'atom_y', 'atom_z', 'atom_b_factor']:
        np.testing.assert_array_equal(
            struc[field], expected, err_msg=f'{i} {field}'
        )
      np.testing.assert_array_equal(
          struc.atom_name, np.array(['CA'] * 4), err_msg=f'{i} atom_name'
      )

  def test_stack_unstack_roundtrip(self):
    mmcif_str = test_utils.read_test_cif_file('1wib')
    struc = structure.from_mmcif(
        mmcif_str, model_id=structure.ModelID.ALL, include_bonds=False
    )
    strucs = [struc.slice_leading_dims[1], struc.slice_leading_dims[5]]
    stacked = structure.stack(strucs)
    unstacked = stacked.unstack()
    self.assertLen(unstacked, len(strucs))
    for orig, new in zip(strucs, unstacked, strict=True):
      self.assertStructuresEqual(orig, new)

  def test_unstack_stack_roundtrip(self):
    mmcif_str = test_utils.read_test_cif_file('1wib')
    struc = structure.from_mmcif(
        mmcif_str, model_id=structure.ModelID.ALL, include_bonds=False
    )
    unstacked = struc.unstack()
    stacked = structure.stack(unstacked)
    self.assertStructuresEqual(struc, stacked)

  def test_assign_unique_chain_ids(self):
    chain_ids = [
        ['C', 'D'],
        ['A', 'B'],
        ['D', 'A', 'E'],
        ['F', 'C', 'A'],
    ]
    chain_ids = [np.array(chain_id) for chain_id in chain_ids]
    strucs = [
        self.make_test_structure(
            chain_id=np.array(chain_id, dtype=object),
            res_id=np.ones(chain_id.shape, dtype=np.int32),
        )
        for chain_id in chain_ids
    ]
    renamed_strucs = structure_lib._assign_unique_chain_ids(strucs)
    self.assertLen(renamed_strucs, len(strucs))
    seen_chains = set()
    for struc, renamed_struc in zip(strucs, renamed_strucs, strict=True):
      self.assertLen(struc.chain_id, len(renamed_struc.chain_id))
      renamed_chains = renamed_struc.chains
      self.assertEmpty(seen_chains.intersection(renamed_chains))
      seen_chains.update(renamed_chains)

  def test_concat_raises(self):
    chain_ids = [['A', 'B'], ['C'], ['D', 'B']]
    strucs = [
        self.make_test_structure(
            chain_id=np.array(chain_id, dtype=object),
            res_name=np.array(['ALA'] * len(chain_id), dtype=object),
            res_id=np.arange(len(chain_id), dtype=np.int32),
        )
        for chain_id in chain_ids
    ]
    with self.subTest('empty'):
      with self.assertRaises(ValueError):
        _ = structure.concat([])
    with self.subTest('duplicate_chain_ids'):
      with self.assertRaises(ValueError):
        _ = structure.concat(strucs, assign_unique_chain_ids=False)

  def test_concat_fixes_entity_id(self):
    strucs = [
        self.make_test_structure(
            chain_id=np.array(['A', 'A'], dtype=object),
            res_name=np.array(['ALA', 'ALA'], dtype=object),
            res_id=np.arange(2, dtype=np.int32),
        ),
        self.make_test_structure(
            chain_id=np.array(['B', 'B'], dtype=object),
            res_name=np.array(['PRO', 'PRO'], dtype=object),
            res_id=np.arange(2, dtype=np.int32),
        ),
    ]
    concatted = structure.concat(
        strucs, assign_unique_chain_ids=False, assign_unique_entity_ids=True
    )
    self.assertCountEqual(concatted.chains_table.entity_id, ['1', '2'])

  @parameterized.parameters(
      dict(
          chain_ids=[['A']],
          names=['a'],
          assign_unique_chain_ids=True,
          expected_chains=['A'],
          expected_name='a',
      ),
      dict(
          chain_ids=[['B'], ['A']],
          names=['foo', 'bar'],
          assign_unique_chain_ids=False,
          expected_chains=['B', 'A'],
          expected_name='foo_bar',
      ),
      dict(
          chain_ids=[['A', 'C'], ['B'], ['D', 'E']],
          names=['all', 'bar', None],
          assign_unique_chain_ids=False,
          expected_chains=['A', 'C', 'B', 'D', 'E'],
          expected_name='all_bar_None',
      ),
      dict(
          chain_ids=[['Z'], ['C'], ['B', 'C']],
          names=[None, 'a', None],
          assign_unique_chain_ids=True,
          expected_chains=['A', 'B', 'C', 'D'],
          expected_name='None_a_None',
      ),
  )
  def test_concat(
      self,
      chain_ids,
      names,
      assign_unique_chain_ids,
      expected_chains,
      expected_name,
  ):
    assert len(chain_ids) == len(names)
    strucs = []
    for chain_id, name in zip(chain_ids, names, strict=True):
      num_chains = len(chain_id)
      struc = self.make_test_structure(
          name=name,
          chain_id=np.repeat(chain_id, 3).astype(object),
          res_name=np.array(['ALA', 'PRO', 'HIS'] * num_chains, dtype=object),
          res_id=np.array([2, 5, 6] * num_chains, dtype=np.int32),
          atom_key=np.arange(num_chains * 3, dtype=np.int64),
          atom_x=np.zeros((4, 2, 3 * num_chains), dtype=np.float32),
          atom_name=np.array(['CA', 'CA', 'CA'] * num_chains, dtype=object),
          all_residues={
              c: [
                  ('GLY', 1),
                  ('ALA', 2),
                  ('GLY', 3),
                  ('GLY', 4),
                  ('PRO', 5),
                  ('HIS', 6),
              ]
              for c in chain_id
          },
          bonds=structure.Bonds(
              key=np.array([1]),
              type=np.array(['type'], dtype=object),
              role=np.array(['role'], dtype=object),
              from_atom_key=np.array([1]),
              dest_atom_key=np.array([0]),
          ),
      )
      strucs.append(struc)
    concatted = structure.concat(
        strucs, assign_unique_chain_ids=assign_unique_chain_ids
    )
    with self.subTest('name'):
      self.assertEqual(expected_name, concatted.name)
    with self.subTest('chains'):
      self.assertSequenceEqual(expected_chains, concatted.chains)
      self.assertSameElements(expected_chains, concatted.all_residues.keys())
    with self.subTest('filter_to_each_chain_and_check_equality'):
      if assign_unique_chain_ids:
        strucs = structure_lib._assign_unique_chain_ids(strucs)
      for struc in strucs:
        for c in struc.chains:
          orig_chain_struc = struc.filter(chain_id=c)
          concat_chain_struc = concatted.filter(chain_id=c)
          self.assertStructuresEqual(
              orig_chain_struc,
              concat_chain_struc,
              ignore_fields=[
                  'chain_id',
                  'name',
                  'release_date',
                  'resolution',
                  'structure_method',
                  'author_naming_scheme',
                  'all_residues',
                  'bioassembly_data',
              ],
          )

  def test_concat_mixed_all_res_present(self):
    strucs = [
        self.make_test_structure(
            chain_id=np.array(['A', 'A'], dtype=object),
            res_name=np.array(['ALA', 'GLY'], dtype=object),
            res_id=np.array([1, 2], dtype=np.int32),
            all_residues={'A': [('ALA', 1), ('GLY', 2), ('MSE', 3)]},
        ),
        self.make_test_structure(
            chain_id=np.array(['A', 'A'], dtype=object),
            res_name=np.array(['ALA', 'GLY'], dtype=object),
            res_id=np.array([1, 2], dtype=np.int32),
        ),
    ]
    concatted = structure.concat(strucs)
    self.assertDictEqual(
        concatted.all_residues,
        {
            'A': [('ALA', 1), ('GLY', 2), ('MSE', 3)],
            'B': [('ALA', 1), ('GLY', 2)],
        },
    )

  def test_concat_mixed_bonds_present(self):
    strucs = [
        self.make_test_structure(
            chain_id=np.array(['A', 'A', 'A'], dtype=object),
            res_id=np.array([1, 2, 3], dtype=np.int32),
            res_name=np.array(['ALA', 'GLY', 'MSE'], dtype=object),
            atom_key=np.array([0, 1, 2], dtype=np.int64),
            bonds=structure.Bonds(
                key=np.array([1]),
                type=np.array(['type'], dtype=object),
                role=np.array(['role'], dtype=object),
                from_atom_key=np.array([1]),
                dest_atom_key=np.array([0]),
            ),
        ),
        self.make_test_structure(
            chain_id=np.array(['A', 'A', 'A'], dtype=object),
            res_id=np.array([1, 2, 3], dtype=np.int32),
            res_name=np.array(['ALA', 'GLY', 'MSE'], dtype=object),
            atom_key=np.array([0, 1, 2]),
        ),
        self.make_test_structure(
            chain_id=np.array(['A', 'A', 'A'], dtype=object),
            res_id=np.array([1, 2, 3], dtype=np.int32),
            res_name=np.array(['ALA', 'GLY', 'MSE'], dtype=object),
        ),
    ]
    concatted = structure.concat(strucs)
    self.assertEqual(concatted.bonds.size, 1)
    self.assertLen(concatted.atom_key, 9)

  @parameterized.named_parameters(
      ('default_name', None, '2PV7_1WIB'),
      ('give_name', 'cats', 'cats'),
  )
  def test_concat_real_data(self, name, expected_name):
    strucs = [
        structure.from_mmcif(
            test_utils.read_test_cif_file('2pv7'), include_bonds=False
        ),
        structure.from_mmcif(
            test_utils.read_test_cif_file('1wib'), include_bonds=False
        ),
    ]
    concatted = structure.concat(strucs, name=name)
    renamed = structure_lib._assign_unique_chain_ids(strucs)

    with self.subTest('filter_to_each_chain_and_check_equality'):
      for struc in renamed:
        for c in struc.chains:
          orig_chain_struc = struc.filter(chain_id=c)
          concat_chain_struc = concatted.filter(chain_id=c)
          self.assertStructuresEqual(
              orig_chain_struc,
              concat_chain_struc,
              ignore_fields={
                  'name',
                  'release_date',
                  'resolution',
                  'structure_method',
                  'author_naming_scheme',
                  'all_residues',
                  'bioassembly_data',
                  'chemical_components_data',
                  'bonds',
              },
          )

    with self.subTest('global_fields'):
      self.assertIsNone(concatted.release_date)
      self.assertIsNone(concatted.resolution)
      self.assertIsNone(concatted.structure_method)
      self.assertIsNone(concatted.bioassembly_data)
      self.assertEqual(concatted.name, expected_name)

    with self.subTest('to_mmcif'):
      _ = concatted.to_mmcif()

    cif = concatted.to_mmcif_dict()
    with self.subTest('entities'):
      self.assertEqual(cif['_entity.id'], ['1', '2', '3', '4', '5', '6', '7'])
      self.assertEqual(
          cif['_struct_asym.id'], ['A', 'B', 'C', 'D', 'E', 'F', 'G']
      )
      self.assertEqual(
          cif['_struct_asym.entity_id'], ['1', '2', '3', '4', '5', '6', '7']
      )

    with self.subTest('author_chain_ids'):
      # Author chains should also be renamed and match internal chain IDs.
      self.assertEqual(
          cif['_atom_site.label_asym_id'], cif['_atom_site.auth_asym_id']
      )

  def test_enumerate_atoms(self):
    struc = self.make_test_structure(
        chain_id=np.array(['A', 'A', 'B', 'B', 'B'], dtype=object),
        res_id=np.array([6, 7, 1, 5, 5], dtype=np.int32),
    )
    self.assertSequenceEqual(
        [0, 1, 2, 3, 3],
        [resi for resi, _ in structure.enumerate_residues(struc.iter_atoms())],
    )

  def test_enumerate_atoms_with_all_residues(self):
    all_residues = {
        'A': [('UNK', rid) for rid in range(1, 7)],
        'B': [('UNK', rid) for rid in range(1, 8)],
        'C': [('UNK', rid) for rid in range(1, 3)],
        'D': [('UNK', rid) for rid in range(1, 3)],
    }
    struc = self.make_test_structure(
        chain_id=np.array(
            ['A', 'A', 'A', 'B', 'B', 'B', 'C', 'C', 'C', 'D'], dtype=object
        ),
        res_id=np.array([2, 3, 5, 5, 6, 7, 1, 1, 2, 1], dtype=np.int32),
        res_name=np.array(['UNK'] * 10, dtype=object),
        all_residues=all_residues,
    )
    struc = struc.filter(
        struc.chain_id != 'D',
        cascade_delete=structure.CascadeDelete.NONE,
    )
    res_idxs = [
        resi
        for resi, _ in structure.enumerate_residues(
            struc.iter_atoms(), all_residues=all_residues
        )
    ]
    self.assertSequenceEqual(res_idxs, [1, 2, 4, 10, 11, 12, 13, 13, 14])

  @parameterized.parameters(
      dict(all_res={'A': [1, 2, 3, 5]}),
      dict(all_res={'A': [1, 2, 3, 4]}),
  )
  def test_enumerate_atoms_with_bad_all_residues(self, all_res):
    all_residues = {
        c: [('', i) for i in res_ids] for c, res_ids in all_res.items()
    }
    with self.assertRaises(Exception):
      struc = self.make_test_structure(
          chain_id=np.array(['A', 'A', 'A', 'A'], dtype=object),
          res_id=np.array([2, 2, 4, 5], dtype=np.int32),
          res_name=np.array(['ALA', 'ALA', 'GLY', 'HIS'], dtype=object),
          all_residues=all_residues,
      )
      _ = list(
          structure.enumerate_residues(
              struc.iter_atoms(), all_residues=all_residues
          )
      )

  def test_empty_structure(self):
    struc = structure.make_empty_structure()
    self.assertEqual(struc.num_atoms, 0)
    self.assertEqual(struc.num_residues(count_unresolved=False), 0)
    self.assertEqual(struc.num_chains, 0)
    self.assertEmpty(list(struc.iter_atoms()))
    self.assertEmpty(list(struc.iter_residues()))
    self.assertEmpty(list(struc.iter_chains()))
    self.assertEqual(struc.order_and_drop_atoms_to_match(struc).num_atoms, 0)
    self.assertEmpty(structure.multichain_residue_index(struc))
    self.assertEmpty(struc.is_ligand_mask)

  @parameterized.parameters([
      (
          [
              mmcif_names.DNA_CHAIN,
              mmcif_names.RNA_CHAIN,
              mmcif_names.PROTEIN_CHAIN,
              mmcif_names.PROTEIN_CHAIN,
          ],
          [False, False, True, True],
      ),
      (
          [
              mmcif_names.DNA_CHAIN,
              mmcif_names.RNA_CHAIN,
              mmcif_names.NON_POLYMER_CHAIN,
              mmcif_names.PROTEIN_CHAIN,
          ],
          [False, False, False, True],
      ),
      (
          [
              mmcif_names.NON_POLYMER_CHAIN,
              mmcif_names.BRANCHED_CHAIN,
              mmcif_names.MACROLIDE_CHAIN,
              mmcif_names.NON_POLYMER_CHAIN,
          ],
          [False, False, False, False],
      ),
  ])
  def test_is_protein_mask(self, chain_types, expected_mask):
    struc = self.make_test_structure(
        res_id=np.array([1, 1, 1, 1], dtype=np.int32),
        chain_id=np.array(['A', 'B', 'C', 'D'], dtype=object),
        chain_type=np.array(chain_types, dtype=object),
    )
    self.assertSequenceEqual(struc.is_protein_mask.tolist(), expected_mask)

  @parameterized.parameters([
      (
          [
              mmcif_names.DNA_CHAIN,
              mmcif_names.RNA_CHAIN,
              mmcif_names.PROTEIN_CHAIN,
              mmcif_names.PROTEIN_CHAIN,
          ],
          [False, False, False, False],
      ),
      (
          [
              mmcif_names.DNA_CHAIN,
              mmcif_names.RNA_CHAIN,
              mmcif_names.NON_POLYMER_CHAIN,
              mmcif_names.PROTEIN_CHAIN,
          ],
          [False, False, True, False],
      ),
      (
          [
              mmcif_names.NON_POLYMER_CHAIN,
              mmcif_names.BRANCHED_CHAIN,
              mmcif_names.MACROLIDE_CHAIN,
              mmcif_names.NON_POLYMER_CHAIN,
          ],
          [True, True, True, True],
      ),
  ])
  def test_is_ligand_mask(self, chain_types, expected_mask):
    struc = self.make_test_structure(
        res_id=np.array([1, 1, 1, 1], dtype=np.int32),
        chain_id=np.array(['A', 'B', 'C', 'D'], dtype=object),
        chain_type=np.array(chain_types, dtype=object),
    )
    self.assertSequenceEqual(struc.is_ligand_mask.tolist(), expected_mask)

  @parameterized.parameters([
      (
          [
              mmcif_names.DNA_CHAIN,
              mmcif_names.RNA_CHAIN,
              mmcif_names.PROTEIN_CHAIN,
              mmcif_names.PROTEIN_CHAIN,
          ],
          [True, False, False, False],
      ),
      (
          [
              mmcif_names.DNA_CHAIN,
              mmcif_names.RNA_CHAIN,
              mmcif_names.NON_POLYMER_CHAIN,
              mmcif_names.PROTEIN_CHAIN,
          ],
          [True, False, False, False],
      ),
      (
          [
              mmcif_names.NON_POLYMER_CHAIN,
              mmcif_names.BRANCHED_CHAIN,
              mmcif_names.MACROLIDE_CHAIN,
              mmcif_names.NON_POLYMER_CHAIN,
          ],
          [False, False, False, False],
      ),
  ])
  def test_is_dna_mask(self, chain_types, expected_mask):
    struc = self.make_test_structure(
        res_id=np.array([1, 1, 1, 1], dtype=np.int32),
        chain_id=np.array(['A', 'B', 'C', 'D'], dtype=object),
        chain_type=np.array(chain_types, dtype=object),
    )
    self.assertSequenceEqual(struc.is_dna_mask.tolist(), expected_mask)

  @parameterized.parameters([
      (
          [
              mmcif_names.DNA_CHAIN,
              mmcif_names.RNA_CHAIN,
              mmcif_names.PROTEIN_CHAIN,
              mmcif_names.PROTEIN_CHAIN,
          ],
          [False, True, False, False],
      ),
      (
          [
              mmcif_names.DNA_CHAIN,
              mmcif_names.RNA_CHAIN,
              mmcif_names.NON_POLYMER_CHAIN,
              mmcif_names.PROTEIN_CHAIN,
          ],
          [False, True, False, False],
      ),
      (
          [
              mmcif_names.NON_POLYMER_CHAIN,
              mmcif_names.BRANCHED_CHAIN,
              mmcif_names.MACROLIDE_CHAIN,
              mmcif_names.NON_POLYMER_CHAIN,
          ],
          [False, False, False, False],
      ),
  ])
  def test_is_rna_mask(self, chain_types, expected_mask):
    struc = self.make_test_structure(
        res_id=np.array([1, 1, 1, 1], dtype=np.int32),
        chain_id=np.array(['A', 'B', 'C', 'D'], dtype=object),
        chain_type=np.array(chain_types, dtype=object),
    )
    self.assertSequenceEqual(struc.is_rna_mask.tolist(), expected_mask)

  @parameterized.parameters([
      (
          [
              mmcif_names.DNA_CHAIN,
              mmcif_names.RNA_CHAIN,
              mmcif_names.PROTEIN_CHAIN,
              mmcif_names.PROTEIN_CHAIN,
          ],
          [True, True, False, False],
      ),
      (
          [
              mmcif_names.DNA_CHAIN,
              mmcif_names.RNA_CHAIN,
              mmcif_names.NON_POLYMER_CHAIN,
              mmcif_names.PROTEIN_CHAIN,
          ],
          [True, True, False, False],
      ),
      (
          [
              mmcif_names.NON_POLYMER_CHAIN,
              mmcif_names.BRANCHED_CHAIN,
              mmcif_names.MACROLIDE_CHAIN,
              mmcif_names.NON_POLYMER_CHAIN,
          ],
          [False, False, False, False],
      ),
  ])
  def test_is_nucleic_mask(self, chain_types, expected_mask):
    struc = self.make_test_structure(
        res_id=np.array([1, 1, 1, 1], dtype=np.int32),
        chain_id=np.array(['A', 'B', 'C', 'D'], dtype=object),
        chain_type=np.array(chain_types, dtype=object),
    )
    self.assertSequenceEqual(struc.is_nucleic_mask.tolist(), expected_mask)

  @parameterized.parameters([
      (
          [
              mmcif_names.WATER,
              mmcif_names.WATER,
              mmcif_names.MACROLIDE_CHAIN,
              mmcif_names.BRANCHED_CHAIN,
          ],
          [True, True, False, False],
      ),
      (
          [
              mmcif_names.DNA_CHAIN,
              mmcif_names.RNA_CHAIN,
              mmcif_names.PROTEIN_CHAIN,
              mmcif_names.WATER,
          ],
          [False, False, False, True],
      ),
      (
          [
              mmcif_names.NON_POLYMER_CHAIN,
              mmcif_names.DNA_CHAIN,
              mmcif_names.RNA_CHAIN,
              mmcif_names.PROTEIN_CHAIN,
          ],
          [False, False, False, False],
      ),
  ])
  def test_is_water_mask(self, chain_types, expected_mask):
    struc = self.make_test_structure(
        res_id=np.array([1, 1, 1, 1], dtype=np.int32),
        chain_id=np.array(['A', 'B', 'C', 'D'], dtype=object),
        chain_type=np.array(chain_types, dtype=object),
    )
    self.assertSequenceEqual(struc.is_water_mask.tolist(), expected_mask)

  def test_masks_with_different_chain_ordering(self):
    chain_type = [
        mmcif_names.PROTEIN_CHAIN,
        mmcif_names.DNA_CHAIN,
        mmcif_names.RNA_CHAIN,
        mmcif_names.NON_POLYMER_CHAIN,
        mmcif_names.WATER,
        mmcif_names.WATER,
    ]
    chain_id = np.array(['A', 'B', 'C', 'D', 'E', 'F'], dtype=object)
    struc = self.make_test_structure(
        res_id=np.array([1, 1, 1, 1, 1, 1], dtype=np.int32),
        chain_id=chain_id,
        chain_type=np.array(chain_type, dtype=object),
    )
    perm = np.array([1, 0, 3, 2, 5, 4], dtype=np.int64)
    permuted = struc.reorder_chains(chain_id[perm])  # pyrefly: ignore[bad-argument-type]

    protein_mask = np.array([True, False, False, False, False, False])
    dna_mask = np.array([False, True, False, False, False, False])
    rna_mask = np.array([False, False, True, False, False, False])
    ligand_mask = np.array([False, False, False, True, False, False])
    water_mask = np.array([False, False, False, False, True, True])
    nucleic_mask = np.array([False, True, True, False, False, False])

    np.testing.assert_array_equal(struc.is_protein_mask, protein_mask)
    np.testing.assert_array_equal(struc.is_dna_mask, dna_mask)
    np.testing.assert_array_equal(struc.is_rna_mask, rna_mask)
    np.testing.assert_array_equal(struc.is_ligand_mask, ligand_mask)
    np.testing.assert_array_equal(struc.is_water_mask, water_mask)
    np.testing.assert_array_equal(struc.is_nucleic_mask, nucleic_mask)

    np.testing.assert_array_equal(permuted.is_protein_mask, protein_mask[perm])
    np.testing.assert_array_equal(permuted.is_dna_mask, dna_mask[perm])
    np.testing.assert_array_equal(permuted.is_rna_mask, rna_mask[perm])
    np.testing.assert_array_equal(permuted.is_ligand_mask, ligand_mask[perm])
    np.testing.assert_array_equal(permuted.is_water_mask, water_mask[perm])
    np.testing.assert_array_equal(permuted.is_nucleic_mask, nucleic_mask[perm])

  @parameterized.named_parameters(
      dict(
          testcase_name='mixed_polymers',
          chain_type=[
              *[mmcif_names.PROTEIN_CHAIN] * 2,
              *[mmcif_names.RNA_CHAIN] * 2,
              *[mmcif_names.DNA_CHAIN] * 2,
              *[mmcif_names.NON_POLYMER_CHAIN] * 2,
          ],
          expected_atom_name=['CA', 'OXT', 'C1', 'OXT', 'OP3'],
      ),
      dict(
          testcase_name='protein_without_oxt_and_ligand_oxt',
          chain_type=[
              *[mmcif_names.NON_POLYMER_CHAIN] * 4,
              *[mmcif_names.PROTEIN_CHAIN] * 2,
              *[mmcif_names.NON_POLYMER_CHAIN] * 2,
          ],
          expected_atom_name=[
              'CA',
              'OXT',
              'OP3',
              'OXT',
              'OP3',
              'C1',
              'OXT',
              'OP3',
          ],
      ),
      dict(
          testcase_name='ligands_only',
          chain_type=[mmcif_names.NON_POLYMER_CHAIN] * 8,
          expected_atom_name=[
              'CA',
              'OXT',
              'OP3',
              'OXT',
              'OP3',
              'C1',
              'OXT',
              'OP3',
          ],
      ),
  )
  def test_without_terminal_oxygens(self, chain_type, expected_atom_name):
    struc = self.make_test_structure(
        chain_id=np.array(list('AABBCCDD'), dtype=object),
        res_id=np.ones(8, dtype=np.int32),
        chain_type=np.array(chain_type, dtype=object),
        atom_name=np.array(
            ['CA', 'OXT', 'OP3', 'OXT', 'OP3', 'C1', 'OXT', 'OP3'], dtype=object
        ),
    )
    np.testing.assert_array_equal(
        struc.without_terminal_oxygens().atom_name,
        np.array(expected_atom_name, dtype=object),
        strict=True,
    )

  @parameterized.parameters([True, False])
  def test_drop_non_standard_protein_atoms(self, drop_oxt):
    chain_type = [
        *[mmcif_names.PROTEIN_CHAIN] * 3,
        *[mmcif_names.PROTEIN_CHAIN] * 4,
        *[mmcif_names.NON_POLYMER_CHAIN] * 1,
        *[mmcif_names.RNA_CHAIN] * 5,
    ]
    atom_name = [
        *['N', 'X', 'OXT'],
        *['n', 'CA', 'OXT', 'Y'],
        'Y',
        *['Z1', 'Z2', 'Z3', 'Z4', 'Z5'],
    ]

    struc = self.make_test_structure(
        chain_id=np.array(list('AAABBBBCDDDDD'), dtype=object),
        res_id=np.ones(13, dtype=np.int32),
        chain_type=np.array(chain_type, dtype=object),
        atom_name=np.array(atom_name, dtype=object),
    )
    filtered = struc.drop_non_standard_protein_atoms(drop_oxt=drop_oxt)

    self.assertEqual(filtered.num_chains, 4)
    self.assertEqual(filtered.num_residues(count_unresolved=False), 4)
    self.assertEqual(filtered.num_atoms, 8 if drop_oxt else 10)

    np.testing.assert_array_equal(
        filtered.atom_name,
        np.array([
            'N',
            *['OXT'] * (not drop_oxt),
            'CA',
            *['OXT'] * (not drop_oxt),
            'Y',
            *['Z1', 'Z2', 'Z3', 'Z4', 'Z5'],
        ]),
    )

  def test_drop_non_standard_atoms(self):
    # pyformat: disable
    atom_names = [
        'N', 'CA', 'CB', 'SG', 'SD', 'CE', 'CZ', 'OH', 'C', 'O', 'N', 'CA', 'C',
        'O', 'CB', 'SG', 'ZN', 'CA', 'C1', 'C100'
    ]
    # pyformat: enable
    struc = self.make_test_structure(
        chain_id=np.array(
            [*['A'] * 10, *['B'] * 6, *['C'] * 2, *['D'] * 2], dtype=object
        ),
        res_id=np.ones(20, dtype=np.int32),
        res_name=np.array(
            [*['CME'] * 10, *['CYS'] * 6, *['ZN'] * 2, *['LWDM:)'] * 2],
            dtype=object,
        ),
        chain_type=np.array(
            ['polypeptide(L)'] * 16 + ['non-polymer'] * 4, dtype=object
        ),
        atom_name=np.array(atom_names, dtype=object),
    )

    filtered = struc.fix_non_standard_polymer_res()  # CME -> CYS.
    filtered = filtered.drop_non_standard_atoms(
        ccd=chemical_components.Ccd(), drop_unk=False, drop_non_ccd=False
    )

    self.assertEqual(filtered.num_chains, 4)
    self.assertEqual(filtered.num_residues(count_unresolved=False), 4)
    self.assertEqual(filtered.num_atoms, 15)

    np.testing.assert_array_equal(
        filtered.atom_name,
        np.array(
            [
                *('N', 'CA', 'CB', 'SG', 'C', 'O'),  # CYS.
                *('N', 'CA', 'C', 'O', 'CB', 'SG'),  # CME -> CYS.
                'ZN',  # Ligand kept, CA not.
                *('C1', 'C100'),  # Keep non-ccd.
            ],
            dtype=object,
        ),
    )

  def test_chain_boundaries_empty_structure(self):
    struc = structure_lib.make_empty_structure()
    struc = struc.copy_and_update(
        chains=structure.Chains.from_defaults(id=np.array(['A'], dtype=object)),
    )
    boundaries = struc.chain_boundaries
    self.assertEmpty(boundaries)

  def test_chain_boundaries_empty_chain_sandwich(self):
    struc = structure.Structure(
        chains=structure.Chains.from_defaults(
            key=np.array([0, 1, 2], dtype=np.int64),
            id=np.array(['A', 'B', 'C'], dtype=object),
        ),
        residues=structure.Residues.from_defaults(
            key=np.array([0, 1, 2], dtype=np.int64),
            chain_key=np.array([0, 1, 2], dtype=np.int64),
            id=np.array([1, 1, 1], dtype=np.int32),
            name=np.array(['ALA', 'UNK', 'PRO'], dtype=object),
        ),
        atoms=structure.Atoms.from_defaults(
            chain_key=np.array([0, 2], dtype=np.int64),
            res_key=np.array([0, 2], dtype=np.int64),
            name=np.array(['CA', 'CA'], dtype=object),
        ),
        bonds=structure.Bonds.make_empty(),
    )

    boundaries = struc.chain_boundaries
    np.testing.assert_array_equal(boundaries, np.array([0, 1]))

  def test_chain_boundaries_one_chain(self):
    struc = self.make_test_structure(
        chain_id=np.array(['A', 'A', 'A'], dtype=object),
        res_id=np.array([1, 2, 3], dtype=np.int32),
        res_name=np.array(['CYS', 'GLY', 'ALA'], dtype=object),
    )

    boundaries = struc.chain_boundaries
    self.assertLen(boundaries, 1)
    np.testing.assert_array_equal(np.array([0]), boundaries)

  def test_chain_boundaries_two_chains(self):
    struc = self.make_test_structure(
        chain_id=np.array(['A', 'A', 'A', 'B'], dtype=object),
        res_id=np.array([1, 2, 3, 4], dtype=np.int32),
        res_name=np.array(['CYS', 'GLY', 'ALA', 'HIS'], dtype=object),
    )

    boundaries = struc.chain_boundaries
    np.testing.assert_array_equal(np.array([0, 3]), boundaries)

  def test_find_chains_with_unknown_sequence(self):
    struc = self.make_test_structure(
        chain_id=np.array(['A', 'A', 'B', 'C', 'D', 'E', 'E'], dtype=object),
        res_id=np.array([1, 2, 3, 4, 5, 6, 7], dtype=np.int32),
        res_name=np.array(
            ['CYS', 'UNK', 'UNK', 'G', 'N', 'UNK', 'UNL'], dtype=object
        ),
        chain_type=np.array(
            ['polypeptide(L)'] * 3
            + ['polyribonucleotide'] * 2
            + ['non-polymer'] * 2,
            dtype=object,
        ),
        atom_name=np.array(
            ['CA', 'CA', 'CA', 'P', 'P', 'C', 'C'], dtype=object
        ),
    )
    expected_removed = 3
    expected_chains = ['A', 'C']

    chains_with_unknown_sequences = struc.find_chains_with_unknown_sequence()
    n_removed = len(chains_with_unknown_sequences)
    filtered = struc.filter_out(chain_id=chains_with_unknown_sequences)

    with self.subTest('partly_resolved_chains_not_returned'):
      self.assertNotIn('A', chains_with_unknown_sequences)
    with self.subTest('fully_resolved_chains_not_returned'):
      self.assertNotIn('C', chains_with_unknown_sequences)
    with self.subTest('mixed_unknown_types_returned'):
      self.assertIn('E', chains_with_unknown_sequences)
    self.assertEqual(n_removed, expected_removed)
    np.testing.assert_array_equal(
        np.array(filtered.chains), np.array(expected_chains)
    )

  def test_find_chains_with_unknown_sequence_finds_unresolved_chains(self):
    struc = structure.Structure(
        chains=structure.Chains.from_defaults(
            key=np.array([0, 1, 2, 3], dtype=np.int64),
            id=np.array(['A', 'B', 'C', 'D'], dtype=object),
        ),
        residues=structure.Residues.from_defaults(
            key=np.array([0, 1, 2, 3, 4, 5, 6], dtype=np.int64),
            chain_key=np.array([0, 0, 0, 1, 2, 2, 3], dtype=np.int64),
            id=np.array([1, 2, 3, 1, 1, 2, 1], dtype=np.int32),
            name=np.array(
                ['ALA', 'ILE', 'GLY', 'UNK', 'UNK', 'PRO', 'ALA'],
                dtype=object,
            ),
        ),
        atoms=structure.Atoms.from_defaults(
            chain_key=np.array([0, 0, 2], dtype=np.int64),
            res_key=np.array([0, 2, 4], dtype=np.int64),
            name=np.array(['CA', 'CA', 'CA'], dtype=object),
        ),
        bonds=structure.Bonds.make_empty(),
    )

    chains_with_unknown_sequences = struc.find_chains_with_unknown_sequence()

    self.assertEqual(chains_with_unknown_sequences, ['B'])

  def test_find_chains_with_unknown_sequence_with_empty_structure(self):
    struc = structure.make_empty_structure()
    struc = struc.copy_and_update(
        chains=structure.Chains.from_defaults(id=np.array(['A'], dtype=object)),
    )

    chains_with_unknown_sequences = struc.find_chains_with_unknown_sequence()

    self.assertEqual(chains_with_unknown_sequences, ['A'])

  def test_resolved_residues_mask(self):
    struc = structure.from_mmcif(
        test_utils.read_test_cif_file('2pv7'), include_bonds=False
    )
    self.assertEqual(
        struc.num_residues(count_unresolved=True),
        len(struc.present_residues) + len(struc.unresolved_residues),
    )

    # Check the test isn't trivially passing.
    self.assertGreater(struc.num_residues(count_unresolved=True), 0)
    self.assertGreater(struc.num_residues(count_unresolved=False), 0)

    self.assertLen(
        struc.present_residues,
        struc.num_residues(count_unresolved=False),
    )

    self.assertLen(
        struc.unresolved_residues,
        struc.num_residues(count_unresolved=True)
        - struc.num_residues(count_unresolved=False),
    )

  def test_spaces_in_name(self):
    name = 'name with spaces'
    struc = structure.from_sequences_and_bonds(
        name=name,
        sequences=['A'],
        chain_types=[mmcif_names.PROTEIN_CHAIN],
        sequence_formats=[structure.SequenceFormat.FASTA],
        bonded_atom_pairs=None,
        ccd=chemical_components.Ccd(),
    )
    mmcif = struc.to_mmcif()
    self.assertEqual(structure.from_mmcif(mmcif).name, name)

    cif = cif_dict.from_string(mmcif)
    self.assertEqual(cif['data_'], ['name-with-spaces'])
    self.assertEqual(cif['_entry.id'], ['name with spaces'])

  def test_fully_unresolved_chain_mmcif_roundtrip(self):
    struc = self.make_test_structure(
        chain_id=np.array(['A', 'A', 'B', 'B'], dtype=object),
        res_id=np.array([1, 2, 3, 4], dtype=np.int32),
        res_name=np.array(['CYS', 'GLU', 'ARG', 'PRO'], dtype=object),
        chain_type=np.array([mmcif_names.PROTEIN_CHAIN] * 4, dtype=object),
        atom_name=np.array(['CA', 'CA', 'CA', 'CA'], dtype=object),
    )

    # Make all residues in chain A unresolved.
    struc = struc.filter_out(
        res_id=[1, 2], cascade_delete=structure.CascadeDelete.NONE
    )
    self.assertEqual(struc.num_chains, 2)
    self.assertEqual(struc.num_residues(count_unresolved=True), 4)
    self.assertEqual(struc.num_residues(count_unresolved=False), 2)
    self.assertEqual(struc.num_atoms, 2)
    struc = structure.from_mmcif(struc.to_mmcif())
    self.assertEqual(struc.num_chains, 2)
    self.assertEqual(struc.num_residues(count_unresolved=True), 4)
    self.assertEqual(struc.num_residues(count_unresolved=False), 2)
    self.assertEqual(struc.num_atoms, 2)

  def test_order_and_drop_atoms_to_match_identical_layout_returns_self(self):
    struc = self.make_test_structure(
        chain_id=np.array(['A', 'B'], dtype=object),
        res_id=np.array([1, 2], dtype=np.int32),
        atom_name=np.array(['CA', 'CB'], dtype=object),
    )
    shifted = self.make_test_structure(
        chain_id=np.array(['A', 'B'], dtype=object),
        res_id=np.array([1, 2], dtype=np.int32),
        atom_name=np.array(['CA', 'CB'], dtype=object),
        atom_x=np.array([10.0, 20.0]),
    )
    self.assertIs(struc.order_and_drop_atoms_to_match(shifted), struc)

  def test_order_and_drop_atoms_to_match_reordered_atoms(self):
    struc = self.make_test_structure(
        chain_id=np.array(['A', 'B'], dtype=object),
        res_id=np.array([1, 2], dtype=np.int32),
        atom_name=np.array(['CA', 'CB'], dtype=object),
    )
    reversed_struc = struc.reorder_chains(['B', 'A'])
    reordered = struc.order_and_drop_atoms_to_match(reversed_struc)
    self.assertIsNot(reordered, struc)
    np.testing.assert_array_equal(reordered.atom_name, reversed_struc.atom_name)

  def test_order_and_drop_atoms_to_match_shifted_and_gapped_chain_keys(self):
    base = self.make_test_structure(
        chain_id=np.array(['A', 'B', 'C'], dtype=object),
        res_id=np.array([1, 1, 1], dtype=np.int32),
        res_name=np.array(['ALA', 'ALA', 'ALA'], dtype=object),
        atom_name=np.array(['CA', 'CA', 'CA'], dtype=object),
        atom_x=np.array([0.0, 5.0, 10.0], dtype=np.float32),
        atom_y=np.zeros(3, dtype=np.float32),
        atom_z=np.zeros(3, dtype=np.float32),
    )
    ab = base.filter(chain_id=('A', 'B'))
    bc_as_ab = base.filter(chain_id=('B', 'C')).rename_chain_ids(
        {'B': 'A', 'C': 'B'}
    )
    ac_as_ab = base.filter(chain_id=('A', 'C')).rename_chain_ids({'C': 'B'})
    matched_contiguous = bc_as_ab.order_and_drop_atoms_to_match(ab)
    np.testing.assert_array_equal(matched_contiguous.coords, bc_as_ab.coords)
    matched_gapped = ac_as_ab.order_and_drop_atoms_to_match(ab)
    np.testing.assert_array_equal(matched_gapped.coords, ac_as_ab.coords)

  def test_filter_by_chains_noop_returns_self(self):
    mono = self.make_test_structure(
        chain_id=np.array(['A', 'A'], dtype=object),
        res_id=np.array([1, 2], dtype=np.int32),
        res_name=np.array(['ALA', 'GLY'], dtype=object),
        atom_name=np.array(['CA', 'CA'], dtype=object),
    )
    self.assertIs(mono.filter(chain_id='A'), mono)
    self.assertIs(
        mono.filter(chain_id='A', cascade_delete=structure.CascadeDelete.FULL),
        mono,
    )

  def test_split_by_chain_empty_and_single_chain_fast_path(self):
    empty_struc = structure.Structure(
        chains=structure.Chains.make_empty(),
        residues=structure.Residues.make_empty(),
        atoms=structure.Atoms.make_empty(),
        bonds=structure.Bonds.make_empty(),
    )
    self.assertEqual(empty_struc.split_by_chain(), [])

    mono = self.make_test_structure(
        chain_id=np.array(['A', 'A'], dtype=object),
        res_id=np.array([1, 2], dtype=np.int32),
        res_name=np.array(['ALA', 'GLY'], dtype=object),
        atom_name=np.array(['CA', 'CA'], dtype=object),
    )
    split = mono.split_by_chain()
    self.assertLen(split, 1)
    self.assertIs(split[0], mono)

    empty_chain_struc = structure.Structure(
        chains=structure.Chains.from_defaults(
            key=np.array([0], dtype=np.int64),
            id=np.array(['A'], dtype=object),
        ),
        residues=structure.Residues.from_defaults(
            key=np.array([0], dtype=np.int64),
            chain_key=np.array([0], dtype=np.int64),
            id=np.array([1], dtype=np.int32),
            name=np.array(['ALA'], dtype=object),
        ),
        atoms=structure.Atoms.make_empty(),
        bonds=structure.Bonds.make_empty(),
    )
    empty_split = empty_chain_struc.split_by_chain()
    self.assertLen(empty_split, 1)
    self.assertEqual(empty_split[0].num_chains, 0)
    self.assertEqual(empty_split[0].num_atoms, 0)

  def test_split_by_chain_intra_and_inter_chain_bonds(self):
    with_bonds = structure.Structure(
        chains=structure.Chains.from_defaults(
            key=np.array([10, 20, 30], dtype=np.int64),
            id=np.array(['A', 'B', 'C'], dtype=object),
        ),
        residues=structure.Residues.from_defaults(
            key=np.array([100, 101, 200, 300], dtype=np.int64),
            chain_key=np.array([10, 20, 20, 30], dtype=np.int64),
            id=np.array([1, 1, 2, 1], dtype=np.int32),
        ),
        atoms=structure.Atoms.from_defaults(
            key=np.array([1000, 1001, 2000, 2001], dtype=np.int64),
            chain_key=np.array([20, 20, 30, 30], dtype=np.int64),
            res_key=np.array([200, 200, 300, 300], dtype=np.int64),
        ),
        bonds=structure.Bonds(
            key=np.array([0, 1, 2], dtype=np.int64),
            type=np.array(['covale', 'covale', 'covale'], dtype=object),
            role=np.array(['?', '?', '?'], dtype=object),
            from_atom_key=np.array([2000, 1000, 1001], dtype=np.int64),
            dest_atom_key=np.array([2001, 2000, 1000], dtype=np.int64),
        ),
    )
    split_chains = with_bonds.split_by_chain()
    self.assertLen(split_chains, 3)
    for part, cid in zip(split_chains, ('A', 'B', 'C'), strict=True):
      self.assertStructuresEqual(part, with_bonds.filter(chain_id=cid))

  def test_split_by_chain_unsorted_keys(self):
    chain_keys = np.array([2, 1], dtype=np.int64)
    res_keys = np.array([0, 1], dtype=np.int64)
    unsorted_keys = structure.Structure(
        chains=structure.Chains.from_defaults(
            key=chain_keys, id=np.array(['B', 'A'], dtype=object)
        ),
        residues=structure.Residues.from_defaults(
            chain_key=chain_keys, id=np.array([1, 1], dtype=np.int32)
        ),
        atoms=structure.Atoms.from_defaults(
            chain_key=chain_keys, res_key=res_keys
        ),
        bonds=structure.Bonds.make_empty(),
    )
    split_chains = unsorted_keys.split_by_chain()
    self.assertLen(split_chains, 2)
    for part, cid in zip(split_chains, ('B', 'A'), strict=True):
      self.assertStructuresEqual(part, unsorted_keys.filter(chain_id=cid))

  @parameterized.named_parameters(
      dict(testcase_name='adjacent_duplicates', chain_ids=('A', 'A')),
      dict(testcase_name='non_adjacent_duplicates', chain_ids=('A', 'B', 'A')),
  )
  def test_split_by_chain_duplicate_chain_ids(self, chain_ids):
    n = len(chain_ids)
    keys = np.arange(n, dtype=np.int64)
    struc = structure.Structure(
        chains=structure.Chains.from_defaults(
            id=np.array(chain_ids, dtype=object)
        ),
        residues=structure.Residues.from_defaults(
            chain_key=keys, id=np.arange(1, n + 1, dtype=np.int32)
        ),
        atoms=structure.Atoms.from_defaults(chain_key=keys, res_key=keys),
        bonds=structure.Bonds.make_empty(),
    )
    self.assertLen(struc.split_by_chain(), n)

  def test_split_by_chain_multiple_empty_and_multistate_chains(self):
    coords = np.arange(12, dtype=np.float32).reshape(2, 2, 3)
    struc = structure.Structure(
        chains=structure.Chains.from_defaults(
            id=np.array(['A', 'B', 'C'], dtype=object)
        ),
        residues=structure.Residues.from_defaults(
            chain_key=np.array([0, 1], dtype=np.int64),
            id=np.array([1, 1], dtype=np.int32),
        ),
        atoms=structure.Atoms.from_defaults(
            chain_key=np.array([1, 1], dtype=np.int64),
            res_key=np.array([1, 1], dtype=np.int64),
            x=coords[..., 0],
            y=coords[..., 1],
            z=coords[..., 2],
        ),
        bonds=structure.Bonds.make_empty(),
    )
    split = struc.split_by_chain()
    self.assertLen(split, 3)
    self.assertStructuresEqual(split[0], struc.filter(chain_id='A'))
    self.assertStructuresEqual(split[2], struc.filter(chain_id='C'))
    self.assertEqual(split[0].coords.shape, (2, 0, 3))
    self.assertEqual(split[1].coords.shape, (2, 2, 3))


if __name__ == '__main__':
  absltest.main()
