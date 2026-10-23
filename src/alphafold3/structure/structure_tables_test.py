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
from alphafold3.structure import structure_tables
import numpy as np


class AtomsTest(parameterized.TestCase):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self._atoms = structure_tables.Atoms(
        key=np.array([1, 2, 3]),
        chain_key=np.array([1, 1, 2]),
        res_key=np.array([1, 2, 1]),
        name=np.array(['CA', 'N', 'O'], dtype=object),
        element=np.array(['C', 'N', 'O'], dtype=object),
        x=np.array(
            [
                [0.1, 0.2, 0.3],
                [0.4, 0.5, 0.6],
            ],
            dtype=np.float32,
        ),
        y=np.array(
            [
                [1.1, 1.2, 1.3],
                [1.4, 1.5, 1.6],
            ],
            dtype=np.float32,
        ),
        z=np.array(
            [
                [2.1, 2.2, 2.3],
                [2.4, 2.5, 2.6],
            ],
            dtype=np.float32,
        ),
        b_factor=np.array([0.1, 0.2, 0.3], dtype=np.float32),
        occupancy=np.array([1.1, 1.2, 1.3], dtype=np.float32),
    )

  def test_make_empty(self):
    empty_atoms = structure_tables.Atoms.make_empty()
    self.assertEqual(empty_atoms.size, 0)

  def test_from_defaults(self):
    atoms = structure_tables.Atoms.from_defaults(
        chain_key=np.array([0, 0, 1], np.int64),
        res_key=np.array([0, 1, 2], np.int64),
    )
    np.testing.assert_array_equal(atoms.name, np.full(3, '?', dtype=object))
    np.testing.assert_array_equal(atoms.element, np.full(3, '?', dtype=object))
    np.testing.assert_array_equal(atoms.x, np.zeros(3, dtype=np.float32))
    np.testing.assert_array_equal(atoms.y, np.zeros(3, dtype=np.float32))
    np.testing.assert_array_equal(atoms.z, np.zeros(3, dtype=np.float32))
    np.testing.assert_array_equal(atoms.b_factor, np.zeros(3, dtype=np.float32))
    np.testing.assert_array_equal(atoms.occupancy, np.ones(3, dtype=np.float32))

  def test_from_defaults_is_validated(self):
    with self.assertRaises(ValueError):
      structure_tables.Atoms.from_defaults(
          chain_key=np.array([0, 0, 1], np.int64),
          res_key=np.array([0, 1], np.int64),
      )

  @parameterized.named_parameters(
      ('nan', float('nan')),
      ('positive_infinity', float('+inf')),
      ('negative_infinity', float('-inf')),
  )
  def test_finite_values(self, value):
    with self.assertRaises(ValueError):
      structure_tables.Atoms(
          key=np.array([1]),
          chain_key=np.array([1]),
          res_key=np.array([1]),
          name=np.array(['CA'], dtype=object),
          element=np.array(['C'], dtype=object),
          x=np.array([value], dtype=np.float32),
          y=np.array([1.0], dtype=np.float32),
          z=np.array([1.0], dtype=np.float32),
          b_factor=np.array([1.0], dtype=np.float32),
          occupancy=np.array([1.0], dtype=np.float32),
      )

  @parameterized.named_parameters(
      ('x', {'x': np.array([float('nan')], dtype=np.float32)}),
      ('y', {'y': np.array([float('nan')], dtype=np.float32)}),
      ('z', {'z': np.array([float('nan')], dtype=np.float32)}),
      ('b_factor', {'b_factor': np.array([float('nan')], dtype=np.float32)}),
      ('occupancy', {'occupancy': np.array([float('nan')], dtype=np.float32)}),
  )
  def test_finite_validation(self, kwargs_update):
    kwargs = {
        'key': np.array([1]),
        'chain_key': np.array([1]),
        'res_key': np.array([1]),
        'name': np.array(['CA'], dtype=object),
        'element': np.array(['C'], dtype=object),
        'x': np.array([1.0], dtype=np.float32),
        'y': np.array([1.0], dtype=np.float32),
        'z': np.array([1.0], dtype=np.float32),
        'b_factor': np.array([1.0], dtype=np.float32),
        'occupancy': np.array([1.0], dtype=np.float32),
    }
    kwargs.update(kwargs_update)
    with self.assertRaises(ValueError):
      structure_tables.Atoms(**kwargs)

  def test_get_value_by_index(self):
    with self.subTest('name'):
      self.assertEqual(self._atoms.get_value_by_index('name', 0), 'CA')
      self.assertEqual(self._atoms.get_value_by_index('name', 1), 'N')
      self.assertEqual(self._atoms.get_value_by_index('name', 2), 'O')
    with self.subTest('coords'):
      np.testing.assert_array_equal(
          self._atoms.get_value_by_index('x', 0),
          np.array([0.1, 0.4], dtype=np.float32),
      )
      np.testing.assert_array_equal(
          self._atoms.get_value_by_index('x', 1),
          np.array([0.2, 0.5], dtype=np.float32),
      )
      np.testing.assert_array_equal(
          self._atoms.get_value_by_index('x', 2),
          np.array([0.3, 0.6], dtype=np.float32),
      )
    with self.subTest('occupancy'):
      self.assertEqual(
          self._atoms.get_value_by_index('occupancy', 0), np.float32(1.1)
      )
      self.assertEqual(
          self._atoms.get_value_by_index('occupancy', 1), np.float32(1.2)
      )
      self.assertEqual(
          self._atoms.get_value_by_index('occupancy', 2), np.float32(1.3)
      )

  def test_copy_and_update_coords(self):
    num_atoms = self._atoms.size
    single_dim = self._atoms.copy_and_update_coords(
        np.zeros([num_atoms, 3], dtype=np.float32)
    )
    with self.subTest('1D'):
      np.testing.assert_array_equal(single_dim.x, 0.0)
      np.testing.assert_array_equal(single_dim.y, 0.0)
      np.testing.assert_array_equal(single_dim.z, 0.0)
      self.assertEqual(single_dim.x.shape, (num_atoms,))
      self.assertEqual(single_dim.y.shape, (num_atoms,))
      self.assertEqual(single_dim.z.shape, (num_atoms,))
    two_dim = self._atoms.copy_and_update_coords(
        np.zeros([6, num_atoms, 3], dtype=np.float32)
    )
    with self.subTest('2D'):
      np.testing.assert_array_equal(two_dim.x, 0.0)
      np.testing.assert_array_equal(two_dim.y, 0.0)
      np.testing.assert_array_equal(two_dim.z, 0.0)
      self.assertEqual(two_dim.x.shape, (6, num_atoms))
      self.assertEqual(two_dim.y.shape, (6, num_atoms))
      self.assertEqual(two_dim.z.shape, (6, num_atoms))
    three_dim = self._atoms.copy_and_update_coords(
        np.zeros([4, 6, num_atoms, 3], dtype=np.float32)
    )
    with self.subTest('3D'):
      np.testing.assert_array_equal(three_dim.x, 0.0)
      np.testing.assert_array_equal(three_dim.y, 0.0)
      np.testing.assert_array_equal(three_dim.z, 0.0)
      self.assertEqual(three_dim.x.shape, (4, 6, num_atoms))
      self.assertEqual(three_dim.y.shape, (4, 6, num_atoms))
      self.assertEqual(three_dim.z.shape, (4, 6, num_atoms))

  def test_copy_and_update_coords_invalid_shape(self):
    with self.subTest('final_dim'):
      with self.assertRaisesRegex(ValueError, 'Expecting 3-dimensional'):
        self._atoms.copy_and_update_coords(np.zeros([3, 5]))
    with self.subTest('atom_dim'):
      with self.assertRaisesRegex(ValueError, 'All columns should have length'):
        self._atoms.copy_and_update_coords(np.zeros([self._atoms.size + 1, 3]))


class ResiduesTest(parameterized.TestCase):

  def test_make_empty(self):
    empty_residues = structure_tables.Residues.make_empty()
    self.assertEqual(empty_residues.size, 0)

  def test_from_defaults(self):
    residues = structure_tables.Residues.from_defaults(
        chain_key=np.array([0, 0, 1], np.int64),
        id=np.array([1, 2, 3], np.int32),
    )
    np.testing.assert_array_equal(
        residues.name, np.full(3, 'UNK', dtype=object)
    )
    np.testing.assert_array_equal(
        residues.auth_seq_id, np.array(['1', '2', '3'], dtype=object)
    )
    np.testing.assert_array_equal(
        residues.insertion_code, np.full(3, '?', dtype=object)
    )

  def test_from_defaults_is_validated(self):
    with self.assertRaises(ValueError):
      structure_tables.Residues.from_defaults(
          chain_key=np.array([0, 0, 1], np.int64),
          id=np.array([1, 2], np.int32),
      )


class ChainsTest(parameterized.TestCase):

  def test_make_empty(self):
    empty_chains = structure_tables.Chains.make_empty()
    self.assertEqual(empty_chains.size, 0)

  def test_from_defaults(self):
    chains = structure_tables.Chains.from_defaults(
        id=np.array(['A', 'B', 'C'], dtype=object),
    )
    np.testing.assert_array_equal(
        chains.type, np.full(3, mmcif_names.PROTEIN_CHAIN, dtype=object)
    )
    np.testing.assert_array_equal(
        chains.auth_asym_id, np.array(['A', 'B', 'C'], dtype=object)
    )
    np.testing.assert_array_equal(
        chains.entity_id, np.array(['1', '2', '3'], dtype=object)
    )
    np.testing.assert_array_equal(
        chains.entity_desc, np.array(['.', '.', '.'], dtype=object)
    )

  def test_from_defaults_is_validated(self):
    with self.assertRaises(ValueError):
      structure_tables.Chains.from_defaults(
          id=np.array(['A', 'B'], dtype=object),
          entity_id=np.array(['1', '2', '3'], dtype=object),
      )


class MmcifWritingTest(absltest.TestCase):

  def test_to_mmcif_sequence_and_entity_tables(self):
    chains = structure_tables.Chains(
        key=np.array([1, 2, 3, 4, 5], dtype=np.int64),
        id=np.array(['A', 'B', 'C', 'D', 'E'], dtype=object),
        type=np.array(
            [
                mmcif_names.PROTEIN_CHAIN,
                mmcif_names.RNA_CHAIN,
                mmcif_names.PROTEIN_CHAIN,
                mmcif_names.BRANCHED_CHAIN,
                mmcif_names.NON_POLYMER_CHAIN,
            ],
            dtype=object,
        ),
        auth_asym_id=np.array(['a', 'b', 'c', 'c', 'a'], dtype=object),
        entity_id=np.array([1, 2, 1, 3, 4]).astype(str).astype(object),
        entity_desc=np.array(
            ['ent_1', 'ent_2', 'ent_1', 'ent_3', 'ent_4'], dtype=object
        ),
    )
    residues = structure_tables.Residues(
        key=np.array([1, 2, 3, 4, 5, 6, 7, 8], dtype=np.int64),
        chain_key=np.array([1, 1, 2, 3, 3, 4, 4, 5], dtype=np.int64),
        id=np.array([1, 2, 1, 1, 2, 1, 2, 1], dtype=np.int64),
        name=np.array(
            ['GLY', 'ALA', 'G', 'GLY', 'ALA', 'NAG', 'NAG', 'UNK'], dtype=object
        ),
        auth_seq_id=np.array(
            ['-2', '-2', '0', '-2', '-2', '8', '9', '7'], dtype=object
        ),
        insertion_code=np.array(
            ['A', 'B', None, 'A', 'B', None, None, None], dtype=object
        ),
    )
    # The following array makes the first protein residues "missing".
    atom_res_key = np.array([2, 3, 5, 6, 7, 8], dtype=np.int64)

    cif = structure_tables.to_mmcif_sequence_and_entity_tables(
        chains, residues, atom_res_key
    )

    with self.subTest('_struct_asym'):
      self.assertEqual(cif['_struct_asym.id'], ['A', 'B', 'C', 'D', 'E'])
      self.assertEqual(cif['_struct_asym.entity_id'], ['1', '2', '1', '3', '4'])

    with self.subTest('_entity'):
      self.assertEqual(cif['_entity.id'], ['1', '2', '3', '4'])
      self.assertEqual(
          cif['_entity.pdbx_description'], ['ent_1', 'ent_2', 'ent_3', 'ent_4']
      )
      self.assertEqual(
          cif['_entity.type'],
          [
              mmcif_names.POLYMER_CHAIN,
              mmcif_names.POLYMER_CHAIN,
              mmcif_names.BRANCHED_CHAIN,
              mmcif_names.NON_POLYMER_CHAIN,
          ],
      )

    with self.subTest('_entity_poly'):
      self.assertEqual(cif['_entity_poly.entity_id'], ['1', '2'])
      self.assertEqual(
          cif['_entity_poly.type'],
          [mmcif_names.PROTEIN_CHAIN, mmcif_names.RNA_CHAIN],
      )
      self.assertEqual(cif['_entity_poly.pdbx_strand_id'], ['a,c', 'b'])

    with self.subTest('_entity_poly_seq'):
      self.assertEqual(cif['_entity_poly_seq.entity_id'], ['1', '1', '2'])
      self.assertEqual(cif['_entity_poly_seq.num'], ['1', '2', '1'])
      self.assertEqual(cif['_entity_poly_seq.mon_id'], ['GLY', 'ALA', 'G'])
      self.assertEqual(cif['_entity_poly_seq.hetero'], ['n', 'n', 'n'])

    with self.subTest('_pdbx_poly_seq_scheme'):
      self.assertEqual(
          cif['_pdbx_poly_seq_scheme.asym_id'], ['A', 'A', 'B', 'C', 'C']
      )
      self.assertEqual(
          cif['_pdbx_poly_seq_scheme.entity_id'], ['1', '1', '2', '1', '1']
      )
      self.assertEqual(
          cif['_pdbx_poly_seq_scheme.seq_id'], ['1', '2', '1', '1', '2']
      )
      self.assertEqual(
          cif['_pdbx_poly_seq_scheme.mon_id'], ['GLY', 'ALA', 'G', 'GLY', 'ALA']
      )
      self.assertEqual(
          cif['_pdbx_poly_seq_scheme.pdb_seq_num'],
          ['-2', '-2', '0', '-2', '-2'],
      )
      self.assertEqual(
          cif['_pdbx_poly_seq_scheme.auth_seq_num'], ['?', '-2', '0', '?', '-2']
      )
      self.assertEqual(
          cif['_pdbx_poly_seq_scheme.pdb_strand_id'], ['a', 'a', 'b', 'c', 'c']
      )
      self.assertEqual(
          cif['_pdbx_poly_seq_scheme.pdb_ins_code'], ['A', 'B', '.', 'A', 'B']
      )
      self.assertEqual(
          cif['_pdbx_poly_seq_scheme.hetero'], ['n', 'n', 'n', 'n', 'n']
      )

    with self.subTest('_pdbx_branch_scheme'):
      self.assertEqual(cif['_pdbx_branch_scheme.asym_id'], ['D', 'D'])
      self.assertEqual(cif['_pdbx_branch_scheme.entity_id'], ['3', '3'])
      self.assertEqual(cif['_pdbx_branch_scheme.mon_id'], ['NAG', 'NAG'])
      self.assertEqual(cif['_pdbx_branch_scheme.num'], ['1', '2'])
      self.assertEqual(cif['_pdbx_branch_scheme.pdb_asym_id'], ['c', 'c'])
      self.assertEqual(cif['_pdbx_branch_scheme.pdb_seq_num'], ['8', '9'])
      self.assertEqual(cif['_pdbx_branch_scheme.auth_asym_id'], ['c', 'c'])
      self.assertEqual(cif['_pdbx_branch_scheme.auth_seq_num'], ['8', '9'])
      self.assertEqual(cif['_pdbx_branch_scheme.pdb_ins_code'], ['.', '.'])
      self.assertEqual(cif['_pdbx_branch_scheme.hetero'], ['n', 'n'])

    with self.subTest('_pdbx_nonpoly_scheme'):
      self.assertEqual(cif['_pdbx_nonpoly_scheme.asym_id'], ['E'])
      self.assertEqual(cif['_pdbx_nonpoly_scheme.entity_id'], ['4'])
      self.assertEqual(cif['_pdbx_nonpoly_scheme.mon_id'], ['UNK'])
      self.assertEqual(cif['_pdbx_nonpoly_scheme.pdb_seq_num'], ['7'])
      self.assertEqual(cif['_pdbx_nonpoly_scheme.auth_seq_num'], ['7'])
      self.assertEqual(cif['_pdbx_nonpoly_scheme.pdb_strand_id'], ['a'])
      self.assertEqual(cif['_pdbx_nonpoly_scheme.pdb_ins_code'], ['.'])


if __name__ == '__main__':
  absltest.main()
