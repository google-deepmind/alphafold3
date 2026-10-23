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
from alphafold3 import structure
from alphafold3.common import resources
from alphafold3.common.testing import data
from alphafold3.cpp import cif_dict
from alphafold3.cpp import mmcif_utils
import numpy as np

_TEST_CIFS = (
    '100d', '170d', '1al4', '1bna', '1h9h', '1ign', '1ob2', '1wib', '2pv7',
    '3co3', '3ibw', '3j7l', '3vyp', '4gxy', '5ffz', '5inz', '5jsu', '6bea',
    '3b9f',
)  # pyformat: disable
# '6d4l' Structure does not find deuterium and hydrogen alt-locs.


_DIFFERENT_ALT_LOC_POS = frozenset({'3b9f'})


_TEST_DATA = data.Data(
    resources.ROOT / 'test_data/miniature_databases/pdb_mmcif'
)


def _read_test_cif_file(name: str) -> cif_dict.CifDict:
  return cif_dict.from_string(_TEST_DATA.load(f'{name}.cif'))


def _atom_site_str_array(
    cif: cif_dict.CifDict, column: str, gather
) -> np.ndarray:
  return cif.get_array(f'_atom_site.{column}', dtype=object, gather=gather)


class MmcifToArraysTest(parameterized.TestCase):

  @parameterized.product(
      cif_name=_TEST_CIFS,
      fix_arginines=[False, True],
      keep_ligands_water=[False, True],
      model_id=[structure.ModelID.FIRST, structure.ModelID.ALL],
  )
  def test_atom_site_index_matches_structure(
      self, cif_name, fix_arginines, keep_ligands_water, model_id
  ):
    cif = _read_test_cif_file(cif_name)
    struc = structure.from_parsed_mmcif(
        cif,
        fix_mse_residues=False,
        fix_arginines=fix_arginines,
        include_water=keep_ligands_water,
        include_other=False,
        model_id=model_id,
    )
    if not keep_ligands_water:
      struc = struc.filter_to_entity_type(protein=True)

    model_id_str = ''
    if model_id == structure.ModelID.FIRST:
      model_id_str = cif.get_array(
          '_atom_site.pdbx_PDB_model_num', dtype=object, gather=slice(1)
      )[0]
    index, wanted_chain_ids, layout = mmcif_utils.filter(
        cif,
        include_nucleotides=keep_ligands_water,
        include_ligands=keep_ligands_water,
        include_water=keep_ligands_water,
        include_other=False,
        model_id=model_id_str,
    )
    with self.subTest('num_models'):
      self.assertEqual(struc.num_models, layout.num_models())
    with self.subTest('num_chains'):
      self.assertEqual(struc.num_chains, layout.num_chains())
    with self.subTest('num_atoms'):
      self.assertEqual(struc.num_atoms, layout.num_atoms())

    label_comp_id = _atom_site_str_array(cif, 'label_comp_id', gather=index[0])
    label_atom_id = _atom_site_str_array(cif, 'label_atom_id', gather=index[0])

    atom_id = cif.get_array('_atom_site.id', dtype=object, gather=index[0])
    atom_x = cif.get_array('_atom_site.Cartn_x', dtype=np.float32, gather=index)
    atom_y = cif.get_array('_atom_site.Cartn_y', dtype=np.float32, gather=index)
    atom_z = cif.get_array('_atom_site.Cartn_z', dtype=np.float32, gather=index)
    mmcif_utils.fix_residues(
        layout,
        comp_id=label_comp_id,
        atom_id=label_atom_id,
        atom_x=atom_x[0],
        atom_y=atom_y[0],
        atom_z=atom_z[0],
        fix_arg=fix_arginines,
    )

    if model_id != structure.ModelID.ALL:
      atom_x = atom_x[0]
      atom_y = atom_y[0]
      atom_z = atom_z[0]

    with self.subTest('chain_ids'):
      self.assertCountEqual(struc.chains, wanted_chain_ids)

    dif = np.where(struc.res_name != label_comp_id)[0]
    with self.subTest('res_name'):
      np.testing.assert_array_equal(
          struc.res_name[dif], label_comp_id[dif], err_msg=str(atom_id[dif])
      )
    with self.subTest('atom_name'):
      np.testing.assert_array_equal(
          struc.atom_name[dif], label_atom_id[dif], err_msg=str(atom_id[dif])
      )

    if cif_name not in _DIFFERENT_ALT_LOC_POS:
      with self.subTest('coords'):
        np.testing.assert_array_equal(struc.atom_x, atom_x)
        np.testing.assert_array_equal(struc.atom_y, atom_y)
        np.testing.assert_array_equal(struc.atom_z, atom_z)

  @parameterized.parameters(
      dict(cif_name='100d', num_atoms=489),
      dict(cif_name='170d', num_atoms=760),
      dict(cif_name='1al4', num_atoms=672),
      dict(cif_name='1bna', num_atoms=566),
      dict(cif_name='1h9h', num_atoms=2005),
      dict(cif_name='1ign', num_atoms=4908),
      dict(cif_name='1ob2', num_atoms=4810),
      dict(cif_name='1wib', num_atoms=1326),
      dict(cif_name='2pv7', num_atoms=4934),
      dict(cif_name='3b9f', num_atoms=5736),
      dict(cif_name='3co3', num_atoms=511),
      dict(cif_name='3ibw', num_atoms=1357),
      dict(cif_name='3j7l', num_atoms=3613),
      dict(cif_name='3vyp', num_atoms=4972),
      dict(cif_name='4gxy', num_atoms=3685),
      dict(cif_name='5ffz', num_atoms=1051),
      dict(cif_name='5inz', num_atoms=968),
      dict(cif_name='5jsu', num_atoms=12435),
      dict(cif_name='6bea', num_atoms=3406),
      dict(cif_name='6d4l', num_atoms=345),
  )
  def test_alt_loc_sizes_are_correct(self, cif_name, num_atoms):
    cif = _read_test_cif_file(cif_name)
    model_id = _atom_site_str_array(cif, 'pdbx_PDB_model_num', gather=slice(1))[
        0
    ]
    index, _, layout = mmcif_utils.filter(
        cif,
        include_nucleotides=True,
        include_ligands=True,
        include_water=True,
        include_other=False,
        model_id=model_id,
    )
    self.assertEqual(index.shape, (layout.num_models(), layout.num_atoms()))
    self.assertEqual(layout.num_atoms(), num_atoms)

  def test_filter_other_chains(self):
    cif = _read_test_cif_file('170d')
    index, _, layout = mmcif_utils.filter(
        cif,
        include_nucleotides=True,
        include_ligands=True,
        include_water=True,
        include_other=True,
        model_id='1',
    )
    self.assertEqual(index.shape, (layout.num_models(), layout.num_atoms()))
    self.assertEqual(layout.num_atoms(), 760)

  @parameterized.parameters(
      dict(cif_name='5inz', chain_residue={('B', '8'), ('C', '18')}),
      dict(cif_name='5jsu', chain_residue={('B', '211')}),
  )
  def test_fix_arg(self, cif_name, chain_residue):
    cif = _read_test_cif_file(cif_name)
    index, _, layout = mmcif_utils.filter(
        cif,
        include_nucleotides=True,
        include_ligands=True,
        include_water=True,
        model_id='1',
    )

    chain_id = _atom_site_str_array(cif, 'label_asym_id', gather=index[0])
    seq_id = _atom_site_str_array(cif, 'label_seq_id', gather=index[0])
    label_comp_id = _atom_site_str_array(cif, 'label_comp_id', gather=index[0])
    label_atom_id_sans = _atom_site_str_array(
        cif, 'label_atom_id', gather=index[0]
    )
    atom_x = cif.get_array('_atom_site.Cartn_x', dtype=np.float32, gather=index)
    atom_y = cif.get_array('_atom_site.Cartn_y', dtype=np.float32, gather=index)
    atom_z = cif.get_array('_atom_site.Cartn_z', dtype=np.float32, gather=index)

    label_atom_id_with = label_atom_id_sans.copy()
    mmcif_utils.fix_residues(
        layout,
        comp_id=label_comp_id,
        atom_id=label_atom_id_with,
        atom_x=atom_x[0],
        atom_y=atom_y[0],
        atom_z=atom_z[0],
        fix_arg=True,
    )
    for residue_index in range(layout.num_residues()):
      res_start, res_end = layout.atom_range(residue_index)
      name = chain_id[res_start], seq_id[res_start]
      if name in chain_residue:
        name_map = dict(
            zip(
                label_atom_id_sans[res_start:res_end],
                label_atom_id_with[res_start:res_end],
                strict=True,
            )
        )
        expected_map = {'NH1': 'NH2', 'HH11': 'HH21', 'HH12': 'HH22'}
        for k, v in name_map.items():
          if k in expected_map:
            self.assertEqual(v, expected_map[k], msg=f'{name} - {k} - {v}')
          elif v in expected_map:
            self.assertEqual(k, expected_map[v], msg=f'{name} - {k} - {v}')
          else:
            self.assertEqual(v, k, msg=f'{name} - {k} - {v}')
      else:
        np.testing.assert_array_equal(
            label_atom_id_sans[res_start:res_end],
            label_atom_id_with[res_start:res_end],
            err_msg=f'{name}',
        )

  def test_6qh9_does_not_crash(self):
    # 6qh9 - last residue has an alt-loc.
    cif = _read_test_cif_file('6qh9')
    index, _, layout = mmcif_utils.filter(
        cif,
        include_nucleotides=True,
        include_ligands=True,
        include_water=False,
        model_id='1',
    )
    label_comp_id = _atom_site_str_array(cif, 'label_comp_id', gather=index[0])
    label_atom_id = _atom_site_str_array(cif, 'label_atom_id', gather=index[0])
    atom_x = cif.get_array('_atom_site.Cartn_x', dtype=np.float32, gather=index)
    atom_y = cif.get_array('_atom_site.Cartn_y', dtype=np.float32, gather=index)
    atom_z = cif.get_array('_atom_site.Cartn_z', dtype=np.float32, gather=index)

    mmcif_utils.fix_residues(
        layout,
        comp_id=label_comp_id,
        atom_id=label_atom_id,
        atom_x=atom_x[0],
        atom_y=atom_y[0],
        atom_z=atom_z[0],
        fix_arg=True,
    )

  @parameterized.parameters(range(20))
  def test_get_model_by_model_id(self, model_index):
    cif = _read_test_cif_file('1wib')
    index_all, _, _ = mmcif_utils.filter(
        cif,
        include_nucleotides=True,
        include_ligands=True,
        include_water=True,
        model_id='',
    )

    index_one, _, layout_one = mmcif_utils.filter(
        cif,
        include_nucleotides=True,
        include_ligands=True,
        include_water=False,
        model_id=str(model_index + 1),
    )
    self.assertEqual(index_one.shape, (1, layout_one.num_atoms()))
    np.testing.assert_array_equal(index_all[model_index], index_one[0])

  def test_unknown_model_id(self):
    cif = _read_test_cif_file('1wib')
    with self.assertRaises(ValueError):
      mmcif_utils.filter(
          cif,
          include_nucleotides=True,
          include_ligands=True,
          include_water=False,
          model_id=str(28),
      )

  @parameterized.parameters(range(20))
  def test_get_model_by_model_id_layout(self, model_index):
    cif = _read_test_cif_file('1wib')
    layout_all = mmcif_utils.read_layout(cif, model_id='')
    self.assertEqual(layout_all.num_models(), 20)

    layout_one = mmcif_utils.read_layout(cif, model_id=str(model_index + 1))
    self.assertEqual(layout_one.num_models(), 1)

    offset = model_index * layout_all.num_atoms()
    self.assertEqual(layout_all.num_chains(), layout_one.num_chains())
    layout_all_residues = np.array(layout_all.residues()) + offset
    layout_one_residues = np.array(layout_one.residues())
    np.testing.assert_array_equal(layout_all_residues, layout_one_residues)

  def test_unknown_model_id_layout(self):
    cif = _read_test_cif_file('1wib')
    with self.assertRaises(ValueError):
      mmcif_utils.filter(
          cif,
          include_nucleotides=True,
          include_ligands=True,
          include_water=False,
          model_id=str(28),
      )

  def test_selected_polymer_residue_mask(self):
    cif = _read_test_cif_file('5by8')
    atom_site_all_models, _, layout = mmcif_utils.filter(
        cif,
        include_nucleotides=True,
        include_ligands=True,
        include_water=False,
        include_other=True,
        model_id='1',
    )

    def _first_model_string_array(col: str) -> np.ndarray:
      return cif.get_array(col, dtype=object, gather=atom_site_all_models[0])

    label_asym_ids = _first_model_string_array('_atom_site.label_asym_id')
    label_seq_ids = _first_model_string_array('_atom_site.label_seq_id')
    label_comp_ids = _first_model_string_array('_atom_site.label_comp_id')

    def _get_poly_seq_scheme_col(col: str) -> np.ndarray:
      return cif.get_array(key=f'_pdbx_poly_seq_scheme.{col}', dtype=object)

    poly_seq_asym_ids = _get_poly_seq_scheme_col('asym_id')
    poly_seq_seq_ids = _get_poly_seq_scheme_col('seq_id')
    poly_seq_mon_ids = _get_poly_seq_scheme_col('mon_id')

    mask = mmcif_utils.selected_polymer_residue_mask(
        layout=layout,
        atom_site_label_asym_ids=label_asym_ids[layout.residue_starts()],
        atom_site_label_seq_ids=label_seq_ids[layout.residue_starts()],
        atom_site_label_comp_ids=label_comp_ids[layout.residue_starts()],
        poly_seq_asym_ids=poly_seq_asym_ids,
        poly_seq_seq_ids=poly_seq_seq_ids,
        poly_seq_mon_ids=poly_seq_mon_ids,
    )

    exp_hetero_res_mask = np.zeros_like(mask, dtype=bool)
    exp_hetero_res_mask[[65, 92, 111, 117, 128, 146, 151, 175, 233, 238]] = True
    num_hetero_res = sum(exp_hetero_res_mask)

    np.testing.assert_array_equal(
        np.array(mask)[exp_hetero_res_mask], [False] * num_hetero_res
    )
    np.testing.assert_array_equal(
        np.array(mask)[~exp_hetero_res_mask], [True] * layout.num_residues()
    )
    self.assertLen(mask, layout.num_residues() + num_hetero_res)

  def test_selected_ligand_residue_mask(self):
    cif = _read_test_cif_file('4uz3')
    atom_site_all_models, _, layout = mmcif_utils.filter(
        cif,
        include_nucleotides=True,
        include_ligands=True,
        include_water=False,
        include_other=True,
        model_id='1',
    )

    def _first_model_string_array(col: str) -> np.ndarray:
      return cif.get_array(col, dtype=object, gather=atom_site_all_models[0])

    label_asym_ids = _first_model_string_array('_atom_site.label_asym_id')
    label_seq_ids = _first_model_string_array('_atom_site.label_seq_id')
    label_comp_ids = _first_model_string_array('_atom_site.label_comp_id')
    auth_seq_ids = _first_model_string_array('_atom_site.auth_seq_id')
    ins_codes = _first_model_string_array('_atom_site.pdbx_PDB_ins_code')

    def _get_nonpoly_scheme_col(col: str) -> np.ndarray:
      return cif.get_array(key=f'_pdbx_nonpoly_scheme.{col}', dtype=object)

    nonpoly_asym_ids = _get_nonpoly_scheme_col('asym_id')
    nonpoly_auth_seq_ids = _get_nonpoly_scheme_col('pdb_seq_num')
    nonpoly_pdb_ins_codes = _get_nonpoly_scheme_col('pdb_ins_code')
    nonpoly_mon_ids = _get_nonpoly_scheme_col('mon_id')

    def _get_branch_scheme_col(col: str) -> np.ndarray:
      return cif.get_array(key=f'_pdbx_branch_scheme.{col}', dtype=object)

    branch_asym_ids = _get_branch_scheme_col('asym_id')
    branch_auth_seq_ids = _get_branch_scheme_col('pdb_seq_num')
    branch_pdb_ins_codes = np.array(['.'] * len(branch_asym_ids), dtype=object)
    branch_mon_ids = _get_branch_scheme_col('mon_id')

    nonpoly_mask, branch_mask = mmcif_utils.selected_ligand_residue_mask(
        layout=layout,
        atom_site_label_asym_ids=label_asym_ids[layout.residue_starts()],
        atom_site_label_seq_ids=label_seq_ids[layout.residue_starts()],
        atom_site_label_comp_ids=label_comp_ids[layout.residue_starts()],
        atom_site_auth_seq_ids=auth_seq_ids[layout.residue_starts()],
        atom_site_pdbx_pdb_ins_codes=ins_codes[layout.residue_starts()],
        nonpoly_asym_ids=nonpoly_asym_ids,
        nonpoly_auth_seq_ids=nonpoly_auth_seq_ids,
        nonpoly_pdb_ins_codes=nonpoly_pdb_ins_codes,
        nonpoly_mon_ids=nonpoly_mon_ids,
        branch_asym_ids=branch_asym_ids,
        branch_auth_seq_ids=branch_auth_seq_ids,
        branch_pdb_ins_codes=branch_pdb_ins_codes,
        branch_mon_ids=branch_mon_ids,
    )

    # No heterogeneous nonpoly residues.
    exp_nonpoly_mask = np.ones_like(nonpoly_mask, dtype=bool)
    np.testing.assert_array_equal(nonpoly_mask, exp_nonpoly_mask)

    # There are two heterogeneous ligand residues on positions 0 and 13.
    exp_branch_mask = np.ones_like(branch_mask, dtype=bool)
    exp_branch_mask[[0, 13]] = False

    np.testing.assert_array_equal(branch_mask, exp_branch_mask)


if __name__ == '__main__':
  absltest.main()
