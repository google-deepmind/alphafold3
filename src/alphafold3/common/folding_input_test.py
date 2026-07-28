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

"""Unit tests for folding_input validation and chain API consistency."""

import json
import pathlib
import tempfile

from absl.testing import absltest
from absl.testing import parameterized

from alphafold3.common import folding_input


_MINIMAL_MMCIF = """\
data_test
_entry.id test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 1 ? 0.000 0.000 0.000 1.00 0.00 1 A 1
#
"""


def _protein_json(
    *,
    templates: list[dict] | None = None,
    extra_sequences: list[dict] | None = None,
) -> str:
  protein = {
      'id': 'A',
      'sequence': 'ACDE',
      'modifications': [],
      'unpairedMsa': '',
      'pairedMsa': '',
  }
  if templates is not None:
    protein['templates'] = templates
  sequences = [{'protein': protein}]
  if extra_sequences:
    sequences.extend(extra_sequences)
  return json.dumps({
      'name': 'test',
      'sequences': sequences,
      'modelSeeds': [1],
      'dialect': 'alphafold3',
      'version': 4,
  })


class TemplateValidationTest(parameterized.TestCase):

  def test_rejects_mismatched_query_and_template_indices(self):
    with self.assertRaisesRegex(ValueError, 'same length'):
      folding_input.Input.from_json(
          _protein_json(
              templates=[{
                  'mmcif': _MINIMAL_MMCIF,
                  'queryIndices': [0, 1],
                  'templateIndices': [0],
              }]
          )
      )

  def test_rejects_template_without_mmcif_or_path(self):
    with self.assertRaisesRegex(ValueError, 'mmcif|mmcifPath'):
      folding_input.Input.from_json(
          _protein_json(
              templates=[{
                  'queryIndices': [0],
                  'templateIndices': [0],
              }]
          )
      )

  def test_rejects_empty_template_mmcif(self):
    with self.assertRaisesRegex(ValueError, 'non-empty'):
      folding_input.Template(mmcif='', query_to_template_map={0: 0})

  def test_accepts_matching_template_indices(self):
    fold_input = folding_input.Input.from_json(
        _protein_json(
            templates=[{
                'mmcif': _MINIMAL_MMCIF,
                'queryIndices': [0, 1],
                'templateIndices': [10, 11],
            }]
        )
    )
    self.assertEqual(
        fold_input.protein_chains[0].templates[0].query_to_template_map,
        {0: 10, 1: 11},
    )


class RnaChainFillMissingFieldsTest(absltest.TestCase):

  def test_preserves_raw_sequence_with_unknown_modification(self):
    # Unknown CCD codes map to 'N' via the `.sequence` property. Filling missing
    # fields must keep the original residue letters in `_sequence`.
    chain = folding_input.RnaChain(
        id='R',
        sequence='ACGU',
        modifications=[('1MA', 2)],
        unpaired_msa=None,
    )
    filled = chain.fill_missing_fields()
    self.assertEqual(filled._sequence, 'ACGU')
    self.assertEqual(filled.unpaired_msa, '')
    self.assertEqual(filled.modifications, (('1MA', 2),))


class DnaChainModificationsPropertyTest(absltest.TestCase):

  def test_modifications_is_a_property(self):
    chain = folding_input.DnaChain(
        id='D', sequence='ACGT', modifications=[('6OG', 1)]
    )
    mods = chain.modifications
    self.assertEqual(mods, (('6OG', 1),))
    self.assertFalse(callable(mods))
    # Iterable like RnaChain.modifications.
    self.assertEqual(list(chain.modifications), [('6OG', 1)])


class LigandValidationTest(parameterized.TestCase):

  def test_rejects_empty_ccd_codes_list(self):
    with self.assertRaisesRegex(ValueError, 'non-empty'):
      folding_input.Ligand(id='L', ccd_ids=[])

  def test_rejects_empty_ccd_codes_in_json(self):
    payload = _protein_json(
        extra_sequences=[{
            'ligand': {'id': 'L', 'ccdCodes': []},
        }]
    )
    with self.assertRaisesRegex(ValueError, 'non-empty'):
      folding_input.Input.from_json(payload)

  def test_rejects_empty_smiles(self):
    with self.assertRaisesRegex(ValueError, 'non-empty'):
      folding_input.Ligand(id='L', smiles='')

  def test_accepts_valid_ccd_ligand(self):
    ligand = folding_input.Ligand(id='L', ccd_ids=['ATP'])
    self.assertEqual(ligand.ccd_ids, ('ATP',))


class LoadFoldInputsFromDirTest(absltest.TestCase):

  def test_empty_directory_raises(self):
    with tempfile.TemporaryDirectory() as tmp_dir:
      with self.assertRaisesRegex(ValueError, 'No JSON files'):
        list(folding_input.load_fold_inputs_from_dir(tmp_dir))

  def test_loads_json_files(self):
    with tempfile.TemporaryDirectory() as tmp_dir:
      path = pathlib.Path(tmp_dir) / 'fold.json'
      path.write_text(_protein_json())
      inputs = list(folding_input.load_fold_inputs_from_dir(tmp_dir))
      self.assertLen(inputs, 1)
      self.assertEqual(inputs[0].name, 'test')


if __name__ == '__main__':
  absltest.main()
