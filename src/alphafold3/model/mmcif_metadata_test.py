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
from alphafold3.constants import chemical_components
from alphafold3.constants import mmcif_names
from alphafold3.model import mmcif_metadata
import numpy as np

_EXP_KEYS = (
    'data_',
    '_atom_site.B_iso_or_equiv',
    '_atom_site.Cartn_x',
    '_atom_site.Cartn_y',
    '_atom_site.Cartn_z',
    '_atom_site.auth_asym_id',
    '_atom_site.auth_seq_id',
    '_atom_site.group_PDB',
    '_atom_site.id',
    '_atom_site.label_alt_id',
    '_atom_site.label_asym_id',
    '_atom_site.label_atom_id',
    '_atom_site.label_comp_id',
    '_atom_site.label_entity_id',
    '_atom_site.label_seq_id',
    '_atom_site.occupancy',
    '_atom_site.pdbx_PDB_ins_code',
    '_atom_site.pdbx_PDB_model_num',
    '_atom_site.type_symbol',
    '_atom_type.symbol',
    '_audit_author.name',
    '_audit_author.pdbx_ordinal',
    '_audit_conform.dict_location',
    '_audit_conform.dict_name',
    '_audit_conform.dict_version',
    '_chem_comp.formula',
    '_chem_comp.formula_weight',
    '_chem_comp.id',
    '_chem_comp.mon_nstd_flag',
    '_chem_comp.name',
    '_chem_comp.pdbx_synonyms',
    '_chem_comp.type',
    '_citation.book_publisher',
    '_citation.country',
    '_citation.id',
    '_citation.journal_full',
    '_citation.journal_id_ASTM',
    '_citation.journal_id_CSD',
    '_citation.journal_id_ISSN',
    '_citation.journal_volume',
    '_citation.page_first',
    '_citation.page_last',
    '_citation.pdbx_database_id_DOI',
    '_citation.pdbx_database_id_PubMed',
    '_citation.title',
    '_citation.year',
    '_citation_author.citation_id',
    '_citation_author.name',
    '_citation_author.ordinal',
    '_entity.id',
    '_entity.pdbx_description',
    '_entity.type',
    '_entity_poly.entity_id',
    '_entity_poly.pdbx_strand_id',
    '_entity_poly.type',
    '_entity_poly_seq.entity_id',
    '_entity_poly_seq.hetero',
    '_entity_poly_seq.mon_id',
    '_entity_poly_seq.num',
    '_entry.id',
    '_ma_data.content_type',
    '_ma_data.id',
    '_ma_data.name',
    '_ma_model_list.data_id',
    '_ma_model_list.model_group_id',
    '_ma_model_list.model_group_name',
    '_ma_model_list.model_id',
    '_ma_model_list.model_name',
    '_ma_model_list.model_type',
    '_ma_model_list.ordinal_id',
    '_ma_protocol_step.method_type',
    '_ma_protocol_step.ordinal_id',
    '_ma_protocol_step.protocol_id',
    '_ma_protocol_step.step_id',
    '_ma_qa_metric.id',
    '_ma_qa_metric.mode',
    '_ma_qa_metric.name',
    '_ma_qa_metric.software_group_id',
    '_ma_qa_metric.type',
    '_ma_qa_metric_global.metric_id',
    '_ma_qa_metric_global.metric_value',
    '_ma_qa_metric_global.model_id',
    '_ma_qa_metric_global.ordinal_id',
    '_ma_qa_metric_local.label_asym_id',
    '_ma_qa_metric_local.label_comp_id',
    '_ma_qa_metric_local.label_seq_id',
    '_ma_qa_metric_local.metric_id',
    '_ma_qa_metric_local.metric_value',
    '_ma_qa_metric_local.model_id',
    '_ma_qa_metric_local.ordinal_id',
    '_ma_software_group.group_id',
    '_ma_software_group.ordinal_id',
    '_ma_software_group.software_id',
    '_ma_target_entity.data_id',
    '_ma_target_entity.entity_id',
    '_ma_target_entity.origin',
    '_ma_target_entity_instance.asym_id',
    '_ma_target_entity_instance.details',
    '_ma_target_entity_instance.entity_id',
    '_pdbx_data_usage.details',
    '_pdbx_data_usage.id',
    '_pdbx_data_usage.type',
    '_pdbx_data_usage.url',
    '_pdbx_poly_seq_scheme.asym_id',
    '_pdbx_poly_seq_scheme.auth_seq_num',
    '_pdbx_poly_seq_scheme.entity_id',
    '_pdbx_poly_seq_scheme.hetero',
    '_pdbx_poly_seq_scheme.mon_id',
    '_pdbx_poly_seq_scheme.pdb_ins_code',
    '_pdbx_poly_seq_scheme.pdb_seq_num',
    '_pdbx_poly_seq_scheme.pdb_strand_id',
    '_pdbx_poly_seq_scheme.seq_id',
    '_software.classification',
    '_software.date',
    '_software.description',
    '_software.name',
    '_software.pdbx_ordinal',
    '_software.type',
    '_software.version',
    '_struct_asym.entity_id',
    '_struct_asym.id',
)


class MmcifMetadataTest(parameterized.TestCase):

  def setUp(self):
    super().setUp()
    self.struc = structure.from_sequences_and_bonds(
        sequences=['DEEPMIND', 'AGCT'],
        chain_types=[mmcif_names.PROTEIN_CHAIN, mmcif_names.DNA_CHAIN],
        sequence_formats=[
            structure.SequenceFormat.FASTA,
            structure.SequenceFormat.FASTA,
        ],
        bonded_atom_pairs=None,
        ccd=chemical_components.Ccd(),
        release_date=None,
    )
    self.golden_cif = data.Data(resources.ROOT / 'model/test_data/').load(
        'af3_exp_output.cif'
    )

  @parameterized.named_parameters(
      ('keep_license', True),
      ('discard_license', False),
  )
  def test_metadata_added(self, keep_license: bool):
    cif_with_metadata = mmcif_metadata.add_metadata_to_mmcif(
        self.struc.to_mmcif_dict(),
        version='TEST_VERSION',
        model_id=b'TEST_MODEL_ID',
        keep_license=keep_license,
    )

    self.assertCountEqual(cif_with_metadata.keys(), _EXP_KEYS)

    self.assertIn('(TEST_MODEL_ID)', cif_with_metadata['_software.version'][0])

    # Check for weird characters like U+2009 (Thin Space).
    cif_string = cif_with_metadata.to_string()
    self.assertTrue(
        cif_string.replace('\n', '').isprintable(),
        msg='The CIF with metadata contains unprintable Unicode characters',
    )
    with self.subTest('license_and_disclaimer'):
      if keep_license:
        self.assertIn('Non-commercial use', cif_string)  # License.
        self.assertIn('not approved for clinical', cif_string)  # Disclaimer.
        self.assertLen(cif_with_metadata['_pdbx_data_usage.details'], 2)
        self.assertLen(cif_with_metadata['_pdbx_data_usage.details'][0], 476)
        self.assertLen(cif_with_metadata['_pdbx_data_usage.details'][1], 286)
      else:
        self.assertNotIn('Non-commercial use', cif_string)  # License.
        self.assertIn('not approved for clinical', cif_string)  # Disclaimer.
        self.assertLen(cif_with_metadata['_pdbx_data_usage.details'], 1)
        self.assertLen(cif_with_metadata['_pdbx_data_usage.details'][0], 286)

  def test_add_legal_comment(self):
    original_cif = self.struc.to_mmcif_dict().to_string()
    cif_with_comment = mmcif_metadata.add_legal_comment(original_cif)
    with self.subTest('comment_added'):
      self.assertStartsWith(
          cif_with_comment,
          '# By using this file you agree to the legally binding terms',
      )
    with self.subTest('cif_unchanged'):
      self.assertEndsWith(cif_with_comment, original_cif)
    with self.subTest('valid_cif'):
      structure.from_mmcif(cif_with_comment)

  def test_golden(self):
    # Structure will remove metadata and comments, check it is added back.
    struc = structure.from_mmcif(self.golden_cif)
    cif_with_metadata = mmcif_metadata.add_metadata_to_mmcif(
        struc.to_mmcif_dict(), version='TEST_VERSION', model_id=b'TEST_MODEL_ID'
    )
    self.assertEqual(cif_with_metadata.to_string(), self.golden_cif)

  def test_per_residue_plddt(self):
    struc = structure.from_mmcif(self.golden_cif)
    cif_with_metadata = mmcif_metadata.add_metadata_to_mmcif(
        struc.to_mmcif_dict(), version='TEST_VERSION', model_id=b'TEST_MODEL_ID'
    )
    per_residue_plddt = cif_with_metadata['_ma_qa_metric_local.metric_value']
    self.assertLen(per_residue_plddt, struc.num_residues(count_unresolved=True))
    plddt_table_residues = zip(
        cif_with_metadata['_ma_qa_metric_local.label_asym_id'],
        cif_with_metadata['_ma_qa_metric_local.label_comp_id'],
        cif_with_metadata['_ma_qa_metric_local.label_seq_id'],
    )
    polymer_residues = zip(
        cif_with_metadata['_pdbx_poly_seq_scheme.asym_id'],
        cif_with_metadata['_pdbx_poly_seq_scheme.mon_id'],
        cif_with_metadata['_pdbx_poly_seq_scheme.seq_id'],
    )
    non_polymer_residues = zip(
        cif_with_metadata['_pdbx_nonpoly_scheme.asym_id'],
        cif_with_metadata['_pdbx_nonpoly_scheme.mon_id'],
        # Non-polymer residues have no internal sequence ID.
        ['.'] * len(cif_with_metadata['_pdbx_nonpoly_scheme.asym_id']),
    )
    self.assertCountEqual(
        list(plddt_table_residues),
        list(polymer_residues) + list(non_polymer_residues),
    )

    met_plddt = [74.15, 82.52, 85.42, 79.50, 73.22, 65.18, 61.99, 56.89]
    glu_plddt = [80.50, 87.21, 91.56, 86.96, 76.14, 65.16, 63.54, 56.95, 58.82]
    self.assertEqual(per_residue_plddt[0], f'{np.mean(met_plddt):.2f}')
    self.assertEqual(per_residue_plddt[1], f'{np.mean(glu_plddt):.2f}')


if __name__ == '__main__':
  absltest.main()
