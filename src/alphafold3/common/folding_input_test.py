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
import copy
import gzip
import io
import json
import lzma
import pathlib
from typing import Any
from unittest import mock

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.common import folding_input
from alphafold3.common import resources
from alphafold3.common.testing import data
from alphafold3.constants import chemical_components
from alphafold3.constants import mmcif_names
from alphafold3.cpp import cif_dict
from etils.epath import testing as epath_testing
import numpy as np
import zstandard as zstd

# The seed to use when mocking _sample_rng_seed().
_SAMPLE_RNG_SEED = 1234567890


class InputTest(parameterized.TestCase):

  def setUp(self):
    super().setUp()

    mock.patch.object(
        folding_input,
        '_sample_rng_seed',
        return_value=_SAMPLE_RNG_SEED,
    ).start()

  def tearDown(self):
    super().tearDown()
    mock.patch.stopall()

  def test_template_equality(self):
    template_1 = folding_input.Template(
        mmcif='irrelevant', query_to_template_map={0: 1, 1: 4, 2: 3}
    )
    template_2 = folding_input.Template(
        mmcif='irrelevant', query_to_template_map={1: 4, 2: 3, 0: 1}
    )

    with self.subTest('eq'):
      self.assertEqual(template_1, template_1)
      self.assertEqual(template_1, template_2)
    with self.subTest('hash'):
      self.assertEqual(hash(template_1), hash(template_1))
      self.assertEqual(hash(template_1), hash(template_2))

  def test_template_mapping_order(self):
    template_1 = folding_input.Template(
        mmcif='irrelevant',
        query_to_template_map={4: 4, 3: 3, 2: 2, 1: 1},
    )
    self.assertEqual(template_1.query_to_template_map, {4: 4, 3: 3, 2: 2, 1: 1})

  @parameterized.named_parameters(
      ('empty_query_to_template_map', {}),
      ('one_entry_query_to_template_map', {1: 1, 2: 2, 3: 3, 4: 4}),
  )
  def test_template_dict_roundtrip(self, query_to_template_map):
    template = folding_input.Template(
        mmcif='irrelevant', query_to_template_map=query_to_template_map
    )
    self.assertEqual(
        template, folding_input.Template.from_dict(template.to_dict())
    )

  def test_templates_from_dict_with_template_indices_none(self):
    template_dict = {
        'mmcif': 'irrelevant',
        'queryIndices': [],
        'templateIndices': None,
    }
    template = folding_input.Template.from_dict(template_dict)
    self.assertEqual(
        template,
        folding_input.Template(mmcif='irrelevant', query_to_template_map={}),
    )

  def test_protein_to_ccd_sequence(self):
    protein_chain = folding_input.ProteinChain(
        id='A', sequence='ABCDEFGHIJ', ptms=[('HY3', 1), ('P1L', 5)]
    )
    self.assertEqual(
        protein_chain.to_ccd_sequence(),
        ['HY3', 'UNK', 'CYS', 'ASP', 'P1L', 'PHE', 'GLY', 'HIS', 'ILE', 'UNK'],
    )

  def test_rna_to_ccd_sequence(self):
    rna_chain = folding_input.RnaChain(
        id='A', sequence='AGCUN', modifications=[('2MG', 2), ('5MC', 4)]
    )
    self.assertEqual(
        rna_chain.to_ccd_sequence(),
        ['A', '2MG', 'C', '5MC', 'N'],
    )

  def test_dna_to_ccd_sequence(self):
    dna_chain = folding_input.DnaChain(
        id='A', sequence='AGCTN', modifications=[('6OG', 1), ('6MA', 2)]
    )
    self.assertEqual(
        dna_chain.to_ccd_sequence(),
        ['6OG', '6MA', 'DC', 'DT', 'DN'],
    )

  def test_protein_to_single_letter_sequence(self):
    protein_chain = folding_input.ProteinChain(
        id='A', sequence='ABCDEFGHIJ', ptms=[('UNK', 1), ('P1L', 5)]
    )
    self.assertEqual(protein_chain.sequence, 'XXCDCFGHIX')

  def test_protein_to_single_letter_sequence_unknowns(self):
    protein_chain = folding_input.ProteinChain(
        id='A', sequence='BJOUXZ', ptms=[]
    )
    self.assertEqual(protein_chain.sequence, 'XXXXXX')

  def test_rna_to_single_letter_sequence(self):
    rna_chain = folding_input.RnaChain(
        id='A', sequence='AGCUNZ', modifications=[('N', 2), ('5MC', 4)]
    )
    self.assertEqual(rna_chain.sequence, 'ANCCNN')

  def test_dna_to_single_letter_sequence(self):
    dna_chain = folding_input.DnaChain(
        id='A', sequence='AGCTNZ', modifications=[('DN', 1), ('6MA', 2)]
    )
    self.assertEqual(dna_chain.sequence, 'NACTNN')

  @parameterized.named_parameters(
      ('modified_sequence', 'CC', [], None, None, None),
      ('modified_ptms', 'AB', [('HY3', 1)], None, None, None),
      ('modified_unpaired_msa', 'AB', [], '>query\nAB\n>u1\n-B', None, None),
      ('modified_paired_msa', 'AB', [], None, '>query\nAB\n>p1\n-B', None),
      ('modified_templates', 'AB', [], None, None, []),
  )
  def test_protein_hash_without_id(
      self, sequence, ptms, unpaired_msa, paired_msa, templates
  ):
    protein_1 = folding_input.ProteinChain(
        id='A',
        sequence='AB',
        ptms=[],
        unpaired_msa=None,
        paired_msa=None,
        templates=None,
    )
    protein_2 = folding_input.ProteinChain(
        id='A',
        sequence=sequence,
        ptms=ptms,
        unpaired_msa=unpaired_msa,
        paired_msa=paired_msa,
        templates=templates,
    )
    self.assertNotEqual(
        protein_1.hash_without_id(), protein_2.hash_without_id()
    )

  @parameterized.named_parameters(
      ('modified_sequence', 'GG', [], None),
      ('modified_modifications', 'AC', [('HY3', 1)], None),
      ('modified_unpaired_msa', 'AC', [], '>query\nAB\n>u1\n-C'),
  )
  def test_rna_hash_without_id(self, sequence, modifications, unpaired_msa):
    rna_1 = folding_input.RnaChain(
        id='A',
        sequence='AC',
        modifications=[],
        unpaired_msa=None,
    )
    rna_2 = folding_input.RnaChain(
        id='A',
        sequence=sequence,
        modifications=modifications,
        unpaired_msa=unpaired_msa,
    )
    self.assertNotEqual(rna_1.hash_without_id(), rna_2.hash_without_id())

  @parameterized.named_parameters(
      ('modified_sequence', 'GG', []),
      ('modified_modifications', 'AC', [('HY3', 1)]),
  )
  def test_dna_hash_without_id(self, sequence, modifications):
    dna_1 = folding_input.DnaChain(
        id='A',
        sequence='AC',
        modifications=[],
    )
    dna_2 = folding_input.DnaChain(
        id='A',
        sequence=sequence,
        modifications=modifications,
    )
    self.assertNotEqual(dna_1.hash_without_id(), dna_2.hash_without_id())

  def test_ligand_hash_without_id(self):
    self.assertNotEqual(
        folding_input.Ligand(id='A', ccd_ids=['ARG']),
        folding_input.Ligand(id='A', ccd_ids=['GLY']),
    )
    self.assertNotEqual(
        folding_input.Ligand(id='A', smiles='C1=CC1'),
        folding_input.Ligand(id='A', smiles='C=C'),
    )

  @parameterized.named_parameters(
      ('templates_None', None),
      ('templates_empty', []),
      (
          'templates_non_empty',
          [folding_input.Template(mmcif='tst', query_to_template_map={0: 1})],
      ),
      (
          'templates_with_empty_map',
          [folding_input.Template(mmcif='tst', query_to_template_map={})],
      ),
  )
  def test_protein_dict_roundtrip(self, templates):
    protein_chain = folding_input.ProteinChain(
        id='A',
        sequence='ABCDEFGHIJ',
        ptms=[('HY3', 1), ('P1L', 5)],
        description='Protein chain',
        paired_msa='>query\nABCDEFGHIJ\n>p1\nABCDEFGHIJ',
        unpaired_msa='>query\nABCDEFGHIJ\n>u1\nABCDEFGHIJ',
        templates=templates,
    )
    self.assertEqual(
        protein_chain,
        protein_chain.from_dict(protein_chain.to_dict()),
    )

  def test_rna_dict_roundtrip(self):
    rna_chain = folding_input.RnaChain(
        id='A',
        sequence='AGCUN',
        modifications=[('2MG', 2), ('5MC', 4)],
        description='RNA chain',
        unpaired_msa='>query\nAGCUN\n>u1\nA-C-N',
    )
    self.assertEqual(rna_chain, rna_chain.from_dict(rna_chain.to_dict()))

  def test_dna_dict_roundtrip(self):
    dna_chain = folding_input.DnaChain(
        id='A',
        sequence='AGCTN',
        modifications=[('6OG', 1), ('6MA', 2)],
        description='DNA chain',
    )
    self.assertEqual(dna_chain, dna_chain.from_dict(dna_chain.to_dict()))

  def test_ligand_dict_roundtrip(self):
    ligand = folding_input.Ligand(id='A', ccd_ids=['NAG', 'NAG'])
    self.assertEqual(ligand, ligand.from_dict(ligand.to_dict()))

  def test_from_json(self):
    test_json = data.Data(resources.ROOT / 'common/test_data/').load(
        'alphafold_input.json'
    )

    af_input = folding_input.Input.from_json(test_json)

    test_template_mmcif = data.Data(resources.ROOT / 'common/test_data/').load(
        'test_template.mmcif'
    )

    exp_input = folding_input.Input(
        name='Test Fold Job Number One',
        chains=[
            folding_input.ProteinChain(
                id='A',
                sequence='PREACHINGS',
                ptms=[('HY3', 1), ('P1L', 5)],
                unpaired_msa='>query\nPREACHINGS\n>unpaired seq 1\nPREA--INGS',
                paired_msa='>query\nPREACHINGS\n>paired seq 1\n--EACHINGS',
                templates=None,
            ),
            folding_input.ProteinChain(
                id='AA',
                sequence='REACHER',
                ptms=[],
                unpaired_msa=None,
                paired_msa=None,
                templates=[
                    folding_input.Template(
                        mmcif=test_template_mmcif,
                        query_to_template_map={
                            0: 0,
                            1: 1,
                            2: 2,
                            4: 3,
                            5: 4,
                            6: 8,
                        },
                    )
                ],
            ),
            folding_input.DnaChain(
                id='C',
                sequence='GATTACA',
                modifications=[('6OG', 1), ('6MA', 2)],
            ),
            folding_input.DnaChain(
                id='D',
                sequence='TGTAATC',
                modifications=[],
            ),
            folding_input.RnaChain(
                id='E',
                sequence='GUAC',
                modifications=[('2MG', 1), ('5MC', 4)],
                unpaired_msa='>query\nGUAC\n>unpaired seq 1\nGUACaa',
            ),
            folding_input.Ligand(id='F', ccd_ids=['ATP']),
            folding_input.Ligand(id='G', ccd_ids=['HEM']),
            folding_input.Ligand(id='H', ccd_ids=['HEM']),
            folding_input.Ligand(id='I', ccd_ids=['MG']),
            folding_input.Ligand(id='II', ccd_ids=['MG']),
            folding_input.Ligand(id='J', ccd_ids=['NAG', 'FUC']),
            folding_input.Ligand(id='JJ', ccd_ids=['NA']),
            folding_input.Ligand(id='X', ccd_ids=['NA']),
            folding_input.Ligand(id='Y', ccd_ids=['NA']),
            folding_input.Ligand(
                id='Z',
                smiles='c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO[P@@](=O)(O)O[P@](=O)(O)OP(=O)(O)O)O)O)N',
            ),
        ],
        rng_seeds=[10, 42],
        bonded_atom_pairs=[
            (('A', 1, 'CA'), ('AA', 1, 'CA')),
            (('A', 1, 'CA'), ('G', 1, 'CHA')),
            (('J', 1, 'O6'), ('J', 2, 'C1')),
        ],
    )

    self.assertEqual(af_input, exp_input)

  @parameterized.named_parameters(
      ('no_compression', lambda d: d),
      ('gzip_compression', gzip.compress),
      ('xz_compression', lzma.compress),
      ('zstd_compression', zstd.compress),
  )
  def test_from_json_with_external_data(self, compress_fn):
    unpaired_msa = '>query\nDEEPMIND'
    paired_msa = '>query\nDEEPMIND'
    template_mmcif = data.Data(resources.ROOT / 'common/test_data/').load(
        'test_template.mmcif'
    )
    unpaired_rna_msa = '>query\nGGG'
    paired_msa_bytes = compress_fn(paired_msa.encode('utf-8'))
    unpaired_msa_bytes = compress_fn(unpaired_msa.encode('utf-8'))
    template_mmcif_bytes = compress_fn(template_mmcif.encode('utf-8'))
    unpaired_rna_msa_bytes = compress_fn(unpaired_rna_msa.encode('utf-8'))

    unpaired_msa_path = self.create_tempfile(
        'unpaired_msa.a3m', content=unpaired_msa_bytes, mode='wb'
    ).full_path
    paired_msa_path = self.create_tempfile(
        'paired_msa.a3m', content=paired_msa_bytes, mode='wb'
    ).full_path
    template_mmcif_path = self.create_tempfile(
        'test_template.mmcif', content=template_mmcif_bytes, mode='wb'
    ).full_path
    unpaired_rna_msa_path = self.create_tempfile(
        'unpaired_rna_msa.a3m', content=unpaired_rna_msa_bytes, mode='wb'
    ).full_path

    test_json = json.dumps({
        'name': 'test_input',
        'modelSeeds': [1],
        'sequences': [
            {
                'protein': {
                    'id': 'A',
                    'sequence': 'DEEPMIND',
                    'unpairedMsaPath': unpaired_msa_path,
                    'pairedMsaPath': paired_msa_path,
                    'templates': [{
                        'mmcifPath': template_mmcif_path,
                        'queryIndices': [0, 1, 2, 3, 4, 5],
                        'templateIndices': [0, 1, 2, 3, 4, 5],
                    }],
                }
            },
            {
                'rna': {
                    'id': 'B',
                    'sequence': 'GGG',
                    'unpairedMsaPath': unpaired_rna_msa_path,
                }
            },
        ],
        'dialect': folding_input.JSON_DIALECT,
        'version': folding_input.JSON_VERSION,
    })

    af_input = folding_input.Input.from_json(test_json)
    self.assertEqual(af_input.protein_chains[0].unpaired_msa, unpaired_msa)
    self.assertEqual(af_input.protein_chains[0].paired_msa, paired_msa)
    self.assertIsNotNone(af_input.protein_chains[0].templates)
    self.assertEqual(
        af_input.protein_chains[0].templates[0].mmcif, template_mmcif
    )

  def test_from_json_with_external_data_relative_path(self):
    unpaired_msa = '>query\nDEEPMIND'

    # Write the MSA file to a temporary directory.
    unpaired_msa_path = self.create_tempfile(
        'unpaired_msa.a3m', content=unpaired_msa, mode='wt'
    ).full_path
    # Now pretend that the JSON file is in the same directory as the MSA file.
    json_path = unpaired_msa_path.replace('unpaired_msa.a3m', 'input.json')
    unpaired_msa_path = 'unpaired_msa.a3m'

    test_json = json.dumps({
        'name': 'test_input',
        'modelSeeds': [1],
        'sequences': [
            {
                'protein': {
                    'id': 'A',
                    'sequence': 'DEEPMIND',
                    'unpairedMsaPath': unpaired_msa_path,
                }
            },
        ],
        'dialect': folding_input.JSON_DIALECT,
        'version': folding_input.JSON_VERSION,
    })

    af_input = folding_input.Input.from_json(
        test_json, json_path=pathlib.Path(json_path)
    )
    self.assertEqual(af_input.protein_chains[0].unpaired_msa, unpaired_msa)

  def test_from_json_external_and_internal_data(self):
    base_json = {
        'name': 'test_input',
        'modelSeeds': [1],
        'sequences': [
            {
                'protein': {
                    'id': 'A',
                    'sequence': 'DEEPMIND',
                    'unpairedMsa': '>query\nDEEPMIND',
                    'pairedMsa': '>query\nDEEPMIND',
                    'templates': [{
                        'mmcif': 'data_TEST',
                        'queryIndices': [0, 1, 2, 3, 4, 5],
                        'templateIndices': [0, 1, 2, 3, 4, 5],
                    }],
                },
            },
            {
                'rna': {
                    'id': 'B',
                    'sequence': 'GGG',
                    'unpairedMsa': '>query\nGGG',
                },
            },
        ],
        'dialect': folding_input.JSON_DIALECT,
        'version': folding_input.JSON_VERSION,
    }

    with self.subTest('protein_unpaired_msa'):
      test_json = copy.deepcopy(base_json)
      test_json['sequences'][0]['protein']['unpairedMsaPath'] = 'unpaired.a3m'
      with self.assertRaisesRegex(ValueError, 'unpairedMsa/unpairedMsaPath'):
        folding_input.Input.from_json(json.dumps(test_json))

    with self.subTest('protein_paired_msa'):
      test_json = copy.deepcopy(base_json)
      test_json['sequences'][0]['protein']['pairedMsaPath'] = 'paired.a3m'
      with self.assertRaisesRegex(ValueError, 'pairedMsa/pairedMsaPath'):
        folding_input.Input.from_json(json.dumps(test_json))

    with self.subTest('protein_templates'):
      test_json = copy.deepcopy(base_json)
      test_json['sequences'][0]['protein']['templates'][0][  # pyrefly: ignore[unsupported-operation]
          'mmcifPath'
      ] = 'template.cif'
      with self.assertRaisesRegex(ValueError, 'mmcif/mmcifPath'):
        folding_input.Input.from_json(json.dumps(test_json))

    with self.subTest('rna_unpaired_msa'):
      test_json = copy.deepcopy(base_json)
      test_json['sequences'][1]['rna']['unpairedMsaPath'] = 'unpaired.a3m'
      with self.assertRaisesRegex(ValueError, 'unpairedMsa/unpairedMsaPath'):
        folding_input.Input.from_json(json.dumps(test_json))

  def test_fill_missing_fields_protein(self):
    protein_chain = folding_input.ProteinChain(
        id='A', sequence='ACDE', ptms=[('UNK', 1)]
    )
    self.assertEqual(
        protein_chain.fill_missing_fields(),
        folding_input.ProteinChain(
            id='A',
            sequence='ACDE',
            ptms=[('UNK', 1)],
            description=None,
            unpaired_msa='',
            paired_msa='',
            templates=[],
        ),
    )

  def test_fill_missing_fields_rna(self):
    rna_chain = folding_input.RnaChain(
        id='A', sequence='ACDE', modifications=[('UNK', 1)]
    )
    self.assertEqual(
        rna_chain.fill_missing_fields(),
        folding_input.RnaChain(
            id='A',
            sequence='ACDE',
            modifications=[('UNK', 1)],
            description=None,
            unpaired_msa='',
        ),
    )

  @parameterized.named_parameters(
      ('fill_missing_fields', True), ('do_not_fill_missing_fields', False)
  )
  def test_to_json(self, fill_missing_fields):
    test_json = data.Data(resources.ROOT / 'common/test_data/').load(
        'alphafold_input.json'
    )
    af_input = folding_input.Input.from_json(test_json)
    if fill_missing_fields:
      af_input = af_input.fill_missing_fields()
    af_input_json = af_input.to_json()
    af_round_trip = folding_input.Input.from_json(af_input_json)
    with self.subTest('round_trip'):
      self.assertEqual(af_input, af_round_trip)
    with self.subTest('no_newlines_in_template_indices'):
      self.assertIn('"queryIndices": [0, 1, 2, 4, 5, 6]', af_input_json)
      self.assertIn('"templateIndices": [0, 1, 2, 3, 4, 8]', af_input_json)

  def test_with_multiple_seeds(self):
    af_input = folding_input.Input(
        name='test_input',
        chains=[folding_input.ProteinChain(id='A', sequence='ACDE', ptms=[])],
        rng_seeds=[1337],
    )
    af_input_with_multiple_seeds = af_input.with_multiple_seeds(3)
    self.assertEqual(af_input_with_multiple_seeds.rng_seeds, (1337, 1338, 1339))

  def test_with_multiple_seeds_with_too_many_existing_seeds(self):
    af_input = folding_input.Input(
        name='test_input',
        chains=[folding_input.ProteinChain(id='A', sequence='ACDE', ptms=[])],
        rng_seeds=[1, 2],
    )
    with self.assertRaises(ValueError):
      af_input.with_multiple_seeds(3)

  @parameterized.parameters(0, 1)
  def test_with_multiple_seeds_bad_number_of_seeds(self, num_seeds):
    af_input = folding_input.Input(
        name='test_input',
        chains=[folding_input.ProteinChain(id='A', sequence='ACDE', ptms=[])],
        rng_seeds=[1],
    )
    with self.assertRaises(ValueError):
      af_input.with_multiple_seeds(num_seeds)

  def test_to_json_sequence_deduplication(self):
    af_input = {
        'name': 'test_input',
        'modelSeeds': [1337],
        'sequences': [
            {'protein': {'id': ['PA', 'PB'], 'sequence': 'RRR'}},
            {'protein': {'id': 'PC', 'sequence': 'RRR'}},
            {'rna': {'id': ['RA', 'RB'], 'sequence': 'UUU'}},
            {'rna': {'id': 'RC', 'sequence': 'UUU'}},
            {'dna': {'id': ['DA', 'DB'], 'sequence': 'AAA'}},
            {'dna': {'id': 'DC', 'sequence': 'AAA'}},
            {'ligand': {'id': ['A', 'B'], 'ccdCodes': ['GLY']}},
            {'ligand': {'id': ['C', 'D'], 'ccdCodes': ['GLY']}},
            {'ligand': {'id': 'E', 'ccdCodes': ['GLY']}},
        ],
        'bondedAtomPairs': None,
        'userCCD': None,
        'dialect': 'alphafold3',
        'version': 3,
    }
    round_trip_input = json.loads(
        folding_input.Input.from_json(json.dumps(af_input)).to_json()
    )
    self.assertLen(round_trip_input['sequences'], 4)
    self.assertLen(round_trip_input['sequences'][0]['protein']['id'], 3)
    self.assertLen(round_trip_input['sequences'][1]['rna']['id'], 3)
    self.assertLen(round_trip_input['sequences'][2]['dna']['id'], 3)
    self.assertLen(round_trip_input['sequences'][3]['ligand']['id'], 5)

  def test_to_json_sequence_deduplication_round_trip(self):
    af_input = {
        'name': 'test_input',
        'modelSeeds': [1337],
        'sequences': [
            {
                'protein': {
                    'id': ['A', 'B', 'C', 'D'],
                    'sequence': 'AAA',
                    'modifications': [],
                    'unpairedMsa': None,
                    'pairedMsa': None,
                    'templates': None,
                }
            },
            {
                'protein': {
                    'id': ['E', 'F'],
                    'sequence': 'EEE',
                    'modifications': [],
                    'unpairedMsa': None,
                    'pairedMsa': None,
                    'templates': None,
                }
            },
            {
                'rna': {
                    'id': 'G',
                    'sequence': 'GGG',
                    'modifications': [],
                    'unpairedMsa': None,
                }
            },
        ],
        'bondedAtomPairs': None,
        'userCCD': None,
        'dialect': 'alphafold3',
        'version': folding_input.JSON_VERSION,
    }
    round_trip_input = json.loads(
        folding_input.Input.from_json(json.dumps(af_input)).to_json()
    )
    self.assertEqual(af_input, round_trip_input)

  def test_from_json_some_chain_ids_given(self):
    test_json = data.Data(resources.ROOT / 'common/test_data/').load(
        'alphafold_input.json'
    )
    test_json = test_json.replace('"id": "A",', '')
    with self.assertRaises(ValueError):
      folding_input.Input.from_json(test_json)

  def test_from_json_chain_id_duplicates(self):
    test_json = data.Data(resources.ROOT / 'common/test_data/').load(
        'alphafold_input.json'
    )
    test_json = test_json.replace('"id": "A",', '"id": "AA",')
    with self.assertRaises(ValueError):
      folding_input.Input.from_json(test_json)

  def test_rna_paired_msa_unsupported(self):
    rna_input = {'rna': {'sequence': 'GG', 'pairedMsa': '>seq 1\nGG'}}
    with self.assertRaises(ValueError):
      folding_input.RnaChain.from_dict(rna_input)

  def test_from_mmcif_round_trip(self):
    fold_input = folding_input.Input(
        name='test_input',
        chains=[
            folding_input.ProteinChain(
                id='A',
                sequence='PREACHINGS',
                ptms=[('HY3', 1), ('P1L', 5)],
            ),
            folding_input.ProteinChain(
                id='B',
                sequence='REACHER',
                ptms=[],
            ),
            folding_input.RnaChain(
                id='C',
                sequence='GUAC',
                modifications=[('2MG', 1), ('5MC', 4)],
            ),
            folding_input.DnaChain(
                id='D',
                sequence='GATTACA',
                modifications=[('6OG', 1), ('6MA', 2)],
            ),
            folding_input.DnaChain(
                id='E',
                sequence='TGTAATC',
                modifications=[],
            ),
            folding_input.Ligand(id='F', ccd_ids=['ATP']),
            folding_input.Ligand(id='G', ccd_ids=['HEM']),
            folding_input.Ligand(id='H', ccd_ids=['MG']),
            folding_input.Ligand(id='I', ccd_ids=['NA']),
            folding_input.Ligand(
                id='J',
                smiles='c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO[P@@](=O)(O)O[P@](=O)(O)OP(=O)(O)O)O)O)N',
            ),
        ],
        # Warning: The seeds are wiped out when converting to mmCIF,
        # so just hack this by setting it to the same as mocked sampled seed
        # value.
        rng_seeds=[_SAMPLE_RNG_SEED],
        bonded_atom_pairs=[
            (('A', 1, 'CA'), ('B', 1, 'CA')),
            (('A', 1, 'CA'), ('F', 1, 'C4')),
        ],
    )

    mmcif_str = fold_input.to_structure(
        ccd=chemical_components.Ccd()
    ).to_mmcif()
    fold_input_from_mmcif = folding_input.Input.from_mmcif(
        mmcif_str, ccd=chemical_components.Ccd()
    )
    self.assertIsNotNone(fold_input_from_mmcif.bonded_atom_pairs)
    self.assertIsInstance(fold_input_from_mmcif.bonded_atom_pairs[0][0][1], int)
    self.assertEqual(fold_input_from_mmcif, fold_input)

  def test_to_structure(self):
    folding_input_with_ptms = folding_input.Input(
        name='test',
        chains=[
            folding_input.ProteinChain(
                id='A',
                sequence='KAGTAGT',
                ptms=[],
                unpaired_msa='',
                paired_msa='',
                templates=[],
            ),
            folding_input.ProteinChain(
                id='B',
                sequence='CKTCGSCGT',
                ptms=[('ALY', 2), ('SEP', 6)],
                unpaired_msa='',
                paired_msa='',
                templates=[],
            ),
            folding_input.RnaChain(
                id='C',
                sequence='ACACU',
                modifications=[('4OC', 4)],
                unpaired_msa='',
            ),
            folding_input.DnaChain(
                id='D',
                sequence='AGTAGT',
                modifications=[],
            ),
            folding_input.DnaChain(
                id='E',
                sequence='TAGATAGA',
                modifications=[('8OG', 3)],
            ),
            folding_input.Ligand(id='F', ccd_ids=['ADP']),
            folding_input.Ligand(id='G', ccd_ids=['CLA']),
            folding_input.Ligand(id='H', ccd_ids=['MG']),
            folding_input.Ligand(id='I', ccd_ids=['NAG', 'FUC']),
        ],
        bonded_atom_pairs=[(('I', 1, 'O6'), ('I', 2, 'C1'))],
        rng_seeds=[42],
    )

    struc = folding_input_with_ptms.to_structure(ccd=chemical_components.Ccd())
    self.assertEqual(struc.name, 'test')
    self.assertEqual(struc.num_chains, 9)
    self.assertEqual(struc.filter_to_entity_type(protein=True).num_chains, 2)
    self.assertEqual(struc.filter_to_entity_type(rna=True).num_chains, 1)
    self.assertEqual(struc.filter_to_entity_type(dna=True).num_chains, 2)
    self.assertEqual(struc.filter_to_entity_type(ligand=True).num_chains, 4)

    self.assertCountEqual(
        struc.chain_res_name_sequence().values(),
        [
            ('LYS', 'ALA', 'GLY', 'THR', 'ALA', 'GLY', 'THR'),
            ('CYS', 'ALY', 'THR', 'CYS', 'GLY', 'SEP', 'CYS', 'GLY', 'THR'),
            ('A', 'C', 'A', '4OC', 'U'),
            ('DA', 'DG', 'DT', 'DA', 'DG', 'DT'),
            ('DT', 'DA', '8OG', 'DA', 'DT', 'DA', 'DG', 'DA'),
            ('ADP',),
            ('CLA',),
            ('MG',),
            ('NAG', 'FUC'),
        ],
    )
    self.assertEqual(
        struc.chains_table.type.tolist(),
        [mmcif_names.PROTEIN_CHAIN] * 2
        + [mmcif_names.RNA_CHAIN] * 1
        + [mmcif_names.DNA_CHAIN] * 2
        + [mmcif_names.NON_POLYMER_CHAIN] * 3
        + [mmcif_names.BRANCHED_CHAIN],
    )
    self.assertLen(struc.bonds_table, 1)
    self.assertEqual(
        struc.bonds_table.type.tolist(), [mmcif_names.COVALENT_BOND]
    )

  @parameterized.named_parameters(
      dict(
          testcase_name='first_bond_first_atom_bad',
          bonded_atom_pairs=[(('P', 2, 'UNK'), ('L', 1, 'C1'))],
          component_name='ASN',
          chain_id='P',
          res_id=2,
      ),
      dict(
          testcase_name='first_bond_second_atom_bad',
          bonded_atom_pairs=[(('P', 1, 'ND2'), ('L', 1, 'UNK'))],
          component_name='NAG',
          chain_id='L',
          res_id=1,
      ),
      dict(
          testcase_name='first_bond_both_atoms_bad',
          bonded_atom_pairs=[(('P', 1, 'UNK'), ('L', 1, 'UNK2'))],
          component_name='ASN',
          chain_id='P',
          res_id=1,
      ),
      dict(
          testcase_name='second_bond_second_atom_bad',
          bonded_atom_pairs=[
              (('P', 1, 'ND2'), ('L', 1, 'C1')),
              (('L', 1, 'O4'), ('L', 2, 'UNK')),
          ],
          component_name='NAG',
          chain_id='L',
          res_id=2,
      ),
  )
  def test_to_structure_bad_bonds(
      self, bonded_atom_pairs, component_name, chain_id, res_id
  ):
    fold_input = folding_input.Input(
        name='bad_bonds_test',
        chains=[
            folding_input.ProteinChain(id='P', sequence='NNNNN', ptms=[]),
            folding_input.Ligand(id='L', ccd_ids=['NAG', 'NAG']),
        ],
        bonded_atom_pairs=bonded_atom_pairs,
        rng_seeds=[1],
    )
    with self.assertRaisesRegex(
        ValueError,
        f'Bonded atom "UNK" was not found in the list of atoms of the chemical'
        rf' component {component_name}.+\(chain_id={chain_id}, res_id={res_id},'
        r' atom_name=UNK\).+',
    ):
      fold_input.to_structure(ccd=chemical_components.Ccd())

  def test_bonds_unset_chain_ids(self):
    test_json = json.dumps(
        {
            'name': 'test_input',
            'modelSeeds': [4, 2],
            'bondedAtomPairs': [[['F', 1, 'C4'], ['G', 1, 'CHA']]],
            'sequences': [
                {'ligand': {'ccdCodes': ['ATP']}},
                {'ligand': {'ccdCodes': ['HEM']}},
            ],
            'dialect': folding_input.JSON_DIALECT,
            'version': folding_input.JSON_VERSION,
        },
        indent=2,
    )
    with self.assertRaisesRegex(ValueError, 'unset IDs'):
      folding_input.Input.from_json(test_json)

  def test_bonds_invalid_chain_ids(self):
    test_json = json.dumps(
        {
            'name': 'test_input',
            'modelSeeds': [4, 2],
            'bondedAtomPairs': [[['A', 1, 'C4'], ['G', 1, 'CHA']]],
            'sequences': [
                {'ligand': {'id': 'F', 'ccdCodes': ['ATP']}},
                {'ligand': {'id': 'G', 'ccdCodes': ['HEM']}},
            ],
            'dialect': folding_input.JSON_DIALECT,
            'version': folding_input.JSON_VERSION,
        },
        indent=2,
    )
    with self.assertRaisesRegex(ValueError, 'Invalid chain ID'):
      folding_input.Input.from_json(test_json)

  def test_bonds_invalid_residue_ids(self):
    test_json = json.dumps(
        {
            'name': 'test_input',
            'modelSeeds': [4, 2],
            'bondedAtomPairs': [[['F', 1, 'C4'], ['G', 2, 'CHA']]],
            'sequences': [
                {'ligand': {'id': 'F', 'ccdCodes': ['ATP']}},
                {'ligand': {'id': 'G', 'ccdCodes': ['HEM']}},
            ],
            'dialect': folding_input.JSON_DIALECT,
            'version': folding_input.JSON_VERSION,
        },
        indent=2,
    )
    with self.assertRaisesRegex(ValueError, 'Invalid residue ID'):
      folding_input.Input.from_json(test_json)

  def test_bonds_not_unique(self):
    test_json = json.dumps(
        {
            'name': 'test_input',
            'modelSeeds': [4, 2],
            'bondedAtomPairs': [
                [['F', 1, 'C4'], ['F', 1, 'CHA']],
                [['F', 1, 'C4'], ['F', 1, 'CHA']],
            ],
            'sequences': [
                {'ligand': {'id': 'F', 'ccdCodes': ['ATP']}},
                {'ligand': {'id': 'G', 'ccdCodes': ['HEM']}},
            ],
            'dialect': folding_input.JSON_DIALECT,
            'version': folding_input.JSON_VERSION,
        },
        indent=2,
    )
    with self.assertRaisesRegex(ValueError, 'Bonds are not unique'):
      folding_input.Input.from_json(test_json)

  def test_bonds_bad_format(self):
    test_json = json.dumps(
        {
            'name': 'test_input',
            'modelSeeds': [4, 2],
            'bondedAtomPairs': [[['F', 'ATP', 'C4'], ['G', 'HEM', 'CHA']]],
            'sequences': [
                {'ligand': {'id': 'F', 'ccdCodes': ['ATP']}},
                {'ligand': {'id': 'G', 'ccdCodes': ['HEM']}},
            ],
            'dialect': folding_input.JSON_DIALECT,
            'version': folding_input.JSON_VERSION,
        },
        indent=2,
    )
    with self.assertRaisesRegex(ValueError, 'must have 3 components'):
      folding_input.Input.from_json(test_json)

  def test_bonds_smiles(self):
    test_json = json.dumps(
        {
            'name': 'test_input',
            'modelSeeds': [4, 2],
            'bondedAtomPairs': [[['F', 1, 'C4'], ['G', 1, 'CHA']]],
            'sequences': [
                {'ligand': {'id': 'F', 'ccdCodes': ['ATP']}},
                {'ligand': {'id': 'G', 'smiles': 'C=C'}},
            ],
            'dialect': folding_input.JSON_DIALECT,
            'version': folding_input.JSON_VERSION,
        },
        indent=2,
    )
    with self.assertRaisesRegex(ValueError, 'unsupported SMILES ligand G'):
      folding_input.Input.from_json(test_json)

  def test_user_ccd_single_component(self):
    ccd = chemical_components.Ccd()
    exp_input = folding_input.Input(
        name='Test Fold Job Number One',
        chains=[folding_input.ProteinChain(id='A', sequence='RRR', ptms=[])],
        rng_seeds=[10],
        user_ccd=cif_dict.CifDict(ccd['ARG']).to_string(),
    )
    actual_input = folding_input.Input.from_json(exp_input.to_json())
    self.assertEqual(actual_input, exp_input)

  def test_user_ccd_multiple_components(self):
    ccd = chemical_components.Ccd()
    arg_component = cif_dict.CifDict(ccd['ARG'])
    gly_component = cif_dict.CifDict(ccd['GLY'])
    exp_input = folding_input.Input(
        name='Test Fold Job Number One',
        chains=[folding_input.ProteinChain(id='A', sequence='RRR', ptms=[])],
        rng_seeds=[10],
        user_ccd=arg_component.to_string() + gly_component.to_string(),
    )
    actual_input = folding_input.Input.from_json(exp_input.to_json())
    self.assertEqual(actual_input, exp_input)

  def test_user_ccd_path(self):
    ccd = chemical_components.Ccd()
    arg_cif = cif_dict.CifDict(ccd['ARG']).to_string()

    user_ccd_path = self.create_tempfile(
        'user_ccd.cif', content=arg_cif, mode='wt'
    ).full_path

    exp_input = folding_input.Input(
        name='Test Fold Job Number One',
        chains=[folding_input.ProteinChain(id='A', sequence='RRR', ptms=[])],
        rng_seeds=[10],
        user_ccd=arg_cif,
    )

    input_json = json.loads(exp_input.to_json())
    assert isinstance(input_json, dict)
    del input_json['userCCD']
    input_json['userCCDPath'] = user_ccd_path
    actual_input = folding_input.Input.from_json(json.dumps(input_json))
    self.assertEqual(actual_input, exp_input)

  def test_user_ccd_missing_keys(self):
    ccd = chemical_components.Ccd()
    bad_arg_component = dict(ccd['ARG'].items())
    gly_component = cif_dict.CifDict(ccd['GLY'])
    del bad_arg_component['_chem_comp_atom.pdbx_model_Cartn_x_ideal']
    bad_arg_component = cif_dict.CifDict(bad_arg_component)
    with self.assertRaisesRegex(
        ValueError,
        "ARG .* these keys: {'_chem_comp_atom.pdbx_model_Cartn_x_ideal'}",
    ):
      folding_input.Input(
          name='Test Fold Job Number One',
          chains=[folding_input.ProteinChain(id='A', sequence='P', ptms=[])],
          rng_seeds=[10],
          user_ccd=bad_arg_component.to_string() + gly_component.to_string(),
      )

  @parameterized.named_parameters(
      dict(
          testcase_name='with_dialect_and_version',
          test_json={
              'name': 'test_input',
              'sequences': [{'ligand': {'id': 'F', 'ccdCodes': ['ATP']}}],
              'dialect': folding_input.JSON_DIALECT,
              'version': folding_input.JSON_VERSION,
          },
          error_message=None,
      ),
      dict(
          testcase_name='unsupported_dialect',
          test_json={
              'name': 'test_input',
              'sequences': [{'ligand': {'id': 'F', 'ccdCodes': ['ATP']}}],
              'dialect': 'alphafold',
              'version': folding_input.JSON_VERSION,
          },
          error_message=(
              'AlphaFold 3 input JSON has unsupported dialect: alphafold,'
              ' expected alphafold3.'
          ),
      ),
      dict(
          testcase_name='unsupported_version',
          test_json={
              'name': 'test_input',
              'sequences': [{'ligand': {'id': 'F', 'ccdCodes': ['ATP']}}],
              'dialect': folding_input.JSON_DIALECT,
              'version': 50,
          },
          error_message=(
              'AlphaFold 3 input JSON has unsupported version: 50, expected one'
              f' of {folding_input.JSON_VERSIONS}.'
          ),
      ),
      dict(
          testcase_name='missing_dialect',
          test_json={
              'name': 'test_input',
              'sequences': [{'ligand': {'id': 'F', 'ccdCodes': ['ATP']}}],
              'version': folding_input.JSON_VERSION,
          },
          error_message=(
              'AlphaFold 3 input JSON must contain `dialect` and `version`'
              ' fields.'
          ),
      ),
      dict(
          testcase_name='missing_version',
          test_json={
              'name': 'test_input',
              'sequences': [{'ligand': {'id': 'F', 'ccdCodes': ['ATP']}}],
              'dialect': folding_input.JSON_DIALECT,
          },
          error_message=(
              'AlphaFold 3 input JSON must contain `dialect` and `version`'
              ' fields.'
          ),
      ),
  )
  def test_version_and_dialect(self, test_json: str, error_message: str | None):
    if error_message is not None:
      with self.assertRaisesWithLiteralMatch(ValueError, error_message):
        folding_input.Input.from_json(json.dumps(test_json, indent=2))

  @parameterized.named_parameters(
      dict(
          testcase_name='protein_simple',
          json_dict={
              'sequence': 'REACHER',
              'count': 1,
          },
          seq_id='A',
          expected_protein_chain=folding_input.ProteinChain(
              id='A', sequence='REACHER', ptms=[]
          ),
      ),
      dict(
          testcase_name='protein_modifications',
          json_dict={
              'sequence': 'PREACHINGS',
              'modifications': [
                  {'ptmType': 'CCD_HY3', 'ptmPosition': 1},
                  {'ptmType': 'CCD_P1L', 'ptmPosition': 5},
              ],
              'count': 1,
          },
          seq_id='A',
          expected_protein_chain=folding_input.ProteinChain(
              id='A',
              sequence='PREACHINGS',
              ptms=[('HY3', 1), ('P1L', 5)],
          ),
      ),
  )
  def test_protein_chain_from_alphafoldserver_dict(
      self,
      json_dict: Mapping[str, Any],
      seq_id: str,
      expected_protein_chain: folding_input.ProteinChain,
  ):
    protein_chain = folding_input.ProteinChain.from_alphafoldserver_dict(
        json_dict, seq_id
    )
    self.assertEqual(protein_chain, expected_protein_chain)

  def test_rna_chain_from_alphafoldserver_dict(self):
    json_dict = {
        'sequence': 'GUAC',
        'modifications': [
            {'modificationType': 'CCD_2MG', 'basePosition': 1},
            {'modificationType': 'CCD_5MC', 'basePosition': 4},
        ],
        'count': 1,
    }
    seq_id = 'A'
    expected_rna_chain = folding_input.RnaChain(
        id=seq_id, sequence='GUAC', modifications=[('2MG', 1), ('5MC', 4)]
    )
    rna_chain = folding_input.RnaChain.from_alphafoldserver_dict(
        json_dict, seq_id
    )
    self.assertEqual(rna_chain, expected_rna_chain)

  def test_dna_chain_from_alphafoldserver_dict(self):
    json_dict = {
        'sequence': 'AGTAGT',
        'modifications': [
            {'modificationType': 'CCD_6OG', 'basePosition': 1},
            {'modificationType': 'CCD_6MA', 'basePosition': 2},
        ],
        'count': 1,
    }
    seq_id = 'A'
    expected_dna_chain = folding_input.DnaChain(
        id=seq_id, sequence='AGTAGT', modifications=[('6OG', 1), ('6MA', 2)]
    )
    dna_chain = folding_input.DnaChain.from_alphafoldserver_dict(
        json_dict, seq_id
    )
    self.assertEqual(dna_chain, expected_dna_chain)

  @parameterized.named_parameters(
      dict(
          testcase_name='ligand',
          json_dict={
              'ligand': 'CCD_ATP',
              'count': 1,
          },
          seq_id='A',
          expected_ligand=folding_input.Ligand(id='A', ccd_ids=['ATP']),
      ),
      dict(
          testcase_name='ion',
          json_dict={
              'ion': 'MG',
              'count': 1,
          },
          seq_id='A',
          expected_ligand=folding_input.Ligand(id='A', ccd_ids=['MG']),
      ),
  )
  def test_ligand_from_alphafoldserver_dict(
      self,
      json_dict: Mapping[str, Any],
      seq_id: str,
      expected_ligand: folding_input.Ligand,
  ):
    ligand = folding_input.Ligand.from_alphafoldserver_dict(json_dict, seq_id)
    self.assertEqual(ligand, expected_ligand)

  @parameterized.named_parameters(
      dict(
          testcase_name='simple_no_glycans',
          test_json={
              'name': 'test 1',
              'modelSeeds': ['2', '3'],
              'sequences': [
                  {
                      'proteinChain': {
                          'sequence': 'PREACHER',
                          'modifications': [
                              {'ptmType': 'CCD_HY3', 'ptmPosition': 1},
                              {'ptmType': 'CCD_P1L', 'ptmPosition': 5},
                          ],
                          'count': 2,
                      }
                  },
                  {'rnaSequence': {'sequence': 'GUAC', 'count': 1}},
                  {
                      'rnaSequence': {
                          'sequence': 'AGCU',
                          'modifications': [
                              {
                                  'modificationType': 'CCD_2MG',
                                  'basePosition': 1,
                              },
                              {
                                  'modificationType': 'CCD_5MC',
                                  'basePosition': 4,
                              },
                          ],
                          'count': 1,
                      }
                  },
                  {'dnaSequence': {'sequence': 'TGTAATC', 'count': 1}},
                  {'ligand': {'ligand': 'CCD_ATP', 'count': 1}},
                  {'ligand': {'ligand': 'CCD_HEM', 'count': 2}},
                  {'ion': {'ion': 'NA', 'count': 3}},
              ],
          },
          expected_folding_input=folding_input.Input(
              name='test 1',
              chains=[
                  folding_input.ProteinChain(
                      id='A', sequence='PREACHER', ptms=[('HY3', 1), ('P1L', 5)]
                  ),
                  folding_input.ProteinChain(
                      id='B', sequence='PREACHER', ptms=[('HY3', 1), ('P1L', 5)]
                  ),
                  folding_input.RnaChain(
                      id='C', sequence='GUAC', modifications=[]
                  ),
                  folding_input.RnaChain(
                      id='D',
                      sequence='AGCU',
                      modifications=[('2MG', 1), ('5MC', 4)],
                  ),
                  folding_input.DnaChain(
                      id='E', sequence='TGTAATC', modifications=[]
                  ),
                  folding_input.Ligand(id='F', ccd_ids=['ATP']),
                  folding_input.Ligand(id='G', ccd_ids=['HEM']),
                  folding_input.Ligand(id='H', ccd_ids=['HEM']),
                  folding_input.Ligand(id='I', ccd_ids=['NA']),
                  folding_input.Ligand(id='J', ccd_ids=['NA']),
                  folding_input.Ligand(id='K', ccd_ids=['NA']),
              ],
              rng_seeds=[2, 3],
          ),
      ),
      dict(
          testcase_name='no_protein_templates',
          test_json={
              'name': 'test 1',
              'modelSeeds': ['2'],
              'sequences': [
                  {
                      'proteinChain': {
                          'sequence': 'AA',
                          'useStructureTemplate': True,
                      }
                  },
                  {
                      'proteinChain': {
                          'sequence': 'BB',
                          'useStructureTemplate': False,
                      }
                  },
              ],
          },
          expected_folding_input=folding_input.Input(
              name='test 1',
              chains=[
                  folding_input.ProteinChain(
                      id='A', sequence='AA', ptms=[], templates=None
                  ),
                  folding_input.ProteinChain(
                      id='B', sequence='BB', ptms=[], templates=[]
                  ),
              ],
              rng_seeds=[2],
          ),
      ),
      dict(
          testcase_name='error_on_glycans',
          test_json={
              'name': 'test 1',
              'modelSeeds': ['4', '5'],
              'sequences': [
                  {
                      'proteinChain': {
                          'sequence': 'PREACHER',
                          'modifications': [
                              {'ptmType': 'CCD_HY3', 'ptmPosition': 1},
                              {'ptmType': 'CCD_P1L', 'ptmPosition': 5},
                          ],
                          'glycans': [
                              {'residues': 'NAG(NAG)(BMA)', 'position': 8},
                              {'residues': 'BMA', 'position': 10},
                          ],
                          'count': 2,
                      }
                  },
              ],
          },
          error_message=(
              'Specifying glycans in the'
              f' `{folding_input.ALPHAFOLDSERVER_JSON_DIALECT}` format is not'
              ' supported.'
          ),
      ),
      dict(
          testcase_name='error_on_max_template_date',
          test_json={
              'name': 'test 1',
              'modelSeeds': ['4', '5'],
              'sequences': [
                  {
                      'proteinChain': {
                          'sequence': 'AAA',
                          'maxTemplateDate': '2024-08-01',
                      }
                  },
              ],
          },
          error_message=(
              'Specifying maxTemplateDate in the'
              f' `{folding_input.ALPHAFOLDSERVER_JSON_DIALECT}` format is not '
              'supported, use the --max_template_date flag instead.'
          ),
      ),
      dict(
          testcase_name='with_dialect_and_version',
          test_json={
              'name': 'test',
              'modelSeeds': [],
              'sequences': [
                  {'proteinChain': {'sequence': 'TEACHINGS', 'count': 1}},
                  {'dnaSequence': {'sequence': 'TAGGACA', 'count': 1}},
              ],
              'dialect': folding_input.ALPHAFOLDSERVER_JSON_DIALECT,
              'version': folding_input.ALPHAFOLDSERVER_JSON_VERSION,
          },
          expected_folding_input=folding_input.Input(
              name='test',
              chains=[
                  folding_input.ProteinChain(
                      id='A', sequence='TEACHINGS', ptms=[]
                  ),
                  folding_input.DnaChain(
                      id='B', sequence='TAGGACA', modifications=[]
                  ),
              ],
              rng_seeds=[_SAMPLE_RNG_SEED],
          ),
      ),
      dict(
          testcase_name='unsupported_dialect',
          test_json={
              'name': 'test',
              'modelSeeds': [],
              'sequences': [
                  {'proteinChain': {'sequence': 'TEACHINGS', 'count': 1}},
                  {'dnaSequence': {'sequence': 'TAGGACA', 'count': 1}},
              ],
              'dialect': 'alphafold3server',
              'version': folding_input.ALPHAFOLDSERVER_JSON_VERSION,
          },
          error_message=(
              'AlphaFold Server input JSON has unsupported dialect:'
              ' alphafold3server, expected alphafoldserver.'
          ),
      ),
      dict(
          testcase_name='unsupported_version',
          test_json={
              'name': 'test',
              'modelSeeds': [],
              'sequences': [
                  {'proteinChain': {'sequence': 'TEACHINGS', 'count': 1}},
                  {'dnaSequence': {'sequence': 'TAGGACA', 'count': 1}},
              ],
              'dialect': folding_input.ALPHAFOLDSERVER_JSON_DIALECT,
              'version': 2,
          },
          error_message=(
              'AlphaFold Server input JSON has unsupported version: 2,'
              ' expected 1.'
          ),
      ),
      dict(
          testcase_name='only_specify_dialect',
          test_json={
              'name': 'test',
              'modelSeeds': [],
              'sequences': [
                  {'proteinChain': {'sequence': 'TEACHINGS', 'count': 1}},
                  {'dnaSequence': {'sequence': 'TAGGACA', 'count': 1}},
              ],
              'dialect': folding_input.ALPHAFOLDSERVER_JSON_DIALECT,
          },
          error_message=(
              'AlphaFold Server input JSON must either contain both `dialect`'
              ' and `version` fields, or neither. If neither is specified, it'
              ' is assumed that'
              f' `dialect="{folding_input.ALPHAFOLDSERVER_JSON_DIALECT}"` and'
              f' `version="{folding_input.ALPHAFOLDSERVER_JSON_VERSION}"`.'
          ),
      ),
      dict(
          testcase_name='only_specify_version',
          test_json={
              'name': 'test',
              'modelSeeds': [],
              'sequences': [
                  {'proteinChain': {'sequence': 'TEACHINGS', 'count': 1}},
                  {'dnaSequence': {'sequence': 'TAGGACA', 'count': 1}},
              ],
              'version': 7,
          },
          error_message=(
              'AlphaFold Server input JSON must either contain both `dialect`'
              ' and `version` fields, or neither. If neither is specified, it'
              ' is assumed that'
              f' `dialect="{folding_input.ALPHAFOLDSERVER_JSON_DIALECT}"` and'
              f' `version="{folding_input.ALPHAFOLDSERVER_JSON_VERSION}"`.'
          ),
      ),
  )
  def test_from_alphafoldserver_fold_job(
      self,
      test_json: Mapping[str, Any],
      expected_folding_input: folding_input.Input | None = None,
      error_message: str | None = None,
  ):
    """Test loading AlphafoldServer JSON to Input."""
    if error_message is not None:
      with self.assertRaisesWithLiteralMatch(ValueError, error_message):
        folding_input.Input.from_alphafoldserver_fold_job(test_json)
    else:
      self.assertEqual(
          folding_input.Input.from_alphafoldserver_fold_job(test_json),
          expected_folding_input,
      )

  @parameterized.named_parameters(
      dict(
          testcase_name='alphafold_server',
          test_json=[
              {
                  'name': 'test0',
                  'modelSeeds': [],
                  'sequences': [
                      {'proteinChain': {'sequence': 'TEACHINGS', 'count': 1}},
                      {'dnaSequence': {'sequence': 'TAGGACA', 'count': 1}},
                  ],
                  'dialect': folding_input.ALPHAFOLDSERVER_JSON_DIALECT,
                  'version': folding_input.ALPHAFOLDSERVER_JSON_VERSION,
              },
              {
                  'name': 'test1',
                  'modelSeeds': [],
                  'sequences': [
                      {'proteinChain': {'sequence': 'PREACHINGS', 'count': 1}},
                      {'dnaSequence': {'sequence': 'GATTACA', 'count': 1}},
                  ],
                  'dialect': folding_input.ALPHAFOLDSERVER_JSON_DIALECT,
                  'version': folding_input.ALPHAFOLDSERVER_JSON_VERSION,
              },
          ],
          expected_fold_inputs=[
              folding_input.Input(
                  name='test0',
                  chains=[
                      folding_input.ProteinChain(
                          id='A', sequence='TEACHINGS', ptms=[]
                      ),
                      folding_input.DnaChain(
                          id='B', sequence='TAGGACA', modifications=[]
                      ),
                  ],
                  rng_seeds=[_SAMPLE_RNG_SEED],
              ),
              folding_input.Input(
                  name='test1',
                  chains=[
                      folding_input.ProteinChain(
                          id='A', sequence='PREACHINGS', ptms=[]
                      ),
                      folding_input.DnaChain(
                          id='B', sequence='GATTACA', modifications=[]
                      ),
                  ],
                  rng_seeds=[_SAMPLE_RNG_SEED],
              ),
          ],
      ),
      dict(
          testcase_name='alphafold3',
          test_json={
              'name': 'test_input',
              'modelSeeds': [1, 2, 3],
              'sequences': [
                  {'ligand': {'id': 'F', 'ccdCodes': ['ATP']}},
                  {'ligand': {'id': 'G', 'ccdCodes': ['HEM']}},
              ],
              'dialect': folding_input.JSON_DIALECT,
              'version': folding_input.JSON_VERSION,
          },
          expected_fold_inputs=[
              folding_input.Input(
                  name='test_input',
                  chains=[
                      folding_input.Ligand(id='F', ccd_ids=['ATP']),
                      folding_input.Ligand(id='G', ccd_ids=['HEM']),
                  ],
                  rng_seeds=[1, 2, 3],
              ),
          ],
      ),
  )
  def test_load_fold_inputs_from_path(
      self, test_json, expected_fold_inputs: Sequence[folding_input.Input]
  ):
    """Test loading fold inputs from a JSON string."""
    temp_file = self.create_tempfile()
    with open(temp_file.full_path, 'w') as f:
      f.write(json.dumps(test_json, indent=2))

    self.assertEqual(
        list(
            folding_input.load_fold_inputs_from_path(
                pathlib.Path(temp_file.full_path)
            )
        ),
        expected_fold_inputs,
    )

  @parameterized.named_parameters(
      dict(
          testcase_name='alphafold_server_not_list',
          test_json={
              'name': 'test0',
              'modelSeeds': [],
              'sequences': [
                  {'proteinChain': {'sequence': 'TEACHINGS', 'count': 1}},
              ],
              'dialect': folding_input.ALPHAFOLDSERVER_JSON_DIALECT,
              'version': folding_input.ALPHAFOLDSERVER_JSON_VERSION,
          },
          error_message='Failed.*AlphaFold 3 dialect.*',
      ),
      dict(
          testcase_name='alphafold3_json_as_list',
          test_json=[
              {
                  'name': 'test_input',
                  'modelSeeds': [1, 2, 3],
                  'sequences': [
                      {'ligand': {'id': 'F', 'ccdCodes': ['ATP']}},
                      {'ligand': {'id': 'G', 'ccdCodes': ['HEM']}},
                  ],
                  'dialect': folding_input.JSON_DIALECT,
                  'version': folding_input.JSON_VERSION,
              },
          ],
          error_message='Failed.*AlphaFold Server dialect.*',
      ),
  )
  def test_load_fold_inputs_from_path_error(self, test_json, error_message):
    """Test loading fold inputs from a JSON string."""
    temp_file = self.create_tempfile()
    with open(temp_file.full_path, 'w') as f:
      f.write(json.dumps(test_json, indent=2))

    with self.assertRaisesRegex(ValueError, error_message):
      list(
          folding_input.load_fold_inputs_from_path(
              pathlib.Path(temp_file.full_path)
          )
      )

  def test_load_fold_inputs_from_dir(self):
    """Test loading fold inputs from a directory."""
    temp_dir = self.create_tempdir()
    temp_dir.create_file(
        'test0.json',
        content=json.dumps([
            {
                'name': 'test0',
                'modelSeeds': [],
                'sequences': [
                    {'proteinChain': {'sequence': 'TEACHINGS', 'count': 1}},
                ],
                'dialect': folding_input.ALPHAFOLDSERVER_JSON_DIALECT,
                'version': folding_input.ALPHAFOLDSERVER_JSON_VERSION,
            },
            {
                'name': 'test1',
                'modelSeeds': [],
                'sequences': [
                    {'proteinChain': {'sequence': 'PREACHINGS', 'count': 1}},
                ],
                'dialect': folding_input.ALPHAFOLDSERVER_JSON_DIALECT,
                'version': folding_input.ALPHAFOLDSERVER_JSON_VERSION,
            },
        ]),
    )
    temp_dir.create_file(
        'test1.json',
        content=json.dumps({
            'name': 'test_input',
            'modelSeeds': [1, 2, 3],
            'sequences': [
                {'ligand': {'id': 'F', 'ccdCodes': ['ATP']}},
            ],
            'dialect': folding_input.JSON_DIALECT,
            'version': folding_input.JSON_VERSION,
        }),
    )
    temp_dir.create_file(
        'ignored_file.txt',
        content=json.dumps({'ignored_job': 'ignored'}),
    )
    self.assertEqual(
        [
            folding_input.Input(
                name='test0',
                chains=[
                    folding_input.ProteinChain(
                        id='A', sequence='TEACHINGS', ptms=[]
                    ),
                ],
                rng_seeds=[_SAMPLE_RNG_SEED],
            ),
            folding_input.Input(
                name='test1',
                chains=[
                    folding_input.ProteinChain(
                        id='A', sequence='PREACHINGS', ptms=[]
                    ),
                ],
                rng_seeds=[_SAMPLE_RNG_SEED],
            ),
            folding_input.Input(
                name='test_input',
                chains=[
                    folding_input.Ligand(id='F', ccd_ids=['ATP']),
                ],
                rng_seeds=[1, 2, 3],
            ),
        ],
        list(
            folding_input.load_fold_inputs_from_dir(
                pathlib.Path(temp_dir.full_path)
            )
        ),
    )

  @parameterized.named_parameters(
      dict(
          testcase_name='alphafold_server_no_rng_seeds',
          test_json={
              'name': 'test_input',
              'modelSeeds': [],
              'sequences': [
                  {'proteinChain': {'sequence': 'TEACHINGS', 'count': 1}},
              ],
              'dialect': folding_input.ALPHAFOLDSERVER_JSON_DIALECT,
              'version': folding_input.ALPHAFOLDSERVER_JSON_VERSION,
          },
          expected_folding_input=folding_input.Input(
              name='test_input',
              chains=[
                  folding_input.ProteinChain(
                      id='A', sequence='TEACHINGS', ptms=[]
                  )
              ],
              rng_seeds=[_SAMPLE_RNG_SEED],
          ),
      ),
      dict(
          testcase_name='alphafold_server_with_rng_seeds',
          test_json={
              'name': 'test_input',
              'modelSeeds': [1, 2, 3],
              'sequences': [
                  {'ligand': {'ligand': 'CCD_ATP', 'count': 1}},
              ],
              'dialect': folding_input.ALPHAFOLDSERVER_JSON_DIALECT,
              'version': folding_input.ALPHAFOLDSERVER_JSON_VERSION,
          },
          expected_folding_input=folding_input.Input(
              name='test_input',
              chains=[folding_input.Ligand(id='A', ccd_ids=['ATP'])],
              rng_seeds=[1, 2, 3],
          ),
      ),
  )
  def test_from_alphafoldserver_fold_job_rng_seeds(
      self,
      test_json: Mapping[str, Any],
      expected_folding_input: folding_input.Input,
  ):
    """Specifically checks behaviour with rng seeds."""
    self.assertEqual(
        folding_input.Input.from_alphafoldserver_fold_job(test_json),
        expected_folding_input,
    )

  def test_folding_input_requires_rng_seed(self):
    """Checks that folding input requires rng seeds."""
    with self.assertRaisesWithLiteralMatch(
        ValueError, 'Input must have at least one RNG seed.'
    ):
      folding_input.Input(
          name='test_input',
          chains=[folding_input.Ligand(id='F', ccd_ids=['ATP'])],
          rng_seeds=[],
      )

  @parameterized.named_parameters(
      ('empty', ''), ('blank', ' '), ('non_permitted', 'ř*')
  )
  def test_folding_input_valid_name(self, name):
    """Checks that folding input requires a name."""
    with self.assertRaisesRegex(ValueError, 'Input name must be non-empty'):
      folding_input.Input(
          name=name,
          chains=[folding_input.Ligand(id='F', ccd_ids=['ATP'])],
          rng_seeds=[1],
      )

  @parameterized.parameters(
      ('Hi1.-_', 'Hi1.-_'), ('Hi 1', 'Hi_1'), ('Seven(ty)', 'Seventy')
  )
  def test_folding_input_sanitised_name(self, name, exp_sanitised_name):
    fold_input = folding_input.Input(
        name=name,
        chains=[folding_input.Ligand(id='F', ccd_ids=['ATP'])],
        rng_seeds=[1],
    )
    self.assertEqual(fold_input.sanitised_name(), exp_sanitised_name)

  def test_folding_input_descriptions(self):
    """Checks round trip of protein/rna/dna/ligand descriptions."""
    fold_input = folding_input.Input(
        name='test_input',
        chains=[
            folding_input.ProteinChain(
                id='A', sequence='RRR', ptms=[], description='Protein desc'
            ),
            folding_input.ProteinChain(
                id='B', sequence='RRR', ptms=[]  # No description
            ),
            folding_input.RnaChain(
                id='C', sequence='AU', modifications=[], description='RNA desc'
            ),
            folding_input.RnaChain(
                id='D', sequence='AU', modifications=[]  # No description
            ),
            folding_input.DnaChain(
                id='E', sequence='AG', modifications=[], description='DNA desc'
            ),
            folding_input.DnaChain(
                id='F', sequence='AG', modifications=[]  # No description
            ),
            folding_input.Ligand(
                id='G', ccd_ids=['ATP'], description='Ligand desc'
            ),
            folding_input.Ligand(id='H', ccd_ids=['ATP']),  # No description
        ],
        rng_seeds=[1],
    )
    fold_input_json = fold_input.to_json()
    reconstructed_fold_input = folding_input.Input.from_json(fold_input_json)

    with self.subTest('descriptions_in_json'):
      self.assertIn('Protein desc', fold_input_json)
      self.assertIn('RNA desc', fold_input_json)
      self.assertIn('DNA desc', fold_input_json)
      self.assertIn('Ligand desc', fold_input_json)
    with self.subTest('input_reconstructed'):
      self.assertEqual(reconstructed_fold_input, fold_input)

  def test_input_json_requires_rng_seed(self):
    """Checks that input JSON requires rng seeds."""
    with self.assertRaisesWithLiteralMatch(
        ValueError,
        'AlphaFold 3 input JSON must specify at least one rng seed in'
        ' `modelSeeds`.',
    ):
      folding_input.Input.from_json(
          json.dumps({
              'name': 'test_input',
              'modelSeeds': [],
              'sequences': [{'ligand': {'id': 'G', 'ccdCodes': ['HEM']}}],
              'dialect': folding_input.JSON_DIALECT,
              'version': folding_input.JSON_VERSION,
          })
      )

  def test_ligand_invalid_smiles(self):
    with self.assertRaisesWithLiteralMatch(
        ValueError, 'Unable to make RDKit Mol from SMILES: invalid_smiles'
    ):
      folding_input.Input.from_json(
          json.dumps({
              'name': 'test_input',
              'modelSeeds': [0],
              'sequences': [
                  {'ligand': {'id': 'A', 'smiles': 'invalid_smiles'}},
              ],
              'dialect': folding_input.JSON_DIALECT,
              'version': folding_input.JSON_VERSION,
          })
      )

  def test_ligand_ccd_codes_string(self):
    """Checks that the CCD codes are not a string."""
    with self.assertRaisesRegex(
        ValueError, 'CCD codes must be a list of strings, got str instead'
    ):
      folding_input.Input.from_json(
          json.dumps({
              'name': 'test_input',
              'modelSeeds': [10],
              'sequences': [{'ligand': {'id': 'G', 'ccdCodes': 'HEM'}}],
              'dialect': folding_input.JSON_DIALECT,
              'version': folding_input.JSON_VERSION,
          })
      )

  def test_ligand_ccd_codes_and_smiles_set(self):
    with self.assertRaisesRegex(
        ValueError, 'Ligand cannot have both CCD code and SMILES set'
    ):
      folding_input.Input.from_json(
          json.dumps({
              'name': 'test_input',
              'modelSeeds': [10],
              'sequences': [
                  {'ligand': {'id': 'G', 'ccdCodes': ['HEM'], 'smiles': 'O=O'}}
              ],
              'dialect': folding_input.JSON_DIALECT,
              'version': folding_input.JSON_VERSION,
          })
      )

  @parameterized.parameters(
      {'unpaired_msa_path': 'unpaired.fa'},
      {'unpaired_msa_path': 'gs://my-bucket/unpaired.fa'},
  )
  def test_load_inputs_from_epath(self, unpaired_msa_path):
    input_json = json.dumps({
        'dialect': 'alphafold3',
        'version': 1,
        'name': 'test_gcs_job',
        'modelSeeds': [1],
        'sequences': [{
            'protein': {
                'id': 'A',
                'sequence': 'MAEVT',
                'unpairedMsaPath': unpaired_msa_path,
            }
        }],
    })
    msa_content = '>seq1\nMAEVT'

    def mock_exists(original_fn, path):
      if str(path).endswith((
          'my-bucket/alphafold_input.json',
          'my-bucket/unpaired.fa',
      )):
        return True
      return original_fn(path)

    def mock_open(original_fn, path, mode):
      if str(path).endswith('my-bucket/alphafold_input.json'):
        return io.StringIO(input_json)
      elif str(path).endswith('my-bucket/unpaired.fa'):
        return io.BytesIO(msa_content.encode())
      return original_fn(path, mode)

    with epath_testing.mock_epath(exists=mock_exists, open=mock_open):
      fold_inputs = list(
          folding_input.load_fold_inputs_from_path(
              'gs://my-bucket/alphafold_input.json'
          )
      )
      self.assertLen(fold_inputs, 1)
      self.assertEqual(fold_inputs[0].name, 'test_gcs_job')
      chain = fold_inputs[0].chains[0]
      assert isinstance(chain, folding_input.ProteinChain)
      self.assertEqual(chain.unpaired_msa, msa_content)

  def test_rna_modifications(self):
    rna_chain = folding_input.RnaChain(
        id='A', sequence='AU', modifications=[('2MG', 1)]
    )
    self.assertSequenceEqual(rna_chain.modifications, [('2MG', 1)])

  def test_dna_modifications(self):
    dna_chain = folding_input.DnaChain(
        id='A', sequence='AG', modifications=[('6OG', 1)]
    )
    self.assertSequenceEqual(dna_chain.modifications, [('6OG', 1)])


class SampleRandomSeedTest(parameterized.TestCase):
  """Separate test class for _sample_rng_seed() where it's not mocked."""

  def test_sample_random_seed(self):
    """Checks that random seed is sampled when not specified."""
    seeds = [folding_input._sample_rng_seed() for _ in range(100)]

    self.assertLessEqual(np.max(seeds), 2**32 - 1)
    self.assertGreaterEqual(np.min(seeds), 0)

    # It's too strong to check that all the seeds are unique, but we can check
    # that they are at least reasonably random.
    self.assertGreater(np.var(seeds), 0)


if __name__ == '__main__':
  absltest.main()
