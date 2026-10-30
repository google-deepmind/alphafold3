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
import dataclasses
import datetime
from typing import Any, Final
from unittest import mock

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3 import structure
from alphafold3.common import resources
from alphafold3.constants import atom_types
from alphafold3.constants import chemical_components
from alphafold3.constants import mmcif_names
from alphafold3.constants import residue_names
from alphafold3.data import msa
from alphafold3.data import msa_config
from alphafold3.data import structure_stores
from alphafold3.data import templates
from alphafold3.structure import test_utils
import chex
import jax
import numpy as np

import shutil

_HMMSEARCH_BINARY_PATH = shutil.which('hmmsearch')
_HMMBUILD_BINARY_PATH = shutil.which('hmmbuild')


_EXAMPLE_A3M = """
>abcd_A/1-8 [subseq from] length:8
ABCDEFGH
>abcd_B/1-16 [subseq from] length:16 insertion
ABCDEFiiiiiiiiGH
>efgh_C/1-4 [subseq from] length:4 deletion
AB----GH
>ijkl_A/1-4 [subseq from] length:4 deletion
AB----GH
>noex_Z/1-4 [subseq from] length:4 deletion
AB----GH
""".lstrip()


_EXAMPLE_SEQRES_DB = """
>abcd_A mol:protein length:8 UNOBTAININ
ABCDEFGH
>abcd_B mol:protein length:16 UNOBTAININ
ABCDEFIIIIIIIIGH
>efgh_C mol:protein length:4 DYPROSIN
ABGH
>ijkl_A mol:protein length:4 ELUDIN
ABGH
>noex_Z mol:protein length:4 NOEXIN
ABGH
""".lstrip()


_STRUCTURE_PATH = (
    resources.ROOT / 'test_data/miniature_databases/pdb_mmcif'
)

_TEMPLATES_V2_CONFIG = msa_config.TemplatesConfig(
    template_tool_config=msa_config.TemplateToolConfig(
        database_path=(
            str(
                resources.ROOT
                / 'test_data/miniature_databases/'
                / 'pdb_seqres_2022_09_28__subsampled_1000.fasta'
            )
        ),
        max_a3m_query_sequences=None,
        hmmsearch_config=msa_config.HmmsearchConfig(
            hmmsearch_binary_path=_HMMSEARCH_BINARY_PATH,
            hmmbuild_binary_path=_HMMBUILD_BINARY_PATH,
            n_cpu=8,
            filter_f1=0.1,
            filter_f2=0.1,
            filter_f3=0.1,
            e_value=100,
            inc_e=100,
            dom_e=100,
            incdom_e=100,
            alphabet='amino',
        ),
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
    ),
    filter_config=msa_config.TemplateFilterConfig.no_op_filter(),
)


_ARTIFICIAL_QUERY_SEQUENCE: Final[str] = 'ABCDEFGH'

_ARTIFICIAL_HITS_DATA = (
    templates.Hit(
        pdb_id='abcd',
        auth_chain_id='A',
        hmmsearch_sequence='ABCDEFGH',
        structure_sequence='ABCDEFGH',
        unresolved_res_indices=[],
        start_index=0,
        end_index=8,
        full_length=8,
        release_date=datetime.date(2010, 1, 1),
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
        query_sequence=_ARTIFICIAL_QUERY_SEQUENCE,
    ),
    templates.Hit(
        pdb_id='abcd',
        auth_chain_id='B',
        hmmsearch_sequence='ABCDEFiiiiiiiiGH',
        structure_sequence='ABCDEFIIIIIIIIGH',
        query_sequence=_ARTIFICIAL_QUERY_SEQUENCE,
        unresolved_res_indices=[],
        start_index=0,
        end_index=16,
        full_length=16,
        release_date=datetime.date(2010, 1, 1),
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
    ),
    templates.Hit(
        pdb_id='efgh',
        auth_chain_id='C',
        hmmsearch_sequence='AB----GH',
        structure_sequence='ABGH',
        query_sequence=_ARTIFICIAL_QUERY_SEQUENCE,
        unresolved_res_indices=[],
        start_index=0,
        end_index=4,
        full_length=4,
        release_date=datetime.date(2020, 1, 1),
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
    ),
    templates.Hit(
        pdb_id='ijkl',
        auth_chain_id='A',
        hmmsearch_sequence='AB----GH',
        structure_sequence='ABGH',
        query_sequence=_ARTIFICIAL_QUERY_SEQUENCE,
        unresolved_res_indices=[],
        start_index=0,
        end_index=4,
        full_length=4,
        release_date=datetime.date(2020, 1, 1),
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
    ),
)

_PDB_HITS_DATA: Sequence[templates.Hit] = (
    templates.Hit(
        pdb_id='5pti',
        auth_chain_id='',
        hmmsearch_sequence='A',
        structure_sequence='A',
        query_sequence=_ARTIFICIAL_QUERY_SEQUENCE,
        unresolved_res_indices=[],
        start_index=0,
        end_index=1,
        full_length=1,
        release_date=datetime.date(1984, 10, 29),
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
    ),
    templates.Hit(
        pdb_id='5pti',
        auth_chain_id='A',
        hmmsearch_sequence='A',
        structure_sequence='A',
        query_sequence=_ARTIFICIAL_QUERY_SEQUENCE,
        unresolved_res_indices=[],
        start_index=0,
        end_index=1,
        full_length=1,
        release_date=datetime.date(1984, 10, 29),
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
    ),
    templates.Hit(
        pdb_id='6ysc',
        auth_chain_id='B',
        hmmsearch_sequence='A',
        structure_sequence='A',
        query_sequence=_ARTIFICIAL_QUERY_SEQUENCE,
        unresolved_res_indices=[],
        start_index=0,
        end_index=1,
        full_length=1,
        release_date=datetime.date(2020, 10, 14),
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
    ),
    templates.Hit(
        pdb_id='7rye',
        auth_chain_id='J',
        hmmsearch_sequence='A',
        structure_sequence='A',
        query_sequence=_ARTIFICIAL_QUERY_SEQUENCE,
        unresolved_res_indices=[],
        start_index=0,
        end_index=1,
        full_length=1,
        release_date=datetime.date(2021, 11, 10),
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
    ),
)


class TemplatesV2Test(test_utils.StructureTestCase):

  def test_from_seq_and_a3m(self):
    query_sequence = 'ABCDEFGH'
    seqres_db = self.create_tempfile()
    with open(seqres_db.full_path, 'wt') as f:
      f.write(_EXAMPLE_SEQRES_DB)
    with mock.patch.object(
        templates,
        '_parse_hit_metadata',
        autospec=True,
        return_value=('2020-01-01', 'ABGH', []),
    ):
      found_templates = templates.Templates.from_seq_and_a3m(
          query_sequence=query_sequence,
          msa_a3m=f'>query\n{query_sequence}\n',
          database_path=seqres_db.full_path,
          hmmsearch_config=_TEMPLATES_V2_CONFIG.template_tool_config.hmmsearch_config,
          max_a3m_query_sequences=_TEMPLATES_V2_CONFIG.template_tool_config.max_a3m_query_sequences,
          max_template_date=datetime.date(3000, 1, 1),
          structure_store=structure_stores.StructureStore(_STRUCTURE_PATH),
      )
    self.assertNotEmpty(found_templates.hits)

  def test_from_hmmsearch_a3m(self):

    def mocked_parse_hit_metadata(
        structure_store, pdb_id: str, auth_chain_id: str
    ) -> tuple[str | None, str | None, Sequence[int] | None]:
      del structure_store
      if pdb_id == 'noex':
        return None, None, None
      elif pdb_id == 'abcd':
        if auth_chain_id == 'A':
          return '2010-01-01', 'ABCDEFGH', []
        elif auth_chain_id == 'B':
          return '2010-01-01', 'ABCDEFIIIIIIIIGH', []
        else:
          raise ValueError(f'Unknown {auth_chain_id=}.')
      elif pdb_id == 'efgh':
        return '2020-01-01', 'ABGH', []
      elif pdb_id == 'ijkl':
        return '2020-01-01', 'ABGH', []
      else:
        raise ValueError(f'Unknown {pdb_id=}.')

    with mock.patch.object(
        templates,
        '_parse_hit_metadata',
        new=mocked_parse_hit_metadata,
    ):
      found_templates = templates.Templates.from_hmmsearch_a3m(
          query_sequence=_ARTIFICIAL_QUERY_SEQUENCE,
          a3m=_EXAMPLE_A3M,
          max_template_date=datetime.date(3000, 1, 1),
          structure_store=structure_stores.StructureStore(_STRUCTURE_PATH),
      )

    for idx, (obs, exp) in enumerate(
        zip(found_templates.hits, _ARTIFICIAL_HITS_DATA, strict=True)
    ):
      with self.subTest(f'Hit{idx}'):
        self.assertEqual(dataclasses.asdict(obs), dataclasses.asdict(exp))

  def test_from_hmmsearch_empty_a3m(self):
    found_templates = templates.Templates.from_hmmsearch_a3m(
        query_sequence='EMPTY',
        a3m='',
        max_template_date=datetime.date(3000, 1, 1),
        structure_store=structure_stores.StructureStore(_STRUCTURE_PATH),
    )
    self.assertEmpty(found_templates.hits)

  @parameterized.parameters(
      (datetime.date(2025, 1, 1), 4, {'abcd', 'efgh', 'ijkl'}),
      (datetime.date(2015, 1, 1), 2, {'abcd'}),
      (datetime.date(2009, 1, 1), 0, set()),
  )
  def test_filter_release_date(
      self, release_date, expected_num_filtered_hits, expected_pdb_ids
  ):
    found_templates = templates.Templates(
        query_sequence=_ARTIFICIAL_QUERY_SEQUENCE,
        query_release_date=release_date,
        hits=_ARTIFICIAL_HITS_DATA,
        max_template_date=datetime.date(3000, 1, 1),
        structure_store=structure_stores.StructureStore(_STRUCTURE_PATH),
    )
    filter_config = _TEMPLATES_V2_CONFIG.filter_config
    filtered_templates = found_templates.filter(
        max_subsequence_ratio=filter_config.max_subsequence_ratio,
        min_align_ratio=filter_config.min_align_ratio,
        min_hit_length=filter_config.min_hit_length,
        deduplicate_sequences=filter_config.deduplicate_sequences,
        max_hits=filter_config.max_hits,
    )
    self.assertEqual(filtered_templates.num_hits, expected_num_filtered_hits)

    template_ids = {h.pdb_id for h in filtered_templates.hits}
    self.assertSetEqual(template_ids, expected_pdb_ids)

  @parameterized.parameters(
      ('max_subsequence_ratio', 0.75, {'abcd_B', 'efgh_C', 'ijkl_A'}, 3),
      ('min_align_ratio', 0.75, {'abcd_A', 'abcd_B'}, 2),
      ('min_hit_length', 10, {'abcd_B'}, 1),
      ('max_hits', 2, {'abcd_A', 'abcd_B'}, 2),
      ('deduplicate_sequences', True, {'abcd_A', 'efgh_C'}, 2),
  )
  def test_filter(
      self,
      filter_name: str,
      filter_value: Any,
      expected_chains: set[str] | None,
      expected_num_filtered_hits: int,
  ):
    found_templates = templates.Templates(
        query_sequence=_ARTIFICIAL_QUERY_SEQUENCE,
        hits=_ARTIFICIAL_HITS_DATA,
        max_template_date=datetime.date(3000, 1, 1),
        structure_store=structure_stores.StructureStore(_STRUCTURE_PATH),
    )

    filter_config = dataclasses.replace(
        _TEMPLATES_V2_CONFIG.filter_config, **{filter_name: filter_value}
    )
    filtered_hits = found_templates.filter(
        max_subsequence_ratio=filter_config.max_subsequence_ratio,
        min_align_ratio=filter_config.min_align_ratio,
        min_hit_length=filter_config.min_hit_length,
        deduplicate_sequences=filter_config.deduplicate_sequences,
        max_hits=filter_config.max_hits,
    )
    observed_chains = {h.full_name for h in filtered_hits.hits}
    self.assertSetEqual(observed_chains, expected_chains)  # pyrefly: ignore[bad-argument-type]
    self.assertEqual(filtered_hits.num_hits, expected_num_filtered_hits)

  @parameterized.parameters(
      (None, '3000-01-01', '3000-01-01', 3),
      ('2015-01-01', '3000-01-01', '2014-11-02', 1),
      ('1000-01-01', '3000-01-01', '1000-01-01', 0),
  )
  def test_get_structures(
      self,
      query_release_date: str | None,
      max_template_date: str,
      cutoff_date: str,
      num_expected_strucs: int,
  ):
    max_template_date = datetime.date.fromisoformat(max_template_date)  # pyrefly: ignore[bad-assignment]
    cutoff_date = datetime.date.fromisoformat(cutoff_date)  # pyrefly: ignore[bad-assignment]

    early_hits = [h for h in _PDB_HITS_DATA if h.release_date <= cutoff_date]  # pyrefly: ignore[unsupported-operation]
    expected_pdb_ids = {h.pdb_id for h in early_hits}

    structure_store = structure_stores.StructureStore(_STRUCTURE_PATH)

    query_release_date = (
        None  # pyrefly: ignore[bad-assignment]
        if query_release_date is None
        else datetime.date.fromisoformat(query_release_date)
    )
    found_templates = templates.Templates(
        query_sequence=_ARTIFICIAL_QUERY_SEQUENCE,
        query_release_date=query_release_date,  # pyrefly: ignore[bad-argument-type]
        hits=early_hits,
        max_template_date=max_template_date,  # pyrefly: ignore[bad-argument-type]
        structure_store=structure_store,
    )
    struc_by_target_name = {
        str(struc.name).lower(): struc for struc in found_templates.structures
    }
    self.assertLen(list(found_templates.structures), num_expected_strucs)
    self.assertEqual(set(struc_by_target_name.keys()), expected_pdb_ids)

    for target_name, loaded_struc in struc_by_target_name.items():
      with self.subTest(target_name):
        mmcif = structure_store.get_mmcif_str(target_name)
        expected_struc = structure.from_mmcif(
            mmcif,
            fix_mse_residues=True,
            fix_arginines=True,
            include_water=False,
            include_other=True,
            include_bonds=False,
        )
        self.assertStructuresEqual(loaded_struc, expected_struc)

  def test_get_structures_with_invalid_cutoff_date(self):

    found_templates = templates.Templates(
        query_sequence=_ARTIFICIAL_QUERY_SEQUENCE,
        hits=_ARTIFICIAL_HITS_DATA,
        max_template_date=datetime.date(1000, 1, 1),
        structure_store=structure_stores.StructureStore(_STRUCTURE_PATH),
    )
    with self.assertRaises(templates.HitDateError):
      _ = list(found_templates.structures)

  @parameterized.parameters(
      (None, '3000-01-01', '3000-01-01'),
      ('4000-01-01', '3000-01-01', '3000-01-01'),
      ('2018-01-01', '3000-01-01', '2017-11-02'),
  )
  def test_release_date_cutoff(
      self,
      query_release_date: str | None,
      max_template_date: str,
      expected_cutoff_date: str,
  ):
    max_template_date = datetime.date.fromisoformat(max_template_date)  # pyrefly: ignore[bad-assignment]
    query_release_date = (
        None  # pyrefly: ignore[bad-assignment]
        if query_release_date is None
        else datetime.date.fromisoformat(query_release_date)
    )

    found_templates = templates.Templates(
        query_sequence=_ARTIFICIAL_QUERY_SEQUENCE,
        query_release_date=query_release_date,  # pyrefly: ignore[bad-argument-type]
        hits=_ARTIFICIAL_HITS_DATA,
        max_template_date=max_template_date,  # pyrefly: ignore[bad-argument-type]
        structure_store=structure_stores.StructureStore(_STRUCTURE_PATH),
    )

    self.assertEqual(found_templates.query_release_date, query_release_date)

    self.assertEqual(
        found_templates.release_date_cutoff,
        datetime.date.fromisoformat((expected_cutoff_date)),
    )

  @parameterized.named_parameters(
      (
          'protein',
          ''.join(list(residue_names.PROTEIN_TYPES_ONE_LETTER_WITH_UNKNOWN)),
          mmcif_names.PROTEIN_CHAIN,
          atom_types.ATOM37_ORDER,
          list(range(len(residue_names.PROTEIN_TYPES_ONE_LETTER_WITH_UNKNOWN))),
      ),
      (
          'dna',
          ''.join(residue_names.DNA_TYPES_ONE_LETTER_WITH_UNKNOWN) * 4,
          mmcif_names.DNA_CHAIN,
          atom_types.ATOM29_ORDER,
          [26, 27, 28, 29, 30] * 4,
      ),
      (
          'rna',
          ''.join(residue_names.RNA_TYPES_ONE_LETTER_WITH_UNKNOWN) * 6,
          mmcif_names.RNA_CHAIN,
          atom_types.ATOM29_ORDER,
          [22, 23, 24, 25, 30] * 6,
      ),
  )
  def test_featurize(
      self,
      query_sequence: str,
      polymer_type: str,
      atom_order: Mapping[str, int],
      expected_template_aatype: Sequence[int],
  ):
    template_sequence = query_sequence
    test_struc = _get_struc(sequence=query_sequence, polymer_type=polymer_type)

    hits = (
        templates.Hit(
            pdb_id='1avg',
            auth_chain_id='A',
            hmmsearch_sequence=template_sequence,
            structure_sequence=template_sequence,
            query_sequence=query_sequence,
            unresolved_res_indices=(),
            start_index=0,
            end_index=len(template_sequence),
            full_length=len(template_sequence),
            release_date=datetime.date(1997, 9, 16),
            chain_poly_type=polymer_type,
        ),
    )

    with mock.patch.object(
        structure, 'from_mmcif', autospec=True, return_value=test_struc
    ):
      found_templates = templates.Templates(
          hits=hits,
          query_sequence=query_sequence,
          structure_store=structure_stores.StructureStore({'1avg': ''}),
          max_template_date=datetime.date(2020, 10, 11),
          query_release_date=None,
      )
      featurized_output = found_templates.featurize(
          include_ligand_features=True
      )

    res_arrays = test_struc.filter(chain_id='A').to_res_arrays(
        include_missing_residues=True, atom_order=atom_order
    )
    positions = res_arrays.atom_positions
    positions_mask = res_arrays.atom_mask
    expected_ligand_features = {
        'B': {
            'ligand_atom_positions': np.ones((1, 3)),
            'ligand_atom_names': np.array(['']),
            'ligand_atom_occupancies': np.array([1]),
            'ccd_id': 'UNKNOWN'.encode(),
        }
    }
    expected_features = {
        'template_all_atom_positions': np.array([positions]),
        'template_all_atom_masks': np.array([positions_mask]),
        'template_sequence': np.array(template_sequence.encode()),
        'template_aatype': np.array([expected_template_aatype], dtype=np.int32),
        'template_domain_names': np.array(['1avg_A'.encode()], dtype=object),
        'template_release_date': np.array(
            ['1997-09-16'.encode()], dtype=object
        ),
        'ligand_features': [expected_ligand_features],
    }
    chex.assert_trees_all_equal(featurized_output, expected_features)

  @parameterized.named_parameters(
      (
          'same_length',
          'ACDEikFnnGH',
          'ACDEIKFXXGH',
          0,
          {0: 0, 1: 1, 2: 2, 3: 3, 4: 6, 5: 9, 6: 10},
          'ACDEFGH',
      ),
      (
          'shifted_start',
          'ACDEikFnnGH',
          'ACDEIKFXXGH',
          1,
          {0: 1, 1: 2, 2: 3, 3: 4, 4: 7, 5: 10, 6: 11},
          'ACDEFGH',
      ),
      (
          'different_length',
          'ACDEikFnnGH',
          'ACFWNGH',
          0,
          {0: 0, 1: 1, 4: 2, 5: 5, 6: 6},
          'AC--FGH',
      ),
  )
  def test_query_to_hit_mapping(
      self,
      hit_sequence: str,
      struc_sequence: str,
      start_index: int,
      expected_query_to_hit_mapping: Mapping[int, int],
      expected_aligned_hit: str,
  ):
    # Pad the seqres and structure with residues X and Y, testing that we can
    # handle hits smaller than the structure sequence in cases where the seqres
    # does not match the structure seq. Padding and alignment shown below:
    # Query :  ACDE  F  GH
    # Hit   :  ACDEIKFNNGH   (Only a subsequence of the actual seqres/struc seq)
    # Seqres: XACDEIKFNNGHYYm
    # Struc1: XACDEIKFXXGHYY
    # Struc2: XAC    FWNGHYY
    query_sequence = 'ACDEFGH'
    struc_sequence = 'X' * start_index + struc_sequence + 'Y' * 2
    hit = templates.Hit(
        pdb_id='1avg',
        auth_chain_id='A',
        hmmsearch_sequence=hit_sequence,
        structure_sequence=struc_sequence,
        unresolved_res_indices=(),
        query_sequence=query_sequence,
        start_index=start_index,
        end_index=start_index + len(hit_sequence),
        full_length=start_index + len(hit_sequence) + 2,
        release_date=datetime.date(1997, 9, 16),
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
    )

    self.assertEqual(hit.query_to_hit_mapping, expected_query_to_hit_mapping)
    aligned_hit = ['-'] * len(query_sequence)
    for query_index in range(len(query_sequence)):
      hit_index = hit.query_to_hit_mapping.get(query_index)
      if hit_index is not None:
        aligned_hit[query_index] = struc_sequence[hit_index]
    self.assertSequenceEqual(''.join(aligned_hit), expected_aligned_hit)

  def test_featurize_protein_multiple_hits(self):
    template_sequence = query_sequence = ''.join(
        list(residue_names.PROTEIN_TYPES_ONE_LETTER_WITH_UNKNOWN)
    )
    test_struc = _get_struc(
        sequence=query_sequence,
        polymer_type=mmcif_names.PROTEIN_CHAIN,
    )

    hit = templates.Hit(
        pdb_id='1avg',
        auth_chain_id='A',
        hmmsearch_sequence=template_sequence,
        structure_sequence=template_sequence,
        unresolved_res_indices=(),
        query_sequence=query_sequence,
        start_index=0,
        end_index=len(template_sequence),
        full_length=len(template_sequence),
        release_date=datetime.date(1997, 9, 16),
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
    )
    hits = (hit, hit)

    with mock.patch.object(
        structure, 'from_mmcif', autospec=True, return_value=test_struc
    ):
      found_templates = templates.Templates(
          hits=hits,
          query_sequence=query_sequence,
          max_template_date=datetime.date(2020, 10, 11),
          structure_store=structure_stores.StructureStore({'1avg': ''}),
          query_release_date=None,
      )
      featurized_output = dict(
          jax.tree.leaves_with_path(found_templates.featurize())
      )

    res_arrays = test_struc.filter(chain_id='A').to_res_arrays(
        include_missing_residues=True, atom_order=atom_types.ATOM37_ORDER
    )
    positions = res_arrays.atom_positions
    positions_mask = res_arrays.atom_mask
    expected_ligand_features = {
        'B': {
            'ligand_atom_positions': np.ones((1, 3), dtype=np.float32),
            'ligand_atom_names': np.array([''], dtype=object),
            'ligand_atom_occupancies': np.array([1], dtype=np.float32),
            'ccd_id': 'UNKNOWN'.encode(),
        }
    }
    expected_template_aatype = list(
        range(len(residue_names.PROTEIN_TYPES_ONE_LETTER_WITH_UNKNOWN))
    )
    expected_features = dict(
        jax.tree.leaves_with_path({
            'template_all_atom_positions': np.array(
                [positions, positions], dtype=np.float64
            ),
            'template_all_atom_masks': np.array(
                [positions_mask, positions_mask], dtype=np.int64
            ),
            'template_sequence': np.array(
                [template_sequence.encode(), template_sequence.encode()]
            ),
            'template_aatype': np.array(
                [expected_template_aatype, expected_template_aatype],
                dtype=np.int32,
            ),
            'template_domain_names': np.array(
                ['1avg_A'.encode(), '1avg_A'.encode()], dtype=object
            ),
            'template_release_date': np.array(
                ['1997-09-16'.encode(), '1997-09-16'.encode()], dtype=object
            ),
            'ligand_features': [
                expected_ligand_features,
                expected_ligand_features,
            ],
        })
    )
    for feature in expected_features:
      with self.subTest(feature):
        np.testing.assert_array_equal(
            featurized_output[feature], expected_features[feature], strict=True
        )

  def test_featurize_protein_not_enough_atoms(self):
    test_struc = _get_struc(
        sequence='A', polymer_type=mmcif_names.PROTEIN_CHAIN, add_ligand=False
    )
    hits = (
        templates.Hit(
            pdb_id='1avg',
            auth_chain_id='A',
            hmmsearch_sequence='A',
            structure_sequence='A',
            unresolved_res_indices=(0,),
            query_sequence=_ARTIFICIAL_QUERY_SEQUENCE,
            start_index=0,
            end_index=1,
            full_length=1,
            release_date=datetime.date(1997, 9, 16),
            chain_poly_type=mmcif_names.PROTEIN_CHAIN,
        ),
    )

    found_templates = templates.Templates(
        hits=hits,
        query_sequence=_ARTIFICIAL_QUERY_SEQUENCE,
        max_template_date=datetime.date(2020, 10, 11),
        structure_store=structure_stores.StructureStore({'1avg': ''}),
        query_release_date=None,
    )
    with mock.patch.object(
        structure, 'from_mmcif', autospec=True, return_value=test_struc
    ):
      with (
          self.subTest('FeaturiseWithoutFilter'),
          self.assertRaises(templates.InvalidTemplateError),
      ):
        found_templates.featurize()
      filter_config = _TEMPLATES_V2_CONFIG.filter_config
      features = found_templates.filter(
          max_subsequence_ratio=filter_config.max_subsequence_ratio,
          min_align_ratio=filter_config.min_align_ratio,
          min_hit_length=filter_config.min_hit_length,
          deduplicate_sequences=filter_config.deduplicate_sequences,
          max_hits=filter_config.max_hits,
      ).featurize()
      chex.assert_trees_all_equal(
          features,
          {
              'template_release_date': np.array([], dtype=object),
              'template_all_atom_masks': np.array([], dtype=np.float64),
              'template_sequence': np.array([], dtype=object),
              'template_domain_names': np.array([], dtype=object),
              'ligand_features': [],
              'template_aatype': np.array([], dtype=np.int32),
              'template_all_atom_positions': np.array([], dtype=np.float64),
          },
      )

  def test_get_hits_with_structures_protein_multiple_hits(self):
    template_sequence = query_sequence = ''.join(
        list(residue_names.PROTEIN_TYPES_ONE_LETTER_WITH_UNKNOWN)
    )
    test_struc = _get_struc(
        sequence=query_sequence,
        polymer_type=mmcif_names.PROTEIN_CHAIN,
    )

    hit = templates.Hit(
        pdb_id='1avg',
        auth_chain_id='A',
        hmmsearch_sequence=template_sequence,
        structure_sequence=template_sequence,
        unresolved_res_indices=(),
        query_sequence=query_sequence,
        start_index=0,
        end_index=len(template_sequence),
        full_length=len(template_sequence),
        release_date=datetime.date(1997, 9, 16),
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
    )
    hits = (hit, hit)

    with mock.patch.object(
        structure, 'from_mmcif', autospec=True, return_value=test_struc
    ):
      found_templates = templates.Templates(
          hits=hits,
          query_sequence=query_sequence,
          max_template_date=datetime.date(2020, 10, 11),
          structure_store=structure_stores.StructureStore({'1avg': ''}),
          query_release_date=None,
      )
      hits_and_strucs = found_templates.get_hits_with_structures()
    self.assertLen(hits_and_strucs, 2)


def _get_struc(
    *,
    sequence: str,
    polymer_type: str,
    add_ligand: bool = True,
) -> structure.Structure:
  sequences = [sequence]
  chain_types = [polymer_type]
  sequence_formats = [structure.SequenceFormat.FASTA]

  if add_ligand:
    sequences.append('UNKNOWN:CC')
    chain_types.append(mmcif_names.NON_POLYMER_CHAIN)
    sequence_formats.append(structure.SequenceFormat.LIGAND_SMILES)

  struc = structure.from_sequences_and_bonds(
      sequences=sequences,
      chain_types=chain_types,
      sequence_formats=sequence_formats,
      bonded_atom_pairs=None,
      ccd=chemical_components.Ccd(),
      name='1avg',
      release_date=datetime.date(1997, 9, 16),
  )
  # Make sure the ligand has the same author chain ID as the polymer chain.
  new_chains = struc.chains_table.copy_and_update(
      auth_asym_id=np.full(struc.num_chains, fill_value='A', dtype=object),
  )
  struc = struc.copy_and_update(chains=new_chains)
  updated_coords = np.ones_like(struc.coords) * np.transpose([struc.res_id])
  struc = struc.copy_and_update_coords(coords=updated_coords)
  return struc


if __name__ == '__main__':
  absltest.main()
