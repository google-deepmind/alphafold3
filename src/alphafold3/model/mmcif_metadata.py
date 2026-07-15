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

"""Adds mmCIF metadata (to be ModelCIF-conformant) and author and legal info."""

from typing import Final

from alphafold3.structure import mmcif
import numpy as np

_LICENSE_URL: Final[str] = (
    'https://github.com/google-deepmind/alphafold3/blob/main/OUTPUT_TERMS_OF_USE.md'
)

_LICENSE: Final[str] = f"""
Non-commercial use only, by using this file you agree to the terms of use found
at {_LICENSE_URL}.
To request access to the AlphaFold 3 model parameters, follow the process set
out at https://github.com/google-deepmind/alphafold3. You may only use these if
received directly from Google. Use is subject to terms of use available at
https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md.
""".strip()

_DISCLAIMER: Final[str] = """\
AlphaFold 3 and its output are not intended for, have not been validated for,
and are not approved for clinical use. They are provided "as-is" without any
warranty of any kind, whether expressed or implied. No warranty is given that
use shall not infringe the rights of any third party.
""".strip()

_MMCIF_PAPER_AUTHORS: Final[tuple[str, ...]] = (
    'Google DeepMind',
    'Isomorphic Labs',
)

# Authors of the mmCIF - we set them to be equal to the authors of the paper.
_MMCIF_AUTHORS: Final[tuple[str, ...]] = _MMCIF_PAPER_AUTHORS


def add_metadata_to_mmcif(
    old_cif: mmcif.Mmcif, version: str, model_id: bytes
) -> mmcif.Mmcif:
  """Adds metadata to a mmCIF to make it ModelCIF-conformant."""
  cif = {}

  # ModelCIF conformation dictionary.
  cif['_audit_conform.dict_name'] = ['mmcif_ma.dic']
  cif['_audit_conform.dict_version'] = ['1.4.5']
  cif['_audit_conform.dict_location'] = [
      'https://raw.githubusercontent.com/ihmwg/ModelCIF/master/dist/mmcif_ma.dic'
  ]

  cif['_pdbx_data_usage.id'] = ['1', '2']
  cif['_pdbx_data_usage.type'] = ['license', 'disclaimer']
  cif['_pdbx_data_usage.details'] = [_LICENSE, _DISCLAIMER]
  cif['_pdbx_data_usage.url'] = [_LICENSE_URL, '?']

  # Structure author details.
  cif['_audit_author.name'] = []
  cif['_audit_author.pdbx_ordinal'] = []
  for author_index, author_name in enumerate(_MMCIF_AUTHORS, start=1):
    cif['_audit_author.name'].append(author_name)
    cif['_audit_author.pdbx_ordinal'].append(str(author_index))

  # Paper author details.
  cif['_citation_author.citation_id'] = []
  cif['_citation_author.name'] = []
  cif['_citation_author.ordinal'] = []
  for author_index, author_name in enumerate(_MMCIF_PAPER_AUTHORS, start=1):
    cif['_citation_author.citation_id'].append('primary')
    cif['_citation_author.name'].append(author_name)
    cif['_citation_author.ordinal'].append(str(author_index))

  # Paper citation details.
  cif['_citation.id'] = ['primary']
  cif['_citation.title'] = [
      'Accurate structure prediction of biomolecular interactions with'
      ' AlphaFold 3'
  ]
  cif['_citation.journal_full'] = ['Nature']
  cif['_citation.journal_volume'] = ['630']
  cif['_citation.page_first'] = ['493']
  cif['_citation.page_last'] = ['500']
  cif['_citation.year'] = ['2024']
  cif['_citation.journal_id_ASTM'] = ['NATUAS']
  cif['_citation.country'] = ['UK']
  cif['_citation.journal_id_ISSN'] = ['0028-0836']
  cif['_citation.journal_id_CSD'] = ['0006']
  cif['_citation.book_publisher'] = ['?']
  cif['_citation.pdbx_database_id_PubMed'] = ['38718835']
  cif['_citation.pdbx_database_id_DOI'] = ['10.1038/s41586-024-07487-w']

  # Type of data in the dataset including data used in the model generation.
  cif['_ma_data.id'] = ['1']
  cif['_ma_data.name'] = ['Model']
  cif['_ma_data.content_type'] = ['model coordinates']

  # Description of number of instances for each entity.
  cif['_ma_target_entity_instance.asym_id'] = old_cif['_struct_asym.id']
  cif['_ma_target_entity_instance.entity_id'] = old_cif[
      '_struct_asym.entity_id'
  ]
  cif['_ma_target_entity_instance.details'] = ['.'] * len(
      cif['_ma_target_entity_instance.entity_id']
  )

  # Details about the target entities.
  cif['_ma_target_entity.entity_id'] = cif[
      '_ma_target_entity_instance.entity_id'
  ]
  cif['_ma_target_entity.data_id'] = ['1'] * len(
      cif['_ma_target_entity.entity_id']
  )
  cif['_ma_target_entity.origin'] = ['.'] * len(
      cif['_ma_target_entity.entity_id']
  )

  # Details of the models being deposited.
  cif['_ma_model_list.ordinal_id'] = ['1']
  cif['_ma_model_list.model_id'] = ['1']
  cif['_ma_model_list.model_group_id'] = ['1']
  cif['_ma_model_list.model_name'] = ['Top ranked model']

  cif['_ma_model_list.model_group_name'] = [
      f'AlphaFold-beta-20231127 ({version})'
  ]
  cif['_ma_model_list.data_id'] = ['1']
  cif['_ma_model_list.model_type'] = ['Ab initio model']

  # Software used.
  cif['_software.pdbx_ordinal'] = ['1']
  cif['_software.name'] = ['AlphaFold']
  cif['_software.version'] = [
      f'AlphaFold-beta-20231127 ({model_id.decode("ascii")})'
  ]
  cif['_software.type'] = ['package']
  cif['_software.description'] = ['Structure prediction']
  cif['_software.classification'] = ['other']
  cif['_software.date'] = ['?']

  # Collection of software into groups.
  cif['_ma_software_group.ordinal_id'] = ['1']
  cif['_ma_software_group.group_id'] = ['1']
  cif['_ma_software_group.software_id'] = ['1']

  # Method description to conform with ModelCIF.
  cif['_ma_protocol_step.ordinal_id'] = ['1', '2', '3']
  cif['_ma_protocol_step.protocol_id'] = ['1', '1', '1']
  cif['_ma_protocol_step.step_id'] = ['1', '2', '3']
  cif['_ma_protocol_step.method_type'] = [
      'coevolution MSA',
      'template search',
      'modeling',
  ]

  # Details of the metrics use to assess model confidence.
  cif['_ma_qa_metric.id'] = ['1', '2']
  cif['_ma_qa_metric.name'] = ['pLDDT', 'pLDDT']
  # Accepted values are distance, energy, normalised score, other, zscore.
  cif['_ma_qa_metric.type'] = ['pLDDT', 'pLDDT']
  cif['_ma_qa_metric.mode'] = ['global', 'local']
  cif['_ma_qa_metric.software_group_id'] = ['1', '1']

  # Global model confidence pLDDT value.
  cif['_ma_qa_metric_global.ordinal_id'] = ['1']
  cif['_ma_qa_metric_global.model_id'] = ['1']
  cif['_ma_qa_metric_global.metric_id'] = ['1']
  # Mean over all atoms, since AlphaFold 3 outputs pLDDT per-atom.
  global_plddt = np.mean(
      [float(v) for v in old_cif['_atom_site.B_iso_or_equiv']]
  )
  cif['_ma_qa_metric_global.metric_value'] = [f'{global_plddt:.2f}']

  # Local (per residue) model confidence pLDDT value.
  cif['_ma_qa_metric_local.ordinal_id'] = []
  cif['_ma_qa_metric_local.model_id'] = []
  cif['_ma_qa_metric_local.label_asym_id'] = []
  cif['_ma_qa_metric_local.label_seq_id'] = []
  cif['_ma_qa_metric_local.label_comp_id'] = []
  cif['_ma_qa_metric_local.metric_id'] = []
  cif['_ma_qa_metric_local.metric_value'] = []

  plddt_grouped_by_res = {}
  for *res, atom_plddt in zip(
      old_cif['_atom_site.label_asym_id'],
      old_cif['_atom_site.label_seq_id'],
      old_cif['_atom_site.label_comp_id'],
      old_cif['_atom_site.B_iso_or_equiv'],
  ):
    plddt_grouped_by_res.setdefault(tuple(res), []).append(float(atom_plddt))

  for ordinal_id, ((chain_id, res_id, res_name), res_plddts) in enumerate(
      plddt_grouped_by_res.items(), start=1
  ):
    res_plddt = np.mean(res_plddts)
    cif['_ma_qa_metric_local.ordinal_id'].append(str(ordinal_id))
    cif['_ma_qa_metric_local.model_id'].append('1')
    cif['_ma_qa_metric_local.label_asym_id'].append(chain_id)
    cif['_ma_qa_metric_local.label_seq_id'].append(res_id)
    cif['_ma_qa_metric_local.label_comp_id'].append(res_name)
    cif['_ma_qa_metric_local.metric_id'].append('2')  # See _ma_qa_metric.id.
    cif['_ma_qa_metric_local.metric_value'].append(f'{res_plddt:.2f}')

  cif['_atom_type.symbol'] = sorted(set(old_cif['_atom_site.type_symbol']))

  return old_cif.copy_and_update(cif)


def add_legal_comment(cif: str) -> str:
  """Adds legal comment at the top of the mmCIF."""
  # fmt: off
  # pylint: disable=line-too-long
  comment = (
      '# By using this file you agree to the legally binding terms of use found at\n'
      f'# {_LICENSE_URL}.\n'
      '# To request access to the AlphaFold 3 model parameters, follow the process set\n'
      '# out at https://github.com/google-deepmind/alphafold3. You may only use these if\n'
      '# received directly from Google. Use is subject to terms of use available at\n'
      '# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md.'
  )
  # pylint: enable=line-too-long
  # fmt: on
  return f'{comment}\n{cif}'
