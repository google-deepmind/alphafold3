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

"""Inexpensive AlphaFold 3 input checks that do not run a search or inference.

Protein MSA/template rules match `DataPipeline.process_protein_chain`. This
module only classifies or rejects that combination. It does not call
Jackhmmer, Nhmmer, or the model.
"""

import enum
from typing import Any


class ProteinInputStatus(enum.Enum):
  """How a protein chain's MSA and template fields are filled in."""

  # Both MSAs unset and templates unset: search for MSA and templates.
  SEARCH_MSA_AND_TEMPLATES = enum.auto()
  # Both MSAs unset and templates explicitly empty: search for MSA only.
  SEARCH_MSA = enum.auto()
  # Both MSAs set and templates unset: search for templates only.
  SEARCH_TEMPLATES = enum.auto()
  # Both MSAs set and templates set, including empty string / empty list.
  SUPPLIED = enum.auto()


def protein_input_status(
    chain_id: str,
    unpaired_msa: str | None,
    paired_msa: str | None,
    templates: Any,
) -> ProteinInputStatus:
  """Returns how protein MSA/template fields should be handled.

  Args:
    chain_id: Chain id, included in the existing partial-input error.
    unpaired_msa: Unpaired A3M, or None when the data pipeline should build it.
    paired_msa: Paired A3M, or None when the data pipeline should build it.
    templates: Template list, an empty list when templates are disabled, or
      None when the data pipeline should search for templates.

  Raises:
    ValueError: If unpaired MSA, paired MSA, and templates are only partly set.
      The message is the one raised by the data pipeline.
  """
  has_unpaired_msa = unpaired_msa is not None
  has_paired_msa = paired_msa is not None
  has_templates = templates is not None

  if not has_unpaired_msa and not has_paired_msa and not templates:
    if templates is None:
      return ProteinInputStatus.SEARCH_MSA_AND_TEMPLATES
    return ProteinInputStatus.SEARCH_MSA
  if has_unpaired_msa and has_paired_msa and not has_templates:
    return ProteinInputStatus.SEARCH_TEMPLATES
  if has_unpaired_msa and has_paired_msa and has_templates:
    return ProteinInputStatus.SUPPLIED
  raise ValueError(
      f'Protein chain {chain_id} has unpaired MSA, paired MSA, or'
      ' templates set only partially. If you want to run the pipeline'
      ' with custom MSA/templates, you need to set all of them. You can'
      ' set MSA to empty string and templates to empty list to signify'
      ' that they should not be used and searched for.'
  )


def describe_protein_input(
    chain_id: str,
    unpaired_msa: str | None,
    paired_msa: str | None,
    templates: Any,
) -> str:
  """Describes whether a protein chain still needs the data pipeline."""
  status = protein_input_status(
      chain_id, unpaired_msa, paired_msa, templates
  )
  if status is ProteinInputStatus.SEARCH_MSA_AND_TEMPLATES:
    return (
        f'protein {chain_id}: MSA and templates require the data pipeline'
    )
  if status is ProteinInputStatus.SEARCH_MSA:
    return (
        f'protein {chain_id}: MSA requires the data pipeline; templates'
        ' disabled'
    )
  if status is ProteinInputStatus.SEARCH_TEMPLATES:
    return (
        f'protein {chain_id}: MSA supplied; templates require the data'
        ' pipeline'
    )
  if not unpaired_msa and not paired_msa:
    msa_note = 'empty'
  else:
    msa_note = 'custom'
  if not templates:
    template_note = 'disabled'
  else:
    template_note = f'{len(templates)} supplied'
  return f'protein {chain_id}: MSA {msa_note}; templates {template_note}'


def describe_rna_msa(chain_id: str, unpaired_msa: str | None) -> str:
  """Describes whether an RNA chain still needs an MSA search."""
  if unpaired_msa is None:
    return f'rna {chain_id}: MSA requires the data pipeline'
  if not unpaired_msa:
    return f'rna {chain_id}: MSA supplied empty'
  return f'rna {chain_id}: MSA supplied'


def preflight_fold_input(fold_input, source: str) -> str:
  """Validates one parsed fold input and returns a text summary.

  JSON schema, chain, bond, SMILES, path, and user-CCD key checks have already
  run while parsing `fold_input`. This adds the pipeline's MSA/template field
  rule, parses supplied A3M with the existing MSA parser, checks supplied
  template mmCIFs with the same parser featurisation uses, and resolves
  explicit CCD references against the local dictionary when that file exists.

  Args:
    fold_input: A parsed `folding_input.Input`.
    source: Input file path to print in the summary.

  Returns:
    A multi-line summary ending in `Result: valid`.

  Raises:
    ValueError: The existing validation error for this input.
  """
  if not fold_input.chains:
    raise ValueError('Fold input has no chains.')

  chain_lines = []
  for chain in fold_input.protein_chains:
    chain_lines.append(
        describe_protein_input(
            chain.id, chain.unpaired_msa, chain.paired_msa, chain.templates
        )
    )
  for chain in fold_input.rna_chains:
    chain_lines.append(describe_rna_msa(chain.id, chain.unpaired_msa))
  for chain in fold_input.ligands:
    if chain.ccd_ids:
      chain_lines.append(
          f'ligand {chain.id}: CCD {", ".join(chain.ccd_ids)}'
      )
    else:
      chain_lines.append(f'ligand {chain.id}: SMILES')

  _validate_supplied_msas(fold_input)
  _validate_supplied_templates(fold_input)
  ccd_line = _validate_ccd_references(fold_input)
  seeds = ', '.join(str(seed) for seed in fold_input.rng_seeds)
  lines = [
      f'Input: {source}',
      f'Job: {fold_input.name}',
      (
          'Chains:'
          f' protein={len(fold_input.protein_chains)},'
          f' rna={len(fold_input.rna_chains)},'
          f' dna={len(fold_input.dna_chains)},'
          f' ligand={len(fold_input.ligands)}'
      ),
      f'Seeds: {seeds}',
      f'Custom CCD: {"yes" if fold_input.user_ccd else "no"}',
      *chain_lines,
      ccd_line,
      'Result: valid',
  ]
  return '\n'.join(lines)


def _validate_supplied_msas(fold_input) -> None:
  """Parses user A3Ms that featurisation will actually read.

  Featurisation skips MSA parsing for chains of 4 or fewer residues. Paired
  protein MSAs are parsed only when more than one distinct protein sequence is
  present. Empty strings are valid MSA-free inputs and are left untouched.
  """
  from alphafold3.constants import mmcif_names
  from alphafold3.data import msa

  protein_sequences = [chain.sequence for chain in fold_input.protein_chains]
  parse_paired = len(set(protein_sequences)) > 1

  def _check(chain_id, kind, sequence, poly_type, a3m):
    try:
      msa.Msa.from_a3m(
          query_sequence=sequence,
          chain_poly_type=poly_type,
          a3m=a3m,
          deduplicate=False,
      )
    except ValueError as e:
      raise ValueError(
          f'{kind} MSA for chain {chain_id} failed validation: {e}'
      ) from e

  for chain in fold_input.protein_chains:
    if len(chain) <= 4:
      continue
    if chain.unpaired_msa:
      _check(
          chain.id,
          'unpaired',
          chain.sequence,
          mmcif_names.PROTEIN_CHAIN,
          chain.unpaired_msa,
      )
    if parse_paired and chain.paired_msa:
      _check(
          chain.id,
          'paired',
          chain.sequence,
          mmcif_names.PROTEIN_CHAIN,
          chain.paired_msa,
      )
  for chain in fold_input.rna_chains:
    if len(chain) <= 4 or not chain.unpaired_msa:
      continue
    _check(
        chain.id,
        'unpaired',
        chain.sequence,
        mmcif_names.RNA_CHAIN,
        chain.unpaired_msa,
    )


def _validate_supplied_templates(fold_input) -> None:
  """Checks supplied template mmCIFs with the featurisation template parser.

  The same `from_mmcif` arguments and `get_polymer_features` call are used.
  Chains of 4 or fewer residues, and later copies of a protein sequence, do
  not reach that call during featurisation, so they are not checked here.
  """
  from alphafold3 import structure
  from alphafold3.constants import mmcif_names
  from alphafold3.data import templates

  seen_sequences = set()
  for chain in fold_input.protein_chains:
    sequence_key = tuple(chain.to_ccd_sequence())
    first_copy = sequence_key not in seen_sequences
    seen_sequences.add(sequence_key)
    if not first_copy or len(chain) <= 4 or not chain.templates:
      continue
    for template_index, template in enumerate(chain.templates, start=1):
      try:
        template_structure = structure.from_mmcif(
            template.mmcif,
            fix_mse_residues=True,
            fix_arginines=True,
            include_bonds=False,
            include_water=False,
            include_other=True,
        )
        templates.get_polymer_features(
            chain=template_structure,
            chain_poly_type=mmcif_names.PROTEIN_CHAIN,
            query_sequence_length=len(chain),
            query_to_hit_mapping=template.query_to_template_map,
        )
      except (ValueError, IndexError, KeyError) as e:
        raise ValueError(
            f'Template {template_index} for protein chain {chain.id} failed'
            f' validation: {e}'
        ) from e


def _explicit_ccd_codes(fold_input) -> list[str]:
  """CCD ids named by the user, in input order, without standard residues."""
  from alphafold3.common import folding_input

  codes = []
  for chain in fold_input.chains:
    if isinstance(chain, folding_input.ProteinChain):
      codes.extend(code for code, _ in chain.ptms)
    elif isinstance(chain, folding_input.RnaChain):
      codes.extend(code for code, _ in chain.modifications)
    elif isinstance(chain, folding_input.DnaChain):
      codes.extend(code for code, _ in chain.modifications)
    elif isinstance(chain, folding_input.Ligand) and chain.ccd_ids:
      codes.extend(chain.ccd_ids)
  return codes


def _user_ccd_component_ids(user_ccd: str | None) -> set[str]:
  if not user_ccd:
    return set()
  from alphafold3.cpp import cif_dict

  return set(cif_dict.parse_multi_data_cif(user_ccd))


def _validate_ccd_references(fold_input) -> str:
  """Resolves user-named CCD codes locally, including bond atom names."""
  codes = _explicit_ccd_codes(fold_input)
  unique_codes = list(dict.fromkeys(codes))
  if not unique_codes and not fold_input.bonded_atom_pairs:
    return 'CCD references: none'

  from alphafold3.constants import chemical_components

  try:
    ccd = chemical_components.Ccd(user_ccd=fold_input.user_ccd)
  except FileNotFoundError as e:
    user_ids = _user_ccd_component_ids(fold_input.user_ccd)
    missing = [code for code in unique_codes if code not in user_ids]
    if missing or fold_input.bonded_atom_pairs:
      unresolved = ', '.join(missing) if missing else 'bonded atom names'
      raise ValueError(
          'Cannot resolve local chemical component references'
          f' ({unresolved}). The CCD dictionary is not available. {e}'
      ) from e
    listed = ', '.join(unique_codes) if unique_codes else 'user CCD only'
    return f'CCD references: defined in user CCD ({listed})'

  # Same structure conversion inference uses later. This rejects unknown CCD
  # codes and bond atom names that are not on that component.
  fold_input.to_structure(ccd)
  if unique_codes:
    return 'CCD references: resolved locally (' + ', '.join(unique_codes) + ')'
  return 'CCD references: resolved locally'
