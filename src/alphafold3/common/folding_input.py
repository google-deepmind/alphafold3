# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Model input dataclass."""

from collections.abc import Collection, Mapping, Sequence
import dataclasses
import json
import logging
import pathlib
import random
import re
import string
from typing import Any, Final, Self, TypeAlias

from alphafold3 import structure
from alphafold3.constants import chemical_components
from alphafold3.constants import mmcif_names
from alphafold3.constants import residue_names
from alphafold3.structure import mmcif as mmcif_lib
import rdkit.Chem as rd_chem


BondAtomId: TypeAlias = tuple[str, int, str]

JSON_DIALECT: Final[str] = 'alphafold3'
JSON_VERSION: Final[int] = 1

ALPHAFOLDSERVER_JSON_DIALECT: Final[str] = 'alphafoldserver'
ALPHAFOLDSERVER_JSON_VERSION: Final[int] = 1


def _validate_keys(actual: Collection[str], expected: Collection[str]):
  """Validates that the JSON doesn't contain any extra unwanted keys."""
  if bad_keys := set(actual) - set(expected):
    raise ValueError(f'Unexpected JSON keys in: {", ".join(sorted(bad_keys))}')


class Template:
  """Structural template input."""

  __slots__ = ('_mmcif', '_query_to_template')

  def __init__(self, mmcif: str, query_to_template_map: Mapping[int, int]):
    """Initializes the template.

    Args:
      mmcif: The structural template in mmCIF format. The mmCIF should have only
        one protein chain.
      query_to_template_map: A mapping from query residue index to template
        residue index.
    """
    self._mmcif = mmcif
    # Needed to make the Template class hashable.
    self._query_to_template = tuple(query_to_template_map.items())

  @property
  def query_to_template_map(self) -> Mapping[int, int]:
    return dict(self._query_to_template)

  @property
  def mmcif(self) -> str:
    return self._mmcif

  def __hash__(self) -> int:
    return hash((self._mmcif, tuple(sorted(self._query_to_template))))

  def __eq__(self, other: Self) -> bool:
    mmcifs_equal = self._mmcif == other._mmcif
    maps_equal = sorted(self._query_to_template) == sorted(
        other._query_to_template
    )
    return mmcifs_equal and maps_equal


@dataclasses.dataclass(frozen=True, slots=True, kw_only=True)
class ProteinChain:
  """Protein chain input.

  Attributes:
    id: Unique protein chain identifier.
    sequence: The amino acid sequence of the chain.
    ptms: A list of tuples containing the post-translational modification type
      and the (1-based) residue index where the modification is applied.
    paired_msa: Paired A3M-formatted MSA for this chain. This MSA is not
      deduplicated and will be used to compute paired features. If None, this
      field is unset and must be filled in by the data pipeline before
      featurisation. If set to an empty string, it will be treated as a custom
      MSA with no sequences.
    unpaired_msa: Unpaired A3M-formatted MSA for this chain. This will be
      deduplicated and used to compute unpaired features. If None, this field is
      unset and must be filled in by the data pipeline before featurisation. If
      set to an empty string, it will be treated as a custom MSA with no
      sequences.
    templates: A list of structural templates for this chain. If None, this
      field is unset and must be filled in by the data pipeline before
      featurisation. The list can be empty or contain up to 20 templates.
  """

  id: str
  sequence: str
  ptms: Sequence[tuple[str, int]]
  paired_msa: str | None = None
  unpaired_msa: str | None = None
  templates: Sequence[Template] | None = None

  def __post_init__(self):
    if not all(res.isalpha() for res in self.sequence):
      raise ValueError(
          f'Protein must contain only digits, got "{self.sequence}"'
      )
    if any(not 0 < mod[1] <= len(self.sequence) for mod in self.ptms):
      raise ValueError(f'Invalid protein modification index: {self.ptms}')

    # Use hashable types for ptms and templates.
    if self.ptms is not None:
      object.__setattr__(self, 'ptms', tuple(self.ptms))
    if self.templates is not None:
      object.__setattr__(self, 'templates', tuple(self.templates))

  @classmethod
  def from_alphafoldserver_dict(
      cls, json_dict: Mapping[str, Any], seq_id: str
  ) -> Self:
    """Constructs ProteinChain from the AlphaFoldServer JSON dict."""
    _validate_keys(
        json_dict.keys(),
        {'sequence', 'glycans', 'modifications', 'count'},
    )
    sequence = json_dict['sequence']

    if 'glycans' in json_dict:
      raise ValueError(
          f'Specifying glycans in the `{ALPHAFOLDSERVER_JSON_DIALECT}` format'
          ' is not currently supported.'
      )

    ptms = [
        (mod['ptmType'].removeprefix('CCD_'), mod['ptmPosition'])
        for mod in json_dict.get('modifications', [])
    ]
    return cls(id=seq_id, sequence=sequence, ptms=ptms)

  @classmethod
  def from_dict(
      cls, json_dict: Mapping[str, Any], seq_id: str | None = None
  ) -> Self:
    """Constructs ProteinChain from the AlphaFold JSON dict."""
    json_dict = json_dict['protein']
    _validate_keys(
        json_dict.keys(),
        {
            'id',
            'sequence',
            'modifications',
            'unpairedMsa',
            'pairedMsa',
            'templates',
        },
    )

    sequence = json_dict['sequence']
    ptms = [
        (mod['ptmType'], mod['ptmPosition'])
        for mod in json_dict.get('modifications', [])
    ]

    unpaired_msa = json_dict.get('unpairedMsa', None)
    paired_msa = json_dict.get('pairedMsa', None)

    raw_templates = json_dict.get('templates', None)

    if raw_templates is None:
      templates = None
    else:
      templates = [
          Template(
              mmcif=template['mmcif'],
              query_to_template_map=dict(
                  zip(template['queryIndices'], template['templateIndices'])
              ),
          )
          for template in raw_templates
      ]

    return cls(
        id=seq_id or json_dict['id'],
        sequence=sequence,
        ptms=ptms,
        paired_msa=paired_msa,
        unpaired_msa=unpaired_msa,
        templates=templates,
    )

  def to_dict(self) -> Mapping[str, Mapping[str, Any]]:
    """Converts ProteinChain to an AlphaFold JSON dict."""
    if self.templates is None:
      templates = None
    else:
      templates = [
          {
              'mmcif': template.mmcif,
              'queryIndices': list(template.query_to_template_map.keys()),
              'templateIndices': (
                  list(template.query_to_template_map.values()) or None
              ),
          }
          for template in self.templates
      ]
    contents = {
        'id': self.id,
        'sequence': self.sequence,
        'modifications': [
            {'ptmType': ptm[0], 'ptmPosition': ptm[1]} for ptm in self.ptms
        ],
        'unpairedMsa': self.unpaired_msa,
        'pairedMsa': self.paired_msa,
        'templates': templates,
    }
    return {'protein': contents}

  def to_ccd_sequence(self) -> Sequence[str]:
    """Converts to a sequence of CCD codes."""
    ccd_coded_seq = [
        residue_names.PROTEIN_COMMON_ONE_TO_THREE.get(res, residue_names.UNK)
        for res in self.sequence
    ]
    for ptm_code, ptm_index in self.ptms:
      ccd_coded_seq[ptm_index - 1] = ptm_code
    return ccd_coded_seq

  def fill_missing_fields(self) -> Self:
    """Fill missing MSA and template fields with default values."""
    return dataclasses.replace(
        self,
        unpaired_msa=self.unpaired_msa or '',
        paired_msa=self.paired_msa or '',
        templates=self.templates or [],
    )


@dataclasses.dataclass(frozen=True, slots=True, kw_only=True)
class RnaChain:
  """RNA chain input.

  Attributes:
    id: Unique RNA chain identifier.
    sequence: The RNA sequence of the chain.
    modifications: A list of tuples containing the modification type and the
      (1-based) residue index where the modification is applied.
    unpaired_msa: Unpaired A3M-formatted MSA for this chain. This will be
      deduplicated and used to compute unpaired features. If None, this field is
      unset and must be filled in by the data pipeline before featurisation. If
      set to an empty string, it will be treated as a custom MSA with no
      sequences.
  """

  id: str
  sequence: str
  modifications: Sequence[tuple[str, int]]
  unpaired_msa: str | None = None

  def __post_init__(self):
    if not all(res.isalpha() for res in self.sequence):
      raise ValueError(f'RNA must contain only digits, got "{self.sequence}"')
    if any(not 0 < mod[1] <= len(self.sequence) for mod in self.modifications):
      raise ValueError(f'Invalid RNA modification index: {self.modifications}')

    # Use hashable types for modifications.
    object.__setattr__(self, 'modifications', tuple(self.modifications))

  @classmethod
  def from_alphafoldserver_dict(
      cls, json_dict: Mapping[str, Any], seq_id: str
  ) -> Self:
    """Constructs RnaChain from the AlphaFoldServer JSON dict."""
    _validate_keys(json_dict.keys(), {'sequence', 'modifications', 'count'})
    sequence = json_dict['sequence']
    modifications = [
        (mod['modificationType'].removeprefix('CCD_'), mod['basePosition'])
        for mod in json_dict.get('modifications', [])
    ]
    return cls(id=seq_id, sequence=sequence, modifications=modifications)

  @classmethod
  def from_dict(
      cls, json_dict: Mapping[str, Any], seq_id: str | None = None
  ) -> Self:
    """Constructs RnaChain from the AlphaFold JSON dict."""
    json_dict = json_dict['rna']
    _validate_keys(
        json_dict.keys(), {'id', 'sequence', 'unpairedMsa', 'modifications'}
    )
    sequence = json_dict['sequence']
    modifications = [
        (mod['modificationType'], mod['basePosition'])
        for mod in json_dict.get('modifications', [])
    ]
    unpaired_msa = json_dict.get('unpairedMsa', None)
    return cls(
        id=seq_id or json_dict['id'],
        sequence=sequence,
        modifications=modifications,
        unpaired_msa=unpaired_msa,
    )

  def to_dict(self) -> Mapping[str, Mapping[str, Any]]:
    """Converts RnaChain to an AlphaFold JSON dict."""
    contents = {
        'id': self.id,
        'sequence': self.sequence,
        'modifications': [
            {'modificationType': mod[0], 'basePosition': mod[1]}
            for mod in self.modifications
        ],
        'unpairedMsa': self.unpaired_msa,
    }
    return {'rna': contents}

  def to_ccd_sequence(self) -> Sequence[str]:
    """Converts to a sequence of CCD codes."""
    mapping = {r: r for r in residue_names.RNA_TYPES}  # Same 1-letter and CCD.
    ccd_coded_seq = [
        mapping.get(res, residue_names.UNK_RNA) for res in self.sequence
    ]
    for ccd_code, modification_index in self.modifications:
      ccd_coded_seq[modification_index - 1] = ccd_code
    return ccd_coded_seq

  def fill_missing_fields(self) -> Self:
    """Fill missing MSA fields with default values."""
    return dataclasses.replace(self, unpaired_msa=self.unpaired_msa or '')


@dataclasses.dataclass(frozen=True, slots=True, kw_only=True)
class DnaChain:
  """Single strand DNA chain input.

  Attributes:
    id: Unique DNA chain identifier.
    sequence: The DNA sequence of the chain.
    modifications: A list of tuples containing the modification type and the
      (1-based) residue index where the modification is applied.
  """

  id: str
  sequence: str
  modifications: Sequence[tuple[str, int]]

  def __post_init__(self):
    if not all(res.isalpha() for res in self.sequence):
      raise ValueError(f'DNA must contain only digits, got "{self.sequence}"')
    if any(not 0 < mod[1] <= len(self.sequence) for mod in self.modifications):
      raise ValueError(f'Invalid DNA modification index: {self.modifications}')

    # Use hashable types for modifications.
    object.__setattr__(self, 'modifications', tuple(self.modifications))

  @classmethod
  def from_alphafoldserver_dict(
      cls, json_dict: Mapping[str, Any], seq_id: str
  ) -> Self:
    """Constructs DnaChain from the AlphaFoldServer JSON dict."""
    _validate_keys(json_dict.keys(), {'sequence', 'modifications', 'count'})
    sequence = json_dict['sequence']
    modifications = [
        (mod['modificationType'].removeprefix('CCD_'), mod['basePosition'])
        for mod in json_dict.get('modifications', [])
    ]
    return cls(id=seq_id, sequence=sequence, modifications=modifications)

  @classmethod
  def from_dict(
      cls, json_dict: Mapping[str, Any], seq_id: str | None = None
  ) -> Self:
    """Constructs DnaChain from the AlphaFold JSON dict."""
    json_dict = json_dict['dna']
    _validate_keys(json_dict.keys(), {'id', 'sequence', 'modifications'})
    sequence = json_dict['sequence']
    modifications = [
        (mod['modificationType'], mod['basePosition'])
        for mod in json_dict.get('modifications', [])
    ]
    return cls(
        id=seq_id or json_dict['id'],
        sequence=sequence,
        modifications=modifications,
    )

  def to_dict(self) -> Mapping[str, Mapping[str, Any]]:
    """Converts DnaChain to an AlphaFold JSON dict."""
    contents = {
        'id': self.id,
        'sequence': self.sequence,
        'modifications': [
            {'modificationType': mod[0], 'basePosition': mod[1]}
            for mod in self.modifications
        ],
    }
    return {'dna': contents}

  def to_ccd_sequence(self) -> Sequence[str]:
    """Converts to a sequence of CCD codes."""
    ccd_coded_seq = [
        residue_names.DNA_COMMON_ONE_TO_TWO.get(res, residue_names.UNK_DNA)
        for res in self.sequence
    ]
    for ccd_code, modification_index in self.modifications:
      ccd_coded_seq[modification_index - 1] = ccd_code
    return ccd_coded_seq


@dataclasses.dataclass(frozen=True, slots=True, kw_only=True)
class Ligand:
  """Ligand input.

  Attributes:
    id: Unique ligand "chain" identifier.
    ccd_ids: The Chemical Component Dictionary or user-defined CCD IDs of the
      chemical components of the ligand. Typically, this is just a single ID,
      but some ligands are composed of multiple components. If that is the case,
      a bond linking these components should be added to the bonded_atom_pairs
      Input field.
    smiles: The SMILES representation of the ligand.
  """

  id: str
  ccd_ids: Sequence[str] | None = None
  smiles: str | None = None

  def __post_init__(self):
    if (self.ccd_ids is None) == (self.smiles is None):
      raise ValueError('Ligand must have one of CCD ID or SMILES set.')

    if self.smiles is not None:
      mol = rd_chem.MolFromSmiles(self.smiles)
      if not mol:
        raise ValueError(f'Unable to make RDKit Mol from SMILES: {self.smiles}')

    # Use hashable types for ccd_ids.
    if self.ccd_ids is not None:
      object.__setattr__(self, 'ccd_ids', tuple(self.ccd_ids))

  @classmethod
  def from_alphafoldserver_dict(
      cls, json_dict: Mapping[str, Any], seq_id: str
  ) -> Self:
    """Constructs Ligand from the AlphaFoldServer JSON dict."""
    # Ligand can be specified either as a ligand, or ion (special-case).
    _validate_keys(json_dict.keys(), {'ligand', 'ion', 'count'})
    if 'ligand' in json_dict:
      return cls(id=seq_id, ccd_ids=[json_dict['ligand'].removeprefix('CCD_')])
    elif 'ion' in json_dict:
      return cls(id=seq_id, ccd_ids=[json_dict['ion']])
    else:
      raise ValueError(f'Unknown ligand type: {json_dict}')

  @classmethod
  def from_dict(
      cls, json_dict: Mapping[str, Any], seq_id: str | None = None
  ) -> Self:
    """Constructs Ligand from the AlphaFold JSON dict."""
    json_dict = json_dict['ligand']
    _validate_keys(json_dict.keys(), {'id', 'ccdCodes', 'smiles'})
    if json_dict.get('ccdCodes') and json_dict.get('smiles'):
      raise ValueError(
          'Ligand cannot have both CCD code and SMILES set at the same time, '
          f'got CCD: {json_dict["ccdCode"]} and SMILES: {json_dict["smiles"]}'
      )

    if 'ccdCodes' in json_dict:
      return cls(id=seq_id or json_dict['id'], ccd_ids=json_dict['ccdCodes'])
    elif 'smiles' in json_dict:
      return cls(id=seq_id or json_dict['id'], smiles=json_dict['smiles'])
    else:
      raise ValueError(f'Unknown ligand type: {json_dict}')

  def to_dict(self) -> Mapping[str, Any]:
    """Converts Ligand to an AlphaFold JSON dict."""
    contents = {'id': self.id}
    if self.ccd_ids is not None:
      contents['ccdCodes'] = self.ccd_ids
    if self.smiles is not None:
      contents['smiles'] = self.smiles
    return {'ligand': contents}


def _sample_rng_seed() -> int:
  """Sample a random seed for AlphaFoldServer job."""
  # See https://alphafoldserver.com/faq#what-are-seeds-and-how-are-they-set.
  return random.randint(0, 2**32 - 1)


@dataclasses.dataclass(frozen=True, slots=True, kw_only=True)
class Input:
  """AlphaFold input.

  Attributes:
    name: The name of the target.
    chains: Protein chains, RNA chains, DNA chains, or ligands.
    protein_chains: Protein chains.
    rna_chains: RNA chains.
    dna_chains: Single strand DNA chains.
    ligands: Ligand (including ion) inputs.
    rng_seeds: Random number generator seeds, one for each model execution.
    bonded_atom_pairs: A list of tuples of atoms that are bonded to each other.
      Each atom is defined by a tuple of (chain_id, res_id, atom_name). Chain
      IDs must be set if there are any bonded atoms. Residue IDs are 1-indexed.
      Atoms in ligands defined by SMILES can't be bonded since SMILES doesn't
      define unique atom names.
    user_ccd: Optional user-defined chemical component dictionary in the CIF
      format. This can be used to provide additional CCD entries that are not
      present in the default CCD and thus define arbitrary new ligands. This is
      more expressive than SMILES since it allows to name all atoms within the
      ligand which in turn makes it possible to define bonds using those atoms.
  """

  name: str
  chains: Sequence[ProteinChain | RnaChain | DnaChain | Ligand]
  rng_seeds: Sequence[int]
  bonded_atom_pairs: Sequence[tuple[BondAtomId, BondAtomId]] | None = None
  user_ccd: str | None = None

  def __post_init__(self):
    if not self.rng_seeds:
      raise ValueError('Input must have at least one RNG seed.')

    if not self.name.strip() or not self.sanitised_name():
      raise ValueError(
          'Input name must be non-empty and contain at least one valid'
          ' character (letters, numbers, dots, dashes, underscores).'
      )

    chain_ids = [c.id for c in self.chains]
    if any(not c.id.isalpha() or c.id.islower() for c in self.chains):
      raise ValueError(
          f'IDs must be alphabetic and upper case, got: {chain_ids}'
      )
    if len(set(chain_ids)) != len(chain_ids):
      raise ValueError('Input JSON contains sequences with duplicate IDs.')

    # Use hashable types for chains, rng_seeds, and bonded_atom_pairs.
    object.__setattr__(self, 'chains', tuple(self.chains))
    object.__setattr__(self, 'rng_seeds', tuple(self.rng_seeds))
    if self.bonded_atom_pairs is not None:
      object.__setattr__(
          self, 'bonded_atom_pairs', tuple(self.bonded_atom_pairs)
      )

  @property
  def protein_chains(self) -> Sequence[ProteinChain]:
    return [chain for chain in self.chains if isinstance(chain, ProteinChain)]

  @property
  def rna_chains(self) -> Sequence[RnaChain]:
    return [chain for chain in self.chains if isinstance(chain, RnaChain)]

  @property
  def dna_chains(self) -> Sequence[DnaChain]:
    return [chain for chain in self.chains if isinstance(chain, DnaChain)]

  @property
  def ligands(self) -> Sequence[Ligand]:
    return [chain for chain in self.chains if isinstance(chain, Ligand)]

  @classmethod
  def from_alphafoldserver_fold_job(cls, fold_job: Mapping[str, Any]) -> Self:
    """Constructs Input from an AlphaFoldServer fold job."""

    # Validate the fold job has the correct format.
    _validate_keys(
        fold_job.keys(),
        {'name', 'modelSeeds', 'sequences', 'dialect', 'version'},
    )
    if 'dialect' not in fold_job and 'version' not in fold_job:
      dialect = ALPHAFOLDSERVER_JSON_DIALECT
      version = ALPHAFOLDSERVER_JSON_VERSION
    elif 'dialect' in fold_job and 'version' in fold_job:
      dialect = fold_job['dialect']
      version = fold_job['version']
    else:
      raise ValueError(
          'AlphaFold Server input JSON must either contain both `dialect` and'
          ' `version` fields, or neither. If neither is specified, it is'
          f' assumed that `dialect="{ALPHAFOLDSERVER_JSON_DIALECT}"` and'
          f' `version="{ALPHAFOLDSERVER_JSON_VERSION}"`.'
      )

    if dialect != ALPHAFOLDSERVER_JSON_DIALECT:
      raise ValueError(
          f'AlphaFold Server input JSON has unsupported dialect: {dialect}, '
          f'expected {ALPHAFOLDSERVER_JSON_DIALECT}.'
      )

    # For now, there is only one AlphaFold Server JSON version.
    if version != ALPHAFOLDSERVER_JSON_VERSION:
      raise ValueError(
          f'AlphaFold Server input JSON has unsupported version: {version}, '
          f'expected {ALPHAFOLDSERVER_JSON_VERSION}.'
      )

    # Parse the chains.
    chains = []
    for sequence in fold_job['sequences']:
      if 'proteinChain' in sequence:
        for _ in range(sequence['proteinChain'].get('count', 1)):
          chains.append(
              ProteinChain.from_alphafoldserver_dict(
                  sequence['proteinChain'],
                  seq_id=mmcif_lib.int_id_to_str_id(len(chains) + 1),
              )
          )
      elif 'rnaSequence' in sequence:
        for _ in range(sequence['rnaSequence'].get('count', 1)):
          chains.append(
              RnaChain.from_alphafoldserver_dict(
                  sequence['rnaSequence'],
                  seq_id=mmcif_lib.int_id_to_str_id(len(chains) + 1),
              )
          )
      elif 'dnaSequence' in sequence:
        for _ in range(sequence['dnaSequence'].get('count', 1)):
          chains.append(
              DnaChain.from_alphafoldserver_dict(
                  sequence['dnaSequence'],
                  seq_id=mmcif_lib.int_id_to_str_id(len(chains) + 1),
              )
          )
      elif 'ion' in sequence:
        for _ in range(sequence['ion'].get('count', 1)):
          chains.append(
              Ligand.from_alphafoldserver_dict(
                  sequence['ion'],
                  seq_id=mmcif_lib.int_id_to_str_id(len(chains) + 1),
              )
          )
      elif 'ligand' in sequence:
        for _ in range(sequence['ligand'].get('count', 1)):
          chains.append(
              Ligand.from_alphafoldserver_dict(
                  sequence['ligand'],
                  seq_id=mmcif_lib.int_id_to_str_id(len(chains) + 1),
              )
          )
      else:
        raise ValueError(f'Unknown sequence type: {sequence}')

    if 'modelSeeds' in fold_job and fold_job['modelSeeds']:
      rng_seeds = [int(seed) for seed in fold_job['modelSeeds']]
    else:
      rng_seeds = [_sample_rng_seed()]

    return cls(name=fold_job['name'], chains=chains, rng_seeds=rng_seeds)

  @classmethod
  def from_json(cls, json_str: str) -> Self:
    """Loads the input from the AlphaFold JSON string."""
    raw_json = json.loads(json_str)

    _validate_keys(
        raw_json.keys(),
        {
            'dialect',
            'version',
            'name',
            'modelSeeds',
            'sequences',
            'bondedAtomPairs',
            'userCCD',
        },
    )

    if 'dialect' not in raw_json or 'version' not in raw_json:
      raise ValueError(
          'AlphaFold 3 input JSON must contain `dialect` and `version` fields.'
      )

    if raw_json['dialect'] != JSON_DIALECT:
      raise ValueError(
          'AlphaFold 3 input JSON has unsupported dialect:'
          f' {raw_json["dialect"]}, expected {JSON_DIALECT}.'
      )

    # For now, there is only one AlphaFold 3 JSON version.
    if raw_json['version'] != JSON_VERSION:
      raise ValueError(
          'AlphaFold 3 input JSON has unsupported version:'
          f' {raw_json["version"]}, expected {JSON_VERSION}.'
      )

    if 'sequences' not in raw_json:
      raise ValueError('AlphaFold 3 input JSON does not contain any sequences.')

    if 'modelSeeds' not in raw_json or not raw_json['modelSeeds']:
      raise ValueError(
          'AlphaFold 3 input JSON must specify at least one rng seed in'
          ' `modelSeeds`.'
      )

    sequences = raw_json['sequences']

    # Make sure sequence IDs are all set.
    raw_sequence_ids = [next(iter(s.values())).get('id') for s in sequences]
    if all(raw_sequence_ids):
      sequence_ids = []
      for sequence_id in raw_sequence_ids:
        if isinstance(sequence_id, list):
          sequence_ids.append(sequence_id)
        else:
          sequence_ids.append([sequence_id])
    else:
      raise ValueError(
          'AlphaFold 3 input JSON contains sequences with unset IDs.'
      )

    flat_seq_ids = []
    for seq_ids in sequence_ids:
      flat_seq_ids.extend(seq_ids)

    chains = []
    for seq_ids, sequence in zip(sequence_ids, sequences, strict=True):
      if len(sequence) != 1:
        raise ValueError(f'Chain {seq_ids} has more than 1 sequence.')
      for seq_id in seq_ids:
        if 'protein' in sequence:
          chains.append(ProteinChain.from_dict(sequence, seq_id=seq_id))
        elif 'rna' in sequence:
          chains.append(RnaChain.from_dict(sequence, seq_id=seq_id))
        elif 'dna' in sequence:
          chains.append(DnaChain.from_dict(sequence, seq_id=seq_id))
        elif 'ligand' in sequence:
          chains.append(Ligand.from_dict(sequence, seq_id=seq_id))
        else:
          raise ValueError(f'Unknown sequence type: {sequence}')

    ligands = [chain for chain in chains if isinstance(chain, Ligand)]
    bonded_atom_pairs = None
    if bonds := raw_json.get('bondedAtomPairs'):
      bonded_atom_pairs = []
      for bond in bonds:
        if len(bond) != 2:
          raise ValueError(f'Bond {bond} must have 2 atoms, got {len(bond)}.')
        bond_beg, bond_end = bond
        if (
            len(bond_beg) != 3
            or not isinstance(bond_beg[0], str)
            or not isinstance(bond_beg[1], int)
            or not isinstance(bond_beg[2], str)
        ):
          raise ValueError(
              f'Atom {bond_beg} in bond {bond} must have 3 components: '
              '(chain_id: str, res_id: int, atom_name: str).'
          )
        if (
            len(bond_end) != 3
            or not isinstance(bond_end[0], str)
            or not isinstance(bond_end[1], int)
            or not isinstance(bond_end[2], str)
        ):
          raise ValueError(
              f'Atom {bond_end} in bond {bond} must have 3 components: '
              '(chain_id: str, res_id: int, atom_name: str).'
          )
        if bond_beg[0] not in flat_seq_ids or bond_end[0] not in flat_seq_ids:
          raise ValueError(f'Invalid chain ID(s) in bond {bond}')
        if bond_beg[1] <= 0 or bond_end[1] <= 0:
          raise ValueError(f'Invalid residue ID(s) in bond {bond}')
        smiles_ligand_ids = set(l.id for l in ligands if l.smiles is not None)
        if bond_beg[0] in smiles_ligand_ids:
          raise ValueError(
              f'Bond {bond} involves an unsupported SMILES ligand {bond_beg[0]}'
          )
        if bond_end[0] in smiles_ligand_ids:
          raise ValueError(
              f'Bond {bond} involves an unsupported SMILES ligand {bond_end[0]}'
          )
        bonded_atom_pairs.append((tuple(bond_beg), tuple(bond_end)))

    return cls(
        name=raw_json['name'],
        chains=chains,
        rng_seeds=[int(seed) for seed in raw_json['modelSeeds']],
        bonded_atom_pairs=bonded_atom_pairs,
        user_ccd=raw_json.get('userCCD'),
    )

  @classmethod
  def from_mmcif(cls, mmcif_str: str, ccd: chemical_components.Ccd) -> Self:
    """Loads the input from an mmCIF string.

    WARNING: Since rng seeds are not stored in mmCIFs, an rng seed is sampled
    in the returned `Input`.

    Args:
      mmcif_str: The mmCIF string.
      ccd: The chemical components dictionary.

    Returns:
      The input in an Input format.
    """

    struc = structure.from_mmcif(
        mmcif_str,
        include_water=False,
        fix_mse_residues=True,
        fix_unknown_dna=True,
        include_bonds=True,
        include_other=False,
    )

    # Create default bioassembly, expanding structures implied by stoichiometry.
    struc = struc.generate_bioassembly(None)

    sequences = struc.chain_single_letter_sequence(
        include_missing_residues=True
    )

    chains = []
    for chain_id, chain_type in zip(
        struc.group_by_chain.chain_id, struc.group_by_chain.chain_type
    ):
      sequence = sequences[chain_id]

      if chain_type in mmcif_names.NON_POLYMER_CHAIN_TYPES:
        residues = list(struc.chain_res_name_sequence()[chain_id])
        if all(ccd.get(res) is not None for res in residues):
          chains.append(Ligand(id=chain_id, ccd_ids=residues))
        elif len(residues) == 1:
          comp_name = residues[0]
          comps = struc.chemical_components_data
          if comps is None:
            raise ValueError(
                'Missing mmCIF chemical components data - this is required for '
                f'a non-CCD ligand {comp_name} defined using SMILES string.'
            )
          chains.append(
              Ligand(id=chain_id, smiles=comps.chem_comp[comp_name].pdbx_smiles)
          )
        else:
          raise ValueError(
              'Multi-component ligand must be defined using CCD IDs, defining'
              ' using SMILES is supported only for single-component ligands. '
              f'Got {residues}'
          )
      else:
        residues = struc.chain_res_name_sequence()[chain_id]
        fixed = struc.chain_res_name_sequence(
            fix_non_standard_polymer_res=True
        )[chain_id]
        modifications = [
            (orig, i + 1)
            for i, (orig, fixed) in enumerate(zip(residues, fixed, strict=True))
            if orig != fixed
        ]

        if chain_type == mmcif_names.PROTEIN_CHAIN:
          chains.append(
              ProteinChain(id=chain_id, sequence=sequence, ptms=modifications)
          )
        elif chain_type == mmcif_names.RNA_CHAIN:
          chains.append(
              RnaChain(
                  id=chain_id, sequence=sequence, modifications=modifications
              )
          )
        elif chain_type == mmcif_names.DNA_CHAIN:
          chains.append(
              DnaChain(
                  id=chain_id, sequence=sequence, modifications=modifications
              )
          )

    bonded_atom_pairs = []
    chain_ids = set(c.id for c in chains)
    for atom_a, atom_b, _ in struc.iter_bonds():
      if atom_a['chain_id'] in chain_ids and atom_b['chain_id'] in chain_ids:
        beg = (atom_a['chain_id'], int(atom_a['res_id']), atom_a['atom_name'])
        end = (atom_b['chain_id'], int(atom_b['res_id']), atom_b['atom_name'])
        bonded_atom_pairs.append((beg, end))

    return cls(
        name=struc.name,
        chains=chains,
        # mmCIFs don't store rng seeds, so we need to sample one here.
        rng_seeds=[_sample_rng_seed()],
        bonded_atom_pairs=bonded_atom_pairs or None,
    )

  def to_structure(self, ccd: chemical_components.Ccd) -> structure.Structure:
    """Converts Input to a Structure.

    WARNING: This method does not preserve the rng seeds.

    Args:
      ccd: The chemical components dictionary.

    Returns:
      The input in a structure.Structure format.
    """
    ids: list[str] = []
    sequences: list[str] = []
    poly_types: list[str] = []
    formats: list[structure.SequenceFormat] = []

    for chain in self.chains:
      ids.append(chain.id)
      match chain:
        case ProteinChain():
          sequences.append('(' + ')('.join(chain.to_ccd_sequence()) + ')')
          poly_types.append(mmcif_names.PROTEIN_CHAIN)
          formats.append(structure.SequenceFormat.CCD_CODES)
        case RnaChain():
          sequences.append('(' + ')('.join(chain.to_ccd_sequence()) + ')')
          poly_types.append(mmcif_names.RNA_CHAIN)
          formats.append(structure.SequenceFormat.CCD_CODES)
        case DnaChain():
          sequences.append('(' + ')('.join(chain.to_ccd_sequence()) + ')')
          poly_types.append(mmcif_names.DNA_CHAIN)
          formats.append(structure.SequenceFormat.CCD_CODES)
        case Ligand():
          if chain.ccd_ids is not None:
            sequences.append('(' + ')('.join(chain.ccd_ids) + ')')
            if len(chain.ccd_ids) == 1:
              poly_types.append(mmcif_names.NON_POLYMER_CHAIN)
            else:
              poly_types.append(mmcif_names.BRANCHED_CHAIN)
            formats.append(structure.SequenceFormat.CCD_CODES)
          elif chain.smiles is not None:
            # Convert to `<unique ligand ID>:<smiles>` format that is expected
            # by structure.from_sequences_and_bonds.
            sequences.append(f'LIG_{chain.id}:{chain.smiles}')
            poly_types.append(mmcif_names.NON_POLYMER_CHAIN)
            formats.append(structure.SequenceFormat.LIGAND_SMILES)
          else:
            raise ValueError('Ligand must have one of CCD ID or SMILES set.')

    # Remap bond chain IDs from chain IDs to chain indices and convert to
    # 0-based residue indexing.
    bonded_atom_pairs = []
    chain_indices = {cid: i for i, cid in enumerate(ids)}
    if self.bonded_atom_pairs is not None:
      for bond_beg, bond_end in self.bonded_atom_pairs:
        bonded_atom_pairs.append((
            (chain_indices[bond_beg[0]], bond_beg[1] - 1, bond_beg[2]),
            (chain_indices[bond_end[0]], bond_end[1] - 1, bond_end[2]),
        ))

    struc = structure.from_sequences_and_bonds(
        sequences=sequences,
        chain_types=poly_types,
        sequence_formats=formats,
        bonded_atom_pairs=bonded_atom_pairs,
        ccd=ccd,
        name=self.sanitised_name(),
        bond_type=mmcif_names.COVALENT_BOND,
        release_date=None,
    )
    # Rename chain IDs to the original ones.
    return struc.rename_chain_ids(dict(zip(struc.chains, ids, strict=True)))

  def to_json(self) -> str:
    """Converts Input to an AlphaFold JSON."""
    alphafold_json = json.dumps(
        {
            'dialect': JSON_DIALECT,
            'version': JSON_VERSION,
            'name': self.name,
            'sequences': [chain.to_dict() for chain in self.chains],
            'modelSeeds': self.rng_seeds,
            'bondedAtomPairs': self.bonded_atom_pairs,
            'userCCD': self.user_ccd,
        },
        indent=2,
    )
    # Remove newlines from the query/template indices arrays. We match the
    # queryIndices/templatesIndices with a non-capturing group. We then match
    # the entire region between the square brackets by looking for lines
    # containing only whitespace, number, or a comma.
    return re.sub(
        r'("(?:queryIndices|templateIndices)": \[)([\s\n\d,]+)(\],?)',
        lambda mtch: mtch[1] + re.sub(r'\n\s+', ' ', mtch[2].strip()) + mtch[3],
        alphafold_json,
    )

  def fill_missing_fields(self) -> Self:
    """Fill missing MSA and template fields with default values."""
    with_missing_fields = [
        c.fill_missing_fields()
        if isinstance(c, (ProteinChain, RnaChain))
        else c
        for c in self.chains
    ]
    return dataclasses.replace(self, chains=with_missing_fields)

  def sanitised_name(self) -> str:
    """Returns sanitised version of the name that can be used as a filename."""
    lower_spaceless_name = self.name.lower().replace(' ', '_')
    allowed_chars = set(string.ascii_lowercase + string.digits + '_-.')
    return ''.join(l for l in lower_spaceless_name if l in allowed_chars)


def check_unique_sanitised_names(fold_inputs: Sequence[Input]) -> None:
  """Checks that the names of the fold inputs are unique."""
  names = [fi.sanitised_name() for fi in fold_inputs]
  if len(set(names)) != len(names):
    raise ValueError(
        f'Fold inputs must have unique sanitised names, got {names}.'
    )


def load_fold_inputs_from_path(json_path: pathlib.Path) -> Sequence[Input]:
  """Loads multiple fold inputs from a JSON string."""
  with open(json_path, 'r') as f:
    json_str = f.read()

  # Parse the JSON string, so we can detect its format.
  raw_json = json.loads(json_str)

  fold_inputs = []
  if isinstance(raw_json, list):
    # AlphaFold Server JSON.
    logging.info(
        'Detected %s is an AlphaFold Server JSON since the top-level is a'
        ' list.',
        json_path,
    )

    logging.info('Loading %d fold jobs from %s', len(raw_json), json_path)
    for fold_job_idx, fold_job in enumerate(raw_json):
      try:
        fold_inputs.append(Input.from_alphafoldserver_fold_job(fold_job))
      except ValueError as e:
        raise ValueError(
            f'Failed to load fold job {fold_job_idx} from {json_path}. The JSON'
            f' at {json_path} was detected to be an AlphaFold Server JSON since'
            ' the top-level is a list.'
        ) from e
  else:
    logging.info(
        'Detected %s is an AlphaFold 3 JSON since the top-level is not a list.',
        json_path,
    )
    # AlphaFold 3 JSON.
    try:
      fold_inputs.append(Input.from_json(json_str))
    except ValueError as e:
      raise ValueError(
          f'Failed to load fold input from {json_path}. The JSON at'
          f' {json_path} was detected to be an AlphaFold 3 JSON since the'
          ' top-level is not a list.'
      ) from e

  check_unique_sanitised_names(fold_inputs)

  return fold_inputs


def load_fold_inputs_from_dir(input_dir: pathlib.Path) -> Sequence[Input]:
  """Loads multiple fold inputs from all JSON files in a given input_dir.

  Args:
    input_dir: The directory containing the JSON files.

  Returns:
    The fold inputs from all JSON files in the input directory.

  Raises:
    ValueError: If the fold inputs have non-unique sanitised names.
  """
  fold_inputs = []
  for file_path in input_dir.glob('*.json'):
    if not file_path.is_file():
      continue

    fold_inputs.extend(load_fold_inputs_from_path(file_path))

  check_unique_sanitised_names(fold_inputs)

  return fold_inputs
