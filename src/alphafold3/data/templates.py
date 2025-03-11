# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""API for retrieving and manipulating template search results."""

from collections.abc import Iterable, Iterator, Mapping, Sequence
import dataclasses
import datetime
import functools
import os
import re
from typing import Any, Final, Self, TypeAlias

from absl import logging
from alphafold3 import structure
from alphafold3.common import resources
from alphafold3.constants import atom_types
from alphafold3.constants import mmcif_names
from alphafold3.constants import residue_names
from alphafold3.data import msa_config
from alphafold3.data import parsers
from alphafold3.data import structure_stores
from alphafold3.data import template_realign
from alphafold3.data.tools import hmmsearch
from alphafold3.structure import mmcif
import numpy as np


_POLYMER_FEATURES: Final[Mapping[str, np.float64 | np.int32 | object]] = {
    'template_aatype': np.int32,
    'template_all_atom_masks': np.float64,
    'template_all_atom_positions': np.float64,
    'template_domain_names': object,
    'template_release_date': object,
    'template_sequence': object,
}

_LIGAND_FEATURES: Final[Mapping[str, Any]] = {
    'ligand_features': Mapping[str, Any]
}


TemplateFeatures: TypeAlias = Mapping[
    str, np.ndarray | bytes | Mapping[str, np.ndarray | bytes]
]
_REQUIRED_METADATA_COLUMNS: Final[Sequence[str]] = (
    'seq_release_date',
    'seq_unresolved_res_num',
    'seq_author_chain_id',
    'seq_sequence',
)


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class _Polymer:
  """Container for alphabet specific (dna, rna, protein) atom information."""

  min_atoms: int
  num_atom_types: int
  atom_order: Mapping[str, int]


_POLYMERS = {
    mmcif_names.PROTEIN_CHAIN: _Polymer(
        min_atoms=5,
        num_atom_types=atom_types.ATOM37_NUM,
        atom_order=atom_types.ATOM37_ORDER,
    ),
    mmcif_names.DNA_CHAIN: _Polymer(
        min_atoms=21,
        num_atom_types=atom_types.ATOM29_NUM,
        atom_order=atom_types.ATOM29_ORDER,
    ),
    mmcif_names.RNA_CHAIN: _Polymer(
        min_atoms=20,
        num_atom_types=atom_types.ATOM29_NUM,
        atom_order=atom_types.ATOM29_ORDER,
    ),
}


def _encode_restype(
    chain_poly_type: str,
    sequence: str,
) -> Sequence[int]:
  """Encodes a sequence of residue names as a sequence of ints.

  Args:
    chain_poly_type: Polymer chain type to determine sequence encoding.
    sequence: Polymer residues. Protein encoded by single letters. RNA and DNA
      encoded by  multi-letter CCD codes.

  Returns:
    A sequence of integers encoding amino acid types for the given chain type.
  """
  if chain_poly_type == mmcif_names.PROTEIN_CHAIN:
    return [
        residue_names.PROTEIN_TYPES_ONE_LETTER_WITH_UNKNOWN_AND_GAP_TO_INT[
            _STANDARDIZED_AA.get(res, res)
        ]
        for res in sequence
    ]

  unk_nucleic = residue_names.UNK_NUCLEIC_ONE_LETTER
  unk_nucleic_idx = residue_names.POLYMER_TYPES_ORDER_WITH_UNKNOWN_AND_GAP[
      unk_nucleic
  ]
  if chain_poly_type == mmcif_names.RNA_CHAIN:
    return [
        residue_names.POLYMER_TYPES_ORDER_WITH_UNKNOWN_AND_GAP.get(
            res, unk_nucleic_idx
        )
        for res in sequence
    ]
  elif chain_poly_type == mmcif_names.DNA_CHAIN:
    # Map UNK DNA to the generic nucleic UNK (N), which happens to also be the
    # same as the RNA UNK.
    return [
        residue_names.POLYMER_TYPES_ORDER_WITH_UNKNOWN_AND_GAP.get(
            residue_names.DNA_COMMON_ONE_TO_TWO.get(res, unk_nucleic),
            unk_nucleic_idx,
        )
        for res in sequence
    ]

  raise NotImplementedError(f'"{chain_poly_type}" unsupported.')


_DAYS_BEFORE_QUERY_DATE: Final[int] = 60
_HIT_DESCRIPTION_REGEX = re.compile(
    r'(?P<pdb_id>[a-z0-9]{4,})_(?P<chain_id>\w+)/(?P<start>\d+)-(?P<end>\d+) '
    r'.* length:(?P<length>\d+)\b.*'
)

_STANDARDIZED_AA = {'B': 'D', 'J': 'X', 'O': 'X', 'U': 'C', 'Z': 'E'}


class Error(Exception):
  """Base class for exceptions."""


class HitDateError(Error):
  """An error indicating that invalid release date was detected."""


class InvalidTemplateError(Error):
  """An error indicating that template is invalid."""


@dataclasses.dataclass(frozen=True, kw_only=True)
class Hit:
  """Template hit metrics derived from the MSA for filtering and featurising.

  Attributes:
    pdb_id: The PDB ID of the hit.
    auth_chain_id: The author chain ID of the hit.
    hmmsearch_sequence: Hit sequence as given in hmmsearch a3m output.
    structure_sequence: Hit sequence as given in PDB structure.
    unresolved_res_indices: Indices of unresolved residues in the structure
      sequence. 0-based.
    query_sequence: The query nucleotide/amino acid sequence.
    start_index: The start index of the sequence relative to the full PDB seqres
      sequence. Inclusive and uses 0-based indexing.
    end_index: The end index of the sequence relative to the full PDB seqres
      sequence. Exclusive and uses 0-based indexing.
    full_length: Length of the full PDB seqres sequence. This can be different
      from the length from the actual sequence we get from the mmCIF and we use
      this to detect whether we need to realign or not.
    release_date: The release date of the PDB corresponding to this hit.
    chain_poly_type: The polymer type of the selected hit structure.
  """

  pdb_id: str
  auth_chain_id: str
  hmmsearch_sequence: str
  structure_sequence: str
  unresolved_res_indices: Sequence[int] | None
  query_sequence: str
  start_index: int
  end_index: int
  full_length: int
  release_date: datetime.date
  chain_poly_type: str

  @functools.cached_property
  def query_to_hit_mapping(self) -> Mapping[int, int]:
    """0-based query index to hit index mapping."""
    query_to_hit_mapping = {}
    hit_index = 0
    query_index = 0
    for residue in self.hmmsearch_sequence:
      # Gap inserted in the template
      if residue == '-':
        query_index += 1
      # Deleted residue in the template (would be a gap in the query).
      elif residue.islower():
        hit_index += 1
      # Normal aligned residue, in both query and template. Add to mapping.
      elif residue.isupper():
        query_to_hit_mapping[query_index] = hit_index
        query_index += 1
        hit_index += 1

    structure_subseq = self.structure_sequence[
        self.start_index : self.end_index
    ]
    if self.matching_sequence != structure_subseq:
      # The seqres sequence doesn't match the structure sequence. Two cases:
      # 1. The sequences have the same length. The sequences are different
      #    because our 3->1 residue code mapping is different from the one PDB
      #    uses. We don't do anything in this case as both sequences have the
      #    same length, so the original query to hit mapping stays valid.
      # 2. The sequences don't have the same length, the one in structure is
      #    shorter. In this case we change the mapping to match the actual
      #    structure sequence using a simple realignment algorithm.
      # This procedure was validated on all PDB seqres (2023_01_12) sequences
      # and handles all cases that can happen.
      if self.full_length != len(self.structure_sequence):
        return template_realign.realign_hit_to_structure(
            hit_sequence=self.matching_sequence,
            hit_start_index=self.start_index,
            hit_end_index=self.end_index,
            full_length=self.full_length,
            structure_sequence=self.structure_sequence,
            query_to_hit_mapping=query_to_hit_mapping,
        )

    # Hmmsearch returns a subsequence and so far indices have been relative to
    # the subsequence. Add an offset to index relative to the full structure
    # sequence.
    return {q: h + self.start_index for q, h in query_to_hit_mapping.items()}

  @property
  def matching_sequence(self) -> str:
    """Returns the matching hit sequence including insertions.

    Make deleted residues uppercase and remove gaps ("-").
    """
    return self.hmmsearch_sequence.upper().replace('-', '')

  @functools.cached_property
  def output_templates_sequence(self) -> str:
    """Returns the final template sequence."""
    result_seq = ['-'] * len(self.query_sequence)
    for query_index, template_index in self.query_to_hit_mapping.items():
      result_seq[query_index] = self.structure_sequence[template_index]
    return ''.join(result_seq)

  @property
  def length_ratio(self) -> float:
    """Ratio of the length of the hit sequence to the query."""
    return len(self.matching_sequence) / len(self.query_sequence)

  @property
  def align_ratio(self) -> float:
    """Ratio of the number of aligned residues to the query length."""
    return len(self.query_to_hit_mapping) / len(self.query_sequence)

  @functools.cached_property
  def is_valid(self) -> bool:
    """Whether hit can be used as a template."""
    if self.unresolved_res_indices is None:
      return False

    return bool(
        set(self.query_to_hit_mapping.values())
        - set(self.unresolved_res_indices)
    )

  @property
  def full_name(self) -> str:
    """A full name of the hit."""
    return f'{self.pdb_id}_{self.auth_chain_id}'

  def __post_init__(self):
    if not self.pdb_id.islower() and not self.pdb_id.isdigit():
      raise ValueError(f'pdb_id must be lowercase {self.pdb_id}')

    if not (0 <= self.start_index <= self.end_index):
      raise ValueError(
          'Start must be non-negative and less than or equal to end index. '
          f'Range: {self.start_index}-{self.end_index}'
      )

    if len(self.matching_sequence) != (self.end_index - self.start_index):
      raise ValueError(
          'Sequence length must be equal to end_index - start_index. '
          f'{len(self.matching_sequence)} != {self.end_index} - '
          f'{self.start_index}'
      )

    if self.full_length < 0:
      raise ValueError(f'Full length must be non-negative: {self.full_length}')

  def keep(
      self,
      *,
      release_date_cutoff: datetime.date | None,
      max_subsequence_ratio: float | None,
      min_hit_length: int | None,
      min_align_ratio: float | None,
  ) -> bool:
    """Returns whether the hit should be kept.

    In addition to filtering on all of the provided parameters, this method also
    excludes hits with unresolved residues.

    Args:
      release_date_cutoff: Maximum release date of the template.
      max_subsequence_ratio: If set, excludes hits which are an exact
        subsequence of the query sequence, and longer than this ratio. Useful to
        avoid ground truth leakage.
      min_hit_length: If set, excludes hits which have fewer residues than this.
      min_align_ratio: If set, excludes hits where the number of residues
        aligned to the query is less than this proportion of the template
        length.
    """
    # Exclude hits which are too recent.
    if (
        release_date_cutoff is not None
        and self.release_date > release_date_cutoff
    ):
      return False

    # Exclude hits which are large duplicates of the query_sequence.
    if (
        max_subsequence_ratio is not None
        and self.length_ratio > max_subsequence_ratio
    ):
      if self.matching_sequence in self.query_sequence:
        return False

    # Exclude hits which are too short.
    if (
        min_hit_length is not None
        and len(self.matching_sequence) < min_hit_length
    ):
      return False

    # Exclude hits with unresolved residues.
    if not self.is_valid:
      return False

    # Exclude hits with too few alignments.
    try:
      if min_align_ratio is not None and self.align_ratio <= min_align_ratio:
        return False
    except template_realign.AlignmentError as e:
      logging.warning('Failed to align %s: %s', self, str(e))
      return False

    return True


def _filter_hits(
    hits: Iterable[Hit],
    release_date_cutoff: datetime.date,
    max_subsequence_ratio: float | None,
    min_align_ratio: float | None,
    min_hit_length: int | None,
    deduplicate_sequences: bool,
    max_hits: int | None,
) -> Sequence[Hit]:
  """Filters hits based on the filter config."""
  filtered_hits = []
  seen_before = set()
  for hit in hits:
    if not hit.keep(
        max_subsequence_ratio=max_subsequence_ratio,
        min_align_ratio=min_align_ratio,
        min_hit_length=min_hit_length,
        release_date_cutoff=release_date_cutoff,
    ):
      continue

    # Remove duplicate templates, keeping the first.
    if deduplicate_sequences:
      if hit.output_templates_sequence in seen_before:
        continue
      seen_before.add(hit.output_templates_sequence)

    filtered_hits.append(hit)
    if max_hits and len(filtered_hits) == max_hits:
      break

  return filtered_hits


@dataclasses.dataclass(init=False)
class Templates:
  """A container for templates that were found for the given query sequence.

  The structure_store is constructed from the config by default. Callers can
  optionally supply a structure_store to the constructor to avoid the cost of
  construction and metadata loading.
  """

  def __init__(
      self,
      *,
      query_sequence: str,
      hits: Sequence[Hit],
      max_template_date: datetime.date,
      structure_store: structure_stores.StructureStore,
      query_release_date: datetime.date | None = None,
  ):
    self._query_sequence = query_sequence
    self._hits = tuple(hits)
    self._max_template_date = max_template_date
    self._query_release_date = query_release_date
    self._hit_structures = {}
    self._structure_store = structure_store

    if any(h.query_sequence != self._query_sequence for h in self.hits):
      raise ValueError('All hits must match the query sequence.')

    if self._hits:
      chain_poly_type = self._hits[0].chain_poly_type
      if any(h.chain_poly_type != chain_poly_type for h in self.hits):
        raise ValueError('All hits must have the same chain_poly_type.')

  @classmethod
  def from_seq_and_a3m(
      cls,
      *,
      query_sequence: str,
      msa_a3m: str,
      max_template_date: datetime.date,
      database_path: os.PathLike[str] | str,
      hmmsearch_config: msa_config.HmmsearchConfig,
      max_a3m_query_sequences: int | None,
      structure_store: structure_stores.StructureStore,
      filter_config: msa_config.TemplateFilterConfig | None = None,
      query_release_date: datetime.date | None = None,
      chain_poly_type: str = mmcif_names.PROTEIN_CHAIN,
  ) -> Self:
    """Creates templates from a run of hmmsearch tool against a custom a3m.

    Args:
      query_sequence: The polymer sequence of the target query.
      msa_a3m: An a3m of related polymers aligned to the query sequence, this is
        used to create an HMM for the hmmsearch run.
      max_template_date: This is used to filter templates for training, ensuring
        that they do not leak ground truth information used in testing sets.
      database_path: A path to the sequence database to search for templates.
      hmmsearch_config: Config with Hmmsearch settings.
      max_a3m_query_sequences: The maximum number of input MSA sequences to use
        to construct the profile which is then used to search for templates.
      structure_store: Structure store to fetch template structures from.
      filter_config: Optional config that controls which and how many hits to
        keep. More performant than constructing and then filtering. If not
        provided, no filtering is done.
      query_release_date: The release_date of the template query, this is used
        to filter templates for training, ensuring that they do not leak
        structure information from the future.
      chain_poly_type: The polymer type of the templates.

    Returns:
      Templates object containing a list of Hits initialised from the
      structure_store metadata and a3m alignments.
    """
    hmmsearch_a3m = run_hmmsearch_with_a3m(
        database_path=database_path,
        hmmsearch_config=hmmsearch_config,
        max_a3m_query_sequences=max_a3m_query_sequences,
        a3m=msa_a3m,
    )
    return cls.from_hmmsearch_a3m(
        query_sequence=query_sequence,
        a3m=hmmsearch_a3m,
        max_template_date=max_template_date,
        query_release_date=query_release_date,
        chain_poly_type=chain_poly_type,
        structure_store=structure_store,
        filter_config=filter_config,
    )

  @classmethod
  def from_hmmsearch_a3m(
      cls,
      *,
      query_sequence: str,
      a3m: str,
      max_template_date: datetime.date,
      structure_store: structure_stores.StructureStore,
      filter_config: msa_config.TemplateFilterConfig | None = None,
      query_release_date: datetime.date | None = None,
      chain_poly_type: str = mmcif_names.PROTEIN_CHAIN,
  ) -> Self:
    """Creates Templates from a Hmmsearch A3M.

    Args:
      query_sequence: The polymer sequence of the target query.
      a3m: Results of Hmmsearch in A3M format. This provides a list of potential
        template alignments and pdb codes.
      max_template_date: This is used to filter templates for training, ensuring
        that they do not leak ground truth information used in testing sets.
      structure_store: Structure store to fetch template structures from.
      filter_config: Optional config that controls which and how many hits to
        keep. More performant than constructing and then filtering. If not
        provided, no filtering is done.
      query_release_date: The release_date of the template query, this is used
        to filter templates for training, ensuring that they do not leak
        structure information from the future.
      chain_poly_type: The polymer type of the templates.

    Returns:
      Templates object containing a list of Hits initialised from the
      structure_store metadata and a3m alignments.
    """

    def hit_generator(a3m: str):
      if not a3m:
        return  # Hmmsearch could return an empty string if there are no hits.

      for hit_seq, hit_desc in parsers.lazy_parse_fasta_string(a3m):
        pdb_id, auth_chain_id, start, end, full_length = _parse_hit_description(
            hit_desc
        )

        release_date, sequence, unresolved_res_ids = _parse_hit_metadata(
            structure_store, pdb_id, auth_chain_id
        )
        if unresolved_res_ids is None:
          continue

        # seq_unresolved_res_num are 1-based, setting to 0-based indices.
        unresolved_indices = [i - 1 for i in unresolved_res_ids]

        yield Hit(
            pdb_id=pdb_id,
            auth_chain_id=auth_chain_id,
            hmmsearch_sequence=hit_seq,
            structure_sequence=sequence,
            query_sequence=query_sequence,
            unresolved_res_indices=unresolved_indices,
            start_index=start - 1,  # Raw value is residue number, not index.
            end_index=end,
            full_length=full_length,
            release_date=datetime.date.fromisoformat(release_date),
            chain_poly_type=chain_poly_type,
        )

    if filter_config is None:
      hits = tuple(hit_generator(a3m))
    else:
      hits = _filter_hits(
          hit_generator(a3m),
          release_date_cutoff=filter_config.max_template_date,
          max_subsequence_ratio=filter_config.max_subsequence_ratio,
          min_align_ratio=filter_config.min_align_ratio,
          min_hit_length=filter_config.min_hit_length,
          deduplicate_sequences=filter_config.deduplicate_sequences,
          max_hits=filter_config.max_hits,
      )

    return Templates(
        query_sequence=query_sequence,
        query_release_date=query_release_date,
        hits=hits,
        max_template_date=max_template_date,
        structure_store=structure_store,
    )

  @property
  def query_sequence(self) -> str:
    return self._query_sequence

  @property
  def hits(self) -> tuple[Hit, ...]:
    return self._hits

  @property
  def query_release_date(self) -> datetime.date | None:
    return self._query_release_date

  @property
  def num_hits(self) -> int:
    return len(self._hits)

  @functools.cached_property
  def release_date_cutoff(self) -> datetime.date:
    if self.query_release_date is None:
      return self._max_template_date
    return min(
        self._max_template_date,
        self.query_release_date
        - datetime.timedelta(days=_DAYS_BEFORE_QUERY_DATE),
    )

  def __repr__(self) -> str:
    return f'Templates({self.num_hits} hits)'

  def filter(
      self,
      *,
      max_subsequence_ratio: float | None,
      min_align_ratio: float | None,
      min_hit_length: int | None,
      deduplicate_sequences: bool,
      max_hits: int | None,
  ) -> Self:
    """Returns a new Templates object with only the hits that pass all filters.

    This also filters on query_release_date and max_template_date.

    Args:
      max_subsequence_ratio: If set, excludes hits which are an exact
        subsequence of the query sequence, and longer than this ratio. Useful to
        avoid ground truth leakage.
      min_align_ratio: If set, excludes hits where the number of residues
        aligned to the query is less than this proportion of the template
        length.
      min_hit_length: If set, excludes hits which have fewer residues than this.
      deduplicate_sequences: Whether to exclude duplicate template sequences,
        keeping only the first. This can be useful in increasing the diversity
        of hits especially in the case of homomer hits.
      max_hits: If set, excludes any hits which exceed this count.
    """
    filtered_hits = _filter_hits(
        hits=self._hits,
        release_date_cutoff=self.release_date_cutoff,
        max_subsequence_ratio=max_subsequence_ratio,
        min_align_ratio=min_align_ratio,
        min_hit_length=min_hit_length,
        deduplicate_sequences=deduplicate_sequences,
        max_hits=max_hits,
    )
    return Templates(
        query_sequence=self.query_sequence,
        query_release_date=self.query_release_date,
        hits=filtered_hits,
        max_template_date=self._max_template_date,
        structure_store=self._structure_store,
    )

  def get_hits_with_structures(
      self,
  ) -> Sequence[tuple[Hit, structure.Structure]]:
    """Returns hits + Structures, Structures filtered to the hit's chain."""
    results = []
    structures = {struc.name.lower(): struc for struc in self.structures}
    for hit in self.hits:
      if not hit.is_valid:
        raise InvalidTemplateError(
            'Hits must be filtered before calling get_hits_with_structures.'
        )
      struc = structures[hit.pdb_id]
      label_chain_id = struc.polymer_auth_asym_id_to_label_asym_id().get(
          hit.auth_chain_id
      )
      results.append((hit, struc.filter(chain_id=label_chain_id)))
    return results

  def featurize(
      self,
      include_ligand_features: bool = True,
  ) -> TemplateFeatures:
    """Featurises the templates and returns a map of feature names to features.

    NB: If you don't do any prefiltering, this method might be slow to run
    as it has to fetch many CIFs and featurize them all.

    Args:
      include_ligand_features: Whether to compute ligand features.

    Returns:
      Template features: A mapping of template feature labels to features, which
        may be numpy arrays, bytes objects, or for the special case of label
        `ligand_features` (if `include_ligand_features` is True), a nested
        feature map of labels to numpy arrays.

    Raises:
      InvalidTemplateError: If hits haven't been filtered before featurization.
    """
    hits_by_pdb_id = {}
    for idx, hit in enumerate(self.hits):
      if not hit.is_valid:
        raise InvalidTemplateError(
            f'Hits must be filtered before featurizing, got unprocessed {hit=}'
        )
      hits_by_pdb_id.setdefault(hit.pdb_id, []).append((idx, hit))

    unsorted_features = []
    for struc in self.structures:
      pdb_id = str(struc.name).lower()
      for idx, hit in hits_by_pdb_id[pdb_id]:
        try:
          label_chain_id = struc.polymer_auth_asym_id_to_label_asym_id()[
              hit.auth_chain_id
          ]
          hit_features = {
              **get_polymer_features(
                  chain=struc.filter(chain_id=label_chain_id),
                  chain_poly_type=hit.chain_poly_type,
                  query_sequence_length=len(hit.query_sequence),
                  query_to_hit_mapping=hit.query_to_hit_mapping,
              ),
          }
          if include_ligand_features:
            hit_features['ligand_features'] = _get_ligand_features(struc)
          unsorted_features.append((idx, hit_features))
        except Error as e:
          raise type(e)(f'Failed to featurise {hit=}') from e

    sorted_features = sorted(unsorted_features, key=lambda x: x[0])
    sorted_features = [feat for _, feat in sorted_features]
    return package_template_features(
        hit_features=sorted_features,
        include_ligand_features=include_ligand_features,
    )

  @property
  def structures(self) -> Iterator[structure.Structure]:
    """Yields template structures for each unique PDB ID among hits.

    If there are multiple hits in the same Structure, the Structure will be
    included only once by this method.

    Yields:
      A Structure object for each unique PDB ID among hits.

    Raises:
      HitDateError: If template's release date exceeds max cutoff date.
    """

    for hit in self.hits:
      if hit.release_date > self.release_date_cutoff:  # pylint: disable=comparison-with-callable
        raise HitDateError(
            f'Invalid release date for hit {hit.pdb_id=}, when release date '
            f'cutoff is {self.release_date_cutoff}.'
        )

    # Get the set of pdbs to load. In particular, remove duplicate PDB IDs.
    targets_to_load = tuple({hit.pdb_id for hit in self.hits})

    for target_name in targets_to_load:
      yield structure.from_mmcif(
          mmcif_string=self._structure_store.get_mmcif_str(target_name),
          fix_mse_residues=True,
          fix_arginines=True,
          include_water=False,
          include_bonds=False,
          include_other=True,  # For non-standard polymer chains.
      )


def _parse_hit_description(description: str) -> tuple[str, str, int, int, int]:
  """Parses the hmmsearch A3M sequence description line."""
  # Example lines (protein, nucleic, no description):
  # >4pqx_A/2-217 [subseq from] mol:protein length:217  Free text
  # >4pqx_A/2-217 [subseq from] mol:na length:217  Free text
  # >5g3r_A/1-55 [subseq from] mol:protein length:352
  if match := re.fullmatch(_HIT_DESCRIPTION_REGEX, description):
    return (
        match['pdb_id'],
        match['chain_id'],
        int(match['start']),
        int(match['end']),
        int(match['length']),
    )
  else:
    raise ValueError(f'Could not parse description "{description}"')


def _parse_hit_metadata(
    structure_store: structure_stores.StructureStore,
    pdb_id: str,
    auth_chain_id: str,
) -> tuple[Any, str | None, Sequence[int] | None]:
  """Parse hit metadata by parsing mmCIF from structure store."""
  try:
    cif = mmcif.from_string(structure_store.get_mmcif_str(pdb_id))
  except structure_stores.NotFoundError:
    logging.warning(
        'Failed to get mmCIF for %s (author chain %s).', pdb_id, auth_chain_id
    )
    return None, None, None
  release_date = mmcif.get_release_date(cif)

  try:
    struc = structure.from_parsed_mmcif(
        cif,
        model_id=structure.ModelID.ALL,
        include_water=True,
        include_other=True,
        include_bonds=False,
    )
  except ValueError:
    struc = structure.from_parsed_mmcif(
        cif,
        model_id=structure.ModelID.FIRST,
        include_water=True,
        include_other=True,
        include_bonds=False,
    )

  sequence = struc.polymer_author_chain_single_letter_sequence(
      include_missing_residues=True,
      protein=True,
      dna=True,
      rna=True,
      other=True,
  )[auth_chain_id]

  unresolved_res_ids = struc.filter(
      chain_auth_asym_id=auth_chain_id
  ).unresolved_residues.id

  return release_date, sequence, unresolved_res_ids


def get_polymer_features(
    *,
    chain: structure.Structure,
    chain_poly_type: str,
    query_sequence_length: int,
    query_to_hit_mapping: Mapping[int, int],
) -> Mapping[str, Any]:
  """Returns features for this polymer chain.

  Args:
    chain: Structure object representing the template. Must be already filtered
      to a single chain.
    chain_poly_type: The chain polymer type (protein, DNA, RNA).
    query_sequence_length: The length of the query sequence.
    query_to_hit_mapping: 0-based query index to hit index mapping.

  Returns:
    A dictionary with polymer features for template_chain_id in the struc.

  Raises:
    ValueError: If the input structure contains more than just a single chain.
  """
  if len(chain.polymer_auth_asym_id_to_label_asym_id()) != 1:
    raise ValueError('The structure must be filtered to a single chain.')

  if chain.name is None:
    raise ValueError('The structure must have a name.')

  if chain.release_date is None:
    raise ValueError('The structure must have a release date.')

  auth_chain_id, label_chain_id = next(
      iter(chain.polymer_auth_asym_id_to_label_asym_id().items())
  )
  chain_sequence = chain.chain_single_letter_sequence()[label_chain_id]

  polymer = _POLYMERS[chain_poly_type]
  positions, positions_mask = chain.to_res_arrays(
      include_missing_residues=True, atom_order=polymer.atom_order
  )
  template_all_atom_positions = np.zeros(
      (query_sequence_length, polymer.num_atom_types, 3), dtype=np.float64
  )
  template_all_atom_masks = np.zeros(
      (query_sequence_length, polymer.num_atom_types), dtype=np.int64
  )

  template_sequence = ['-'] * query_sequence_length
  for query_index, template_index in query_to_hit_mapping.items():
    template_all_atom_positions[query_index] = positions[template_index]
    template_all_atom_masks[query_index] = positions_mask[template_index]
    template_sequence[query_index] = chain_sequence[template_index]

  template_sequence = ''.join(template_sequence)
  template_aatype = _encode_restype(chain_poly_type, template_sequence)
  template_name = f'{chain.name.lower()}_{auth_chain_id}'
  release_date = chain.release_date.strftime('%Y-%m-%d')
  return {
      'template_all_atom_positions': template_all_atom_positions,
      'template_all_atom_masks': template_all_atom_masks,
      'template_sequence': template_sequence.encode(),
      'template_aatype': np.array(template_aatype, dtype=np.int32),
      'template_domain_names': np.array(template_name.encode(), dtype=object),
      'template_release_date': np.array(release_date.encode(), dtype=object),
  }


def _get_ligand_features(
    struc: structure.Structure,
) -> Mapping[str, Mapping[str, np.ndarray | bytes]]:
  """Returns features for the ligands in this structure."""
  ligand_struc = struc.filter_to_entity_type(ligand=True)
  assert ligand_struc.coords is not None
  assert ligand_struc.atom_name is not None
  assert ligand_struc.atom_occupancy is not None

  ligand_features = {}
  for ligand_chain_id in ligand_struc.chains:
    idxs = np.where(ligand_struc.chain_id == ligand_chain_id)[0]
    if idxs.shape[0]:
      ligand_features[ligand_chain_id] = {
          'ligand_atom_positions': ligand_struc.coords[idxs, :].astype(
              np.float32
          ),
          'ligand_atom_names': ligand_struc.atom_name[idxs].astype(object),
          'ligand_atom_occupancies': ligand_struc.atom_occupancy[idxs].astype(
              np.float32
          ),
          'ccd_id': ligand_struc.res_name[idxs][0].encode(),
      }
  return ligand_features


def package_template_features(
    *,
    hit_features: Sequence[Mapping[str, Any]],
    include_ligand_features: bool,
) -> Mapping[str, Any]:
  """Stacks polymer features, adds empty and keeps ligand features unstacked."""

  features_to_include = set(_POLYMER_FEATURES)
  if include_ligand_features:
    features_to_include.update(_LIGAND_FEATURES)

  features = {
      feat: [single_hit_features[feat] for single_hit_features in hit_features]
      for feat in features_to_include
  }

  stacked_features = {}
  for k, v in features.items():
    if k in _POLYMER_FEATURES:
      v = np.stack(v, axis=0) if v else np.array([], dtype=_POLYMER_FEATURES[k])
    stacked_features[k] = v

  return stacked_features


def _resolve_path(path: os.PathLike[str] | str) -> str:
  """Resolves path for data dep paths, stringifies otherwise."""
  # Data dependency paths: db baked into the binary.
  resolved_path = resources.filename(path)
  if os.path.exists(resolved_path):
    return resolved_path
  else:
    # Other paths, e.g. local.
    return str(path)


def run_hmmsearch_with_a3m(
    *,
    database_path: os.PathLike[str] | str,
    hmmsearch_config: msa_config.HmmsearchConfig,
    max_a3m_query_sequences: int | None,
    a3m: str | None,
) -> str:
  """Runs Hmmsearch to get a3m string of hits."""
  searcher = hmmsearch.Hmmsearch(
      binary_path=hmmsearch_config.hmmsearch_binary_path,
      hmmbuild_binary_path=hmmsearch_config.hmmbuild_binary_path,
      database_path=_resolve_path(database_path),
      e_value=hmmsearch_config.e_value,
      inc_e=hmmsearch_config.inc_e,
      dom_e=hmmsearch_config.dom_e,
      incdom_e=hmmsearch_config.incdom_e,
      alphabet=hmmsearch_config.alphabet,
      filter_f1=hmmsearch_config.filter_f1,
      filter_f2=hmmsearch_config.filter_f2,
      filter_f3=hmmsearch_config.filter_f3,
      filter_max=hmmsearch_config.filter_max,
  )
  # STO enables us to annotate query non-gap columns as reference columns.
  sto = parsers.convert_a3m_to_stockholm(a3m, max_a3m_query_sequences)
  return searcher.query_with_sto(sto, model_construction='hand')
