# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Structure class for representing and processing molecular structures."""

import collections
from collections.abc import Callable, Collection, Iterable, Iterator, Mapping, Sequence, Set
import dataclasses
import datetime
import enum
import functools
import itertools
import typing
from typing import Any, ClassVar, Final, Literal, NamedTuple, Self, TypeAlias, TypeVar

from alphafold3.constants import atom_types
from alphafold3.constants import chemical_components
from alphafold3.constants import mmcif_names
from alphafold3.constants import residue_names
from alphafold3.cpp import membership
from alphafold3.cpp import string_array
from alphafold3.structure import bioassemblies
from alphafold3.structure import chemical_components as struc_chem_comps
from alphafold3.structure import mmcif
from alphafold3.structure import structure_tables
from alphafold3.structure import table
import numpy as np

# Controls the default number of decimal places for coordinates when writing to
# mmCIF.
_COORDS_DECIMAL_PLACES: Final[int] = 3


@enum.unique
class CascadeDelete(enum.Enum):
  NONE = 0
  FULL = 1
  CHAINS = 2


# See www.python.org/dev/peps/pep-0484/#support-for-singleton-types-in-unions
class _UnsetSentinel(enum.Enum):
  UNSET = object()


_UnsetType = Literal[_UnsetSentinel.UNSET]
_UNSET = _UnsetSentinel.UNSET


class Bond(NamedTuple):
  """Describes a bond between two atoms."""

  from_atom: Mapping[str, str | int | float | np.ndarray]
  dest_atom: Mapping[str, str | int | float | np.ndarray]
  bond_info: Mapping[str, str | int]


class MissingAtomError(Exception):
  """Error raised when an atom is missing during alignment."""


class MissingAuthorResidueIdError(Exception):
  """Raised when author naming data is missing for a residue.

  This can occur in certain edge cases where missing residue data is provided
  without also providing author IDs for those missing residues.
  """


# AllResidues is a mapping from label_asym_id to a sequence of (label_comp_id,
# label_seq_id) pairs. These represent the full sequence including residues
# that might be missing (e.g. unresolved residues in X-ray data).
AllResidues: TypeAlias = Mapping[str, Sequence[tuple[str, int]]]
AuthorNamingScheme: TypeAlias = structure_tables.AuthorNamingScheme


# External residue ID given to missing residues that don't have an ID
# already provided. In mmCIFs this data is found in _pdbx_poly_seq_scheme.
MISSING_AUTH_SEQ_ID: Final[str] = '.'


# Maps from structure fields to column names in the relevant table.
CHAIN_FIELDS: Final[Mapping[str, str]] = {
    'chain_id': 'id',
    'chain_type': 'type',
    'chain_auth_asym_id': 'auth_asym_id',
    'chain_entity_id': 'entity_id',
    'chain_entity_desc': 'entity_desc',
}


RESIDUE_FIELDS: Final[Mapping[str, str]] = {
    'res_id': 'id',
    'res_name': 'name',
    'res_auth_seq_id': 'auth_seq_id',
    'res_insertion_code': 'insertion_code',
}

ATOM_FIELDS: Final[Mapping[str, str]] = {
    'atom_name': 'name',
    'atom_element': 'element',
    'atom_x': 'x',
    'atom_y': 'y',
    'atom_z': 'z',
    'atom_b_factor': 'b_factor',
    'atom_occupancy': 'occupancy',
    'atom_key': 'key',
}

# Fields in structure.
ARRAY_FIELDS = frozenset({
    'atom_b_factor',
    'atom_element',
    'atom_key',
    'atom_name',
    'atom_occupancy',
    'atom_x',
    'atom_y',
    'atom_z',
    'chain_id',
    'chain_type',
    'res_id',
    'res_name',
})

GLOBAL_FIELDS = frozenset({
    'name',
    'release_date',
    'resolution',
    'structure_method',
    'bioassembly_data',
    'chemical_components_data',
})

# Fields which can be updated in copy_and_update.
_UPDATEABLE_FIELDS: Final[Set[str]] = frozenset({
    'all_residues',
    'atom_b_factor',
    'atom_element',
    'atom_key',
    'atom_name',
    'atom_occupancy',
    'atom_x',
    'atom_y',
    'atom_z',
    'bioassembly_data',
    'bonds',
    'chain_id',
    'chain_type',
    'chemical_components_data',
    'name',
    'release_date',
    'res_id',
    'res_name',
    'resolution',
    'structure_method',
})


def fix_non_standard_polymer_residues(
    res_names: np.ndarray, chain_type: str
) -> np.ndarray:
  """Remaps residue names to the closest standard protein/RNA/DNA residue.

  If residue name is already a standard type, it is not altered.
  If a match cannot be found, returns 'UNK' for protein chainresidues and 'N'
    for RNA/DNA chain residue.

  Args:
     res_names: A numpy array of string residue names (CCD monomer codes). E.g.
       'ARG' (protein), 'DT' (DNA), 'N' (RNA).
     chain_type: The type of the chain, must be PROTEIN_CHAIN, RNA_CHAIN or
       DNA_CHAIN.

  Returns:
    An array remapped so that its elements are all from
    PROTEIN_TYPES_WITH_UNKNOWN | RNA_TYPES | DNA_TYPES | {'N'}.

  Raises:
    ValueError: If chain_type not in PEPTIDE_CHAIN_TYPES or
      {OTHER_CHAIN, RNA_CHAIN, DNA_CHAIN, DNA_RNA_HYBRID_CHAIN}.
  """
  # Map to one letter code, then back to common res_names.
  one_letter_codes = string_array.remap(
      res_names, mapping=residue_names.CCD_NAME_TO_ONE_LETTER, default_value='X'
  )

  if (
      chain_type in mmcif_names.PEPTIDE_CHAIN_TYPES
      or chain_type == mmcif_names.OTHER_CHAIN
  ):
    mapping = residue_names.PROTEIN_COMMON_ONE_TO_THREE
    default_value = 'UNK'
  elif chain_type == mmcif_names.RNA_CHAIN:
    # RNA has single-letter CCD monomer codes.
    mapping = {r: r for r in residue_names.RNA_TYPES}
    default_value = 'N'
  elif chain_type == mmcif_names.DNA_CHAIN:
    mapping = residue_names.DNA_COMMON_ONE_TO_TWO
    default_value = 'N'
  elif chain_type == mmcif_names.DNA_RNA_HYBRID_CHAIN:
    mapping = {r: r for r in residue_names.NUCLEIC_TYPES_WITH_UNKNOWN}
    default_value = 'N'
  else:
    raise ValueError(f'Expected a protein/DNA/RNA chain but got {chain_type}')

  return string_array.remap(
      one_letter_codes, mapping=mapping, default_value=default_value
  )


def _get_change_indices(arr: np.ndarray) -> np.ndarray:
  if arr.size == 0:
    return np.array([], dtype=np.int32)
  else:
    changing_idxs = np.where(arr[1:] != arr[:-1])[0] + 1
    return np.concatenate(([0], changing_idxs), axis=0)


def _unpack_filter_predicates(
    predicate_by_field_name: Mapping[str, table.FilterPredicate],
) -> tuple[
    Mapping[str, table.FilterPredicate],
    Mapping[str, table.FilterPredicate],
    Mapping[str, table.FilterPredicate],
]:
  """Unpacks filter kwargs into predicates for each table."""
  chain_predicates = {}
  res_predicates = {}
  atom_predicates = {}
  for k, pred in predicate_by_field_name.items():
    if col := CHAIN_FIELDS.get(k):
      chain_predicates[col] = pred
    elif col := RESIDUE_FIELDS.get(k):
      res_predicates[col] = pred
    elif col := ATOM_FIELDS.get(k):
      atom_predicates[col] = pred
    else:
      raise ValueError(k)
  return chain_predicates, res_predicates, atom_predicates


_T = TypeVar('_T')


SCALAR_FIELDS: Final[Collection[str]] = frozenset({
    'name',
    'release_date',
    'resolution',
    'structure_method',
    'bioassembly_data',
    'chemical_components_data',
})


TABLE_FIELDS: Final[Collection[str]] = frozenset(
    {'chains', 'residues', 'atoms', 'bonds'}
)


V2_FIELDS: Final[Collection[str]] = frozenset({*SCALAR_FIELDS, *TABLE_FIELDS})


@dataclasses.dataclass(frozen=True, slots=True, kw_only=True)
class StructureTables:
  chains: structure_tables.Chains
  residues: structure_tables.Residues
  atoms: structure_tables.Atoms
  bonds: structure_tables.Bonds


@dataclasses.dataclass(frozen=True, slots=True, kw_only=True)
class ResArrays:
  """Atom-level data arrays with a residue dimension.

  Attributes:
    atom_positions: float32 of shape [num_res, num_atom_type, 3] coordinates.
    atom_mask: float32 of shape [num_res, num_atom_type] indicating if an atom
      is present.
    atom_b_factor: float32 of shape [num_res, num_atom_type] b_factors.
    atom_occupancy: float32 of shape [num_res, num_atom_type] occupancies.
  """

  atom_positions: np.ndarray
  atom_mask: np.ndarray
  atom_b_factor: np.ndarray
  atom_occupancy: np.ndarray


class Structure(table.Database):
  """Structure class for representing and processing molecular structures."""

  tables: ClassVar[Collection[str]] = TABLE_FIELDS

  foreign_keys: ClassVar[Mapping[str, Collection[tuple[str, str]]]] = {
      'residues': (('chain_key', 'chains'),),
      'atoms': (('chain_key', 'chains'), ('res_key', 'residues')),
      'bonds': (('from_atom_key', 'atoms'), ('dest_atom_key', 'atoms')),
  }

  def __init__(
      self,
      *,
      name: str = 'unset',
      release_date: datetime.date | None = None,
      resolution: float | None = None,
      structure_method: str | None = None,
      bioassembly_data: bioassemblies.BioassemblyData | None = None,
      chemical_components_data: (
          struc_chem_comps.ChemicalComponentsData | None
      ) = None,
      chains: structure_tables.Chains,
      residues: structure_tables.Residues,
      atoms: structure_tables.Atoms,
      bonds: structure_tables.Bonds,
      skip_validation: bool = False,
  ):
    # Version number is written to mmCIF and should be incremented when changes
    # are made to mmCIF writing or internals that affect this.
    # b/345221494 Rename this variable when structure_v1 compatibility code
    # is removed.
    self._VERSION = '2.0.0'  # pylint: disable=invalid-name
    self._name = name or 'unset'
    self._release_date = release_date
    self._resolution = resolution
    self._structure_method = structure_method
    self._bioassembly_data = bioassembly_data
    self._chemical_components_data = chemical_components_data

    self._chains = chains
    self._residues = residues
    self._atoms = atoms
    self._bonds = bonds

    if not skip_validation:
      self._validate_table_foreign_keys()
      self._validate_consistent_table_ordering()

  def _validate_table_foreign_keys(self):
    """Validates that all foreign keys are present in the referred tables."""
    residue_keys = set(self._residues.key)
    chain_keys = set(self._chains.key)
    if np.any(membership.isin(self._atoms.res_key, residue_keys, invert=True)):
      raise ValueError(
          'Atom residue keys not in the residues table: '
          f'{set(self._atoms.res_key).difference(self._residues.key)}'
      )
    if np.any(membership.isin(self._atoms.chain_key, chain_keys, invert=True)):
      raise ValueError(
          'Atom chain keys not in the chains table: '
          f'{set(self._atoms.chain_key).difference(self._chains.key)}'
      )
    if np.any(
        membership.isin(self._residues.chain_key, chain_keys, invert=True)
    ):
      raise ValueError(
          'Residue chain keys not in the chains table: '
          f'{set(self._residues.chain_key).difference(self._chains.key)}'
      )

  def _validate_consistent_table_ordering(self):
    """Validates that all tables have the same ordering."""
    atom_chain_keys = self._atoms.chain_key[self.chain_boundaries]
    atom_res_keys = self._atoms.res_key[self.res_boundaries]

    if not np.array_equal(self.present_chains.key, atom_chain_keys):
      raise ValueError(
          f'Atom table chain order\n{atom_chain_keys}\ndoes not match the '
          f'chain table order\n{self._chains.key}'
      )
    if not np.array_equal(self.present_residues.key, atom_res_keys):
      raise ValueError(
          f'Atom table residue order\n{atom_res_keys}\ndoes not match the '
          f'present residue table order\n{self.present_residues.key}'
      )

  def get_table(self, table_name: str) -> table.Table:
    match table_name:
      case 'chains':
        return self.chains_table
      case 'residues':
        return self.residues_table
      case 'atoms':
        return self.atoms_table
      case 'bonds':
        return self.bonds_table
      case _:
        raise ValueError(table_name)

  @property
  def chains_table(self) -> structure_tables.Chains:
    """Chains table."""
    return self._chains

  @property
  def residues_table(self) -> structure_tables.Residues:
    """Residues table."""
    return self._residues

  @property
  def atoms_table(self) -> structure_tables.Atoms:
    """Atoms table."""
    return self._atoms

  @property
  def bonds_table(self) -> structure_tables.Bonds:
    """Bonds table."""
    return self._bonds

  @property
  def name(self) -> str:
    return self._name

  @property
  def release_date(self) -> datetime.date | None:
    return self._release_date

  @property
  def resolution(self) -> float | None:
    return self._resolution

  @property
  def structure_method(self) -> str | None:
    return self._structure_method

  @property
  def bioassembly_data(self) -> bioassemblies.BioassemblyData | None:
    return self._bioassembly_data

  @property
  def chemical_components_data(
      self,
  ) -> struc_chem_comps.ChemicalComponentsData | None:
    return self._chemical_components_data

  @property
  def bonds(self) -> structure_tables.Bonds:
    return self._bonds

  @functools.cached_property
  def author_naming_scheme(self) -> AuthorNamingScheme:
    auth_asym_id = {}
    entity_id = {}
    entity_desc = {}
    auth_seq_id = collections.defaultdict(dict)
    insertion_code = collections.defaultdict(dict)

    for chain_i in range(self._chains.size):
      chain_id = self._chains.id[chain_i]
      auth_asym_id[chain_id] = self._chains.auth_asym_id[chain_i]
      chain_entity_id = self._chains.entity_id[chain_i]
      entity_id[chain_id] = chain_entity_id
      entity_desc[chain_entity_id] = self._chains.entity_desc[chain_i]

    chain_index_by_key = self._chains.index_by_key
    for res_i in range(self._residues.size):
      chain_key = self._residues.chain_key[res_i]
      chain_id = self._chains.id[chain_index_by_key[chain_key]]
      res_id = self._residues.id[res_i]
      res_auth_seq_id = self._residues.auth_seq_id[res_i]
      if res_auth_seq_id == MISSING_AUTH_SEQ_ID:
        continue
      auth_seq_id[chain_id][res_id] = res_auth_seq_id
      ins_code = self._residues.insertion_code[res_i]
      # Compatibility with Structure v1 which used None to represent . or ?.
      insertion_code[chain_id][res_id] = (
          ins_code if ins_code not in {'.', '?'} else None
      )

    return AuthorNamingScheme(
        auth_asym_id=auth_asym_id,
        entity_id=entity_id,
        entity_desc=entity_desc,
        auth_seq_id=dict(auth_seq_id),
        insertion_code=dict(insertion_code),
    )

  @functools.cached_property
  def all_residues(self) -> AllResidues:
    chain_id_by_key = dict(zip(self._chains.key, self._chains.id))
    residue_chain_boundaries = _get_change_indices(self._residues.chain_key)
    boundaries = self._iter_residue_ranges(
        residue_chain_boundaries, count_unresolved=True
    )
    return {
        chain_id_by_key[self._residues.chain_key[start]]: list(
            zip(self._residues.name[start:end], self._residues.id[start:end])
        )
        for start, end in boundaries
    }

  @functools.cached_property
  def label_asym_id_to_entity_id(self) -> Mapping[str, str]:
    return dict(zip(self._chains.id, self._chains.entity_id))

  @functools.cached_property
  def chain_entity_id(self) -> np.ndarray:
    """Returns the entity ID for each atom in the structure."""
    return self.chains_table.apply_array_to_column(
        'entity_id', self._atoms.chain_key
    )

  @functools.cached_property
  def chain_entity_desc(self) -> np.ndarray:
    """Returns the entity description for each atom in the structure."""
    return self.chains_table.apply_array_to_column(
        'entity_desc', self._atoms.chain_key
    )

  @functools.cached_property
  def chain_auth_asym_id(self) -> np.ndarray:
    """Returns the chain auth asym ID for each atom in the structure."""
    return self.chains_table.apply_array_to_column(
        'auth_asym_id', self._atoms.chain_key
    )

  @functools.cached_property
  def chain_id(self) -> np.ndarray:
    chain_index_by_key = self._chains.index_by_key
    return self._chains.id[chain_index_by_key[self._atoms.chain_key]]

  @functools.cached_property
  def chain_type(self) -> np.ndarray:
    chain_index_by_key = self._chains.index_by_key
    return self._chains.type[chain_index_by_key[self._atoms.chain_key]]

  @functools.cached_property
  def res_id(self) -> np.ndarray:
    return self._residues['id', self._atoms.res_key]

  @functools.cached_property
  def res_name(self) -> np.ndarray:
    return self._residues['name', self._atoms.res_key]

  @functools.cached_property
  def res_auth_seq_id(self) -> np.ndarray:
    """Returns the residue auth seq ID for each atom in the structure."""
    return self.residues_table.apply_array_to_column(
        'auth_seq_id', self._atoms.res_key
    )

  @functools.cached_property
  def res_insertion_code(self) -> np.ndarray:
    """Returns the residue insertion code for each atom in the structure."""
    return self.residues_table.apply_array_to_column(
        'insertion_code', self._atoms.res_key
    )

  @property
  def atom_key(self) -> np.ndarray:
    return self._atoms.key

  @property
  def atom_name(self) -> np.ndarray:
    return self._atoms.name

  @property
  def atom_element(self) -> np.ndarray:
    return self._atoms.element

  @property
  def atom_x(self) -> np.ndarray:
    return self._atoms.x

  @property
  def atom_y(self) -> np.ndarray:
    return self._atoms.y

  @property
  def atom_z(self) -> np.ndarray:
    return self._atoms.z

  @property
  def atom_b_factor(self) -> np.ndarray:
    return self._atoms.b_factor

  @property
  def atom_occupancy(self) -> np.ndarray:
    return self._atoms.occupancy

  @functools.cached_property
  def chain_boundaries(self) -> np.ndarray:
    """The indices in the atom fields where each present chain begins."""
    return _get_change_indices(self._atoms.chain_key)

  @functools.cached_property
  def res_boundaries(self) -> np.ndarray:
    """The indices in the atom fields where each present residue begins."""
    return _get_change_indices(self._atoms.res_key)

  @functools.cached_property
  def present_chains(self) -> structure_tables.Chains:
    """Returns table of chains which have at least 1 resolved atom."""
    is_present_mask = np.isin(self._chains.key, self._atoms.chain_key)
    return typing.cast(structure_tables.Chains, self._chains[is_present_mask])

  @functools.cached_property
  def present_residues(self) -> structure_tables.Residues:
    """Returns table of residues which have at least 1 resolved atom."""
    is_present_mask = np.isin(self._residues.key, self._atoms.res_key)
    return typing.cast(
        structure_tables.Residues, self._residues[is_present_mask]
    )

  @functools.cached_property
  def unresolved_residues(self) -> structure_tables.Residues:
    """Returns table of residues which have at least 1 resolved atom."""
    is_unresolved_mask = np.isin(
        self._residues.key, self._atoms.res_key, invert=True
    )
    return typing.cast(
        structure_tables.Residues, self._residues[is_unresolved_mask]
    )

  def __getitem__(self, field: str) -> Any:
    """Gets raw field data using field name as a string."""
    if field in TABLE_FIELDS:
      return self.get_table(field)
    else:
      return getattr(self, field)

  def __getstate__(self) -> dict[str, Any]:
    """Pickle calls this on dump.

    Returns:
      Members with cached properties removed.
    """
    cached_props = {
        k
        for k, v in self.__class__.__dict__.items()
        if isinstance(v, functools.cached_property)
    }
    return {k: v for k, v in self.__dict__.items() if k not in cached_props}

  def __repr__(self):
    return (
        f'Structure({self._name}: {self.num_chains} chains, '
        f'{self.num_residues(count_unresolved=False)} residues, '
        f'{self.num_atoms} atoms)'
    )

  @property
  def num_atoms(self) -> int:
    return self._atoms.size

  def num_residues(self, *, count_unresolved: bool) -> int:
    """Returns the number of residues in this Structure.

    Args:
      count_unresolved: Whether to include unresolved (empty) residues.

    Returns:
      Number of residues in the Structure.
    """
    if count_unresolved:
      return self._residues.size
    else:
      return self.present_residues.size

  @property
  def num_chains(self) -> int:
    return self._chains.size

  @property
  def num_models(self) -> int:
    """The number of models of this Structure."""
    return self._atoms.num_models

  def _atom_mask(self, entities: Set[str]) -> np.ndarray:
    """Boolean label indicating if each atom is from entities or not."""
    mask = np.zeros(self.num_atoms, dtype=bool)
    chain_index_by_key = self._chains.index_by_key
    for start, end in self.iter_chain_ranges():
      chain_index = chain_index_by_key[self._atoms.chain_key[start]]
      chain_type = self._chains.type[chain_index]
      mask[start:end] = chain_type in entities
    return mask

  @functools.cached_property
  def is_protein_mask(self) -> np.ndarray:
    """Boolean label indicating if each atom is from protein or not."""
    return self._atom_mask(entities={mmcif_names.PROTEIN_CHAIN})

  @functools.cached_property
  def is_dna_mask(self) -> np.ndarray:
    """Boolean label indicating if each atom is from DNA or not."""
    return self._atom_mask(entities={mmcif_names.DNA_CHAIN})

  @functools.cached_property
  def is_rna_mask(self) -> np.ndarray:
    """Boolean label indicating if each atom is from RNA or not."""
    return self._atom_mask(entities={mmcif_names.RNA_CHAIN})

  @functools.cached_property
  def is_nucleic_mask(self) -> np.ndarray:
    """Boolean label indicating if each atom is a nucleic acid or not."""
    return self._atom_mask(entities=mmcif_names.NUCLEIC_ACID_CHAIN_TYPES)

  @functools.cached_property
  def is_ligand_mask(self) -> np.ndarray:
    """Boolean label indicating if each atom is a ligand or not."""
    return self._atom_mask(entities=mmcif_names.LIGAND_CHAIN_TYPES)

  @functools.cached_property
  def is_water_mask(self) -> np.ndarray:
    """Boolean label indicating if each atom is from water or not."""
    return self._atom_mask(entities={mmcif_names.WATER})

  def iter_atoms(self) -> Iterator[Mapping[str, Any]]:
    """Iterates over the atoms in the structure."""
    if self._atoms.size == 0:
      return

    current_chain = self._chains.get_row_by_key(
        column_name_map=CHAIN_FIELDS, key=self._atoms.chain_key[0]
    )
    current_chain_key = self._atoms.chain_key[0]
    current_res = self._residues.get_row_by_key(
        column_name_map=RESIDUE_FIELDS, key=self._atoms.res_key[0]
    )
    current_res_key = self._atoms.res_key[0]
    for atom_i in range(self._atoms.size):
      atom_chain_key = self._atoms.chain_key[atom_i]
      atom_res_key = self._atoms.res_key[atom_i]

      if atom_chain_key != current_chain_key:
        chain_index = self._chains.index_by_key[atom_chain_key]
        current_chain = {
            'chain_id': self._chains.id[chain_index],
            'chain_type': self._chains.type[chain_index],
            'chain_auth_asym_id': self._chains.auth_asym_id[chain_index],
            'chain_entity_id': self._chains.entity_id[chain_index],
            'chain_entity_desc': self._chains.entity_desc[chain_index],
        }
        current_chain_key = atom_chain_key
      if atom_res_key != current_res_key:
        res_index = self._residues.index_by_key[atom_res_key]
        current_res = {
            'res_id': self._residues.id[res_index],
            'res_name': self._residues.name[res_index],
            'res_auth_seq_id': self._residues.auth_seq_id[res_index],
            'res_insertion_code': self._residues.insertion_code[res_index],
        }
        current_res_key = atom_res_key

      yield {
          'atom_name': self._atoms.name[atom_i],
          'atom_element': self._atoms.element[atom_i],
          'atom_x': self._atoms.x[..., atom_i],
          'atom_y': self._atoms.y[..., atom_i],
          'atom_z': self._atoms.z[..., atom_i],
          'atom_b_factor': self._atoms.b_factor[..., atom_i],
          'atom_occupancy': self._atoms.occupancy[..., atom_i],
          'atom_key': self._atoms.key[atom_i],
          **current_res,
          **current_chain,
      }

  def iter_residues(
      self,
      include_unresolved: bool = False,
  ) -> Iterator[Mapping[str, Any]]:
    """Iterates over the residues in the structure."""
    res_table = self._residues if include_unresolved else self.present_residues
    if res_table.size == 0:
      return

    current_chain = self._chains.get_row_by_key(
        column_name_map=CHAIN_FIELDS, key=res_table.chain_key[0]
    )
    current_chain_key = res_table.chain_key[0]
    for res_i in range(res_table.size):
      res_chain_key = res_table.chain_key[res_i]

      if res_chain_key != current_chain_key:
        current_chain = self._chains.get_row_by_key(
            column_name_map=CHAIN_FIELDS, key=res_table.chain_key[res_i]
        )
        current_chain_key = res_chain_key

      row = {
          'res_id': res_table.id[res_i],
          'res_name': res_table.name[res_i],
          'res_auth_seq_id': res_table.auth_seq_id[res_i],
          'res_insertion_code': res_table.insertion_code[res_i],
      }
      yield row | current_chain

  def _iter_atom_ranges(
      self, boundaries: Sequence[int] | np.ndarray
  ) -> Iterator[tuple[int, int]]:
    """Iterator for (start, end) pairs from an array of start indices."""
    yield from itertools.pairwise(boundaries)
    # Use explicit length test as boundaries can be a NumPy array.
    if len(boundaries) > 0:  # pylint: disable=g-explicit-length-test
      yield boundaries[-1], self.num_atoms

  def _iter_residue_ranges(
      self,
      boundaries: Sequence[int] | np.ndarray,
      *,
      count_unresolved: bool,
  ) -> Iterator[tuple[int, int]]:
    """Iterator for (start, end) pairs from an array of start indices."""
    yield from itertools.pairwise(boundaries)
    # Use explicit length test as boundaries can be a NumPy array.
    if len(boundaries) > 0:  # pylint: disable=g-explicit-length-test
      yield boundaries[-1], self.num_residues(count_unresolved=count_unresolved)

  def iter_chain_ranges(self) -> Iterator[tuple[int, int]]:
    """Iterates pairs of (chain_start, chain_end) indices.

    Yields:
      Pairs of (start, end) indices for each chain, where end is not inclusive.
      i.e. struc.chain_id[start:end] would be a constant array with length
      equal to the number of atoms in the chain.
    """
    yield from self._iter_atom_ranges(self.chain_boundaries)

  def iter_residue_ranges(self) -> Iterator[tuple[int, int]]:
    """Iterates pairs of (residue_start, residue_end) indices.

    Yields:
      Pairs of (start, end) indices for each residue, where end is not
      inclusive. i.e. struc.res_id[start:end] would be a constant array with
      length equal to the number of atoms in the residue.
    """
    yield from self._iter_atom_ranges(self.res_boundaries)

  def iter_chains(self) -> Iterator[Mapping[str, Any]]:
    """Iterates over the chains in the structure."""
    for chain_i in range(self.present_chains.size):
      yield {
          'chain_id': self.present_chains.id[chain_i],
          'chain_type': self.present_chains.type[chain_i],
          'chain_auth_asym_id': self.present_chains.auth_asym_id[chain_i],
          'chain_entity_id': self.present_chains.entity_id[chain_i],
          'chain_entity_desc': self.present_chains.entity_desc[chain_i],
      }

  def iter_bonds(self) -> Iterator[Bond]:
    """Iterates over the atoms and bond information.

    Example usage:

    ```
    for from_atom, dest_atom, bond_info in struc.iter_bonds():
      print(
          f'From atom: name={from_atom["atom_name"]}, '
          f'chain={from_atom["chain_id"]}, ...'
      )
      # Same for dest_atom
      print(f'Bond info: type={bond_info["type"]}, role={bond_info["role"]}')
    ```

    Yields:
      A `Bond` NamedTuple for each bond in the bonds table.
      These have fields `from_atom`, `dest_atom`, `bond_info` where each
      is a dictionary. The first two have the same keys as the atom dicts
      returned by self.iter_atoms() -- i.e. one key per non-None field.
      The final dict has the same keys as self.bonds.iterrows() -- i.e. one
      key per column in the bonds table.
    """
    from_atom_iter = self._atoms.iterrows(
        row_keys=self._bonds.from_atom_key,
        column_name_map=ATOM_FIELDS,
        chain_key=self._chains.with_column_names(CHAIN_FIELDS),
        res_key=self._residues.with_column_names(RESIDUE_FIELDS),
    )
    dest_atom_iter = self._atoms.iterrows(
        row_keys=self._bonds.dest_atom_key,
        column_name_map=ATOM_FIELDS,
        chain_key=self._chains.with_column_names(CHAIN_FIELDS),
        res_key=self._residues.with_column_names(RESIDUE_FIELDS),
    )

    for from_atom, dest_atom, bond_info in zip(
        from_atom_iter, dest_atom_iter, self._bonds.iterrows(), strict=True
    ):
      yield Bond(from_atom=from_atom, dest_atom=dest_atom, bond_info=bond_info)

  def _apply_atom_index_array(
      self,
      index_arr: np.ndarray,
      chain_boundaries: np.ndarray | None = None,
      res_boundaries: np.ndarray | None = None,
      skip_validation: bool = False,
  ) -> Self:
    """Applies index_arr to the atom table using NumPy-style array indexing.

    Args:
      index_arr: A 1D NumPy array that will be used to index into the atoms
        table. This can either be a boolean array to act as a mask, or an
        integer array to perform a gather operation.
      chain_boundaries: Unused in structure v2.
      res_boundaries: Unused in structure v2.
      skip_validation: Whether to skip the validation step that checks internal
        consistency after applying atom index array. Do not set to True unless
        you are certain the transform is safe, e.g. when the order of atoms is
        guaranteed to not change.

    Returns:
      A new Structure with an updated atoms table.
    """
    del chain_boundaries, res_boundaries

    if index_arr.ndim != 1:
      raise ValueError(
          f'index_arr must be a 1D NumPy array, but has shape {index_arr.shape}'
      )

    if index_arr.dtype == bool and np.all(index_arr):
      return self  # Shortcut: The operation is a no-op, so just return itself.

    atoms = structure_tables.Atoms(
        **{col: self._atoms[col][..., index_arr] for col in self._atoms.columns}
    )
    updated_tables = self._cascade_delete(atoms=atoms)
    return self.copy_and_update(
        atoms=updated_tables.atoms,
        bonds=updated_tables.bonds,
        skip_validation=skip_validation,
    )

  @property
  def group_by_residue(self) -> Self:
    """Returns a Structure with one atom per residue.

    e.g. restypes = struc.group_by_residue['res_id']

    Returns:
      A new Structure with one atom per residue such that per-atom arrays
      such as res_name (i.e. Structure v1 fields) have one element per residue.
    """
    # This use of _apply_atom_index_array is safe because the chain/residue/atom
    # ordering won't change (essentially applying a residue start mask).
    return self._apply_atom_index_array(
        self.res_boundaries, skip_validation=True
    )

  @property
  def group_by_chain(self) -> Self:
    """Returns a Structure where all fields are per-chain.

    e.g. chains = struc.group_by_chain['chain_id']

    Returns:
      A new Structure with one atom per chain such that per-atom arrays
      such as res_name (i.e. Structure v1 fields) have one element per chain.
    """
    # This use of _apply_atom_index_array is safe because the chain/residue/atom
    # ordering won't change (essentially applying a chain start mask).
    return self._apply_atom_index_array(
        self.chain_boundaries, skip_validation=True
    )

  @property
  def with_sorted_chains(self) -> Self:
    """Returns a new structure with the chains are in reverse spreadsheet style.

    This is the usual order to write chains in an mmCIF:
    (A < B < ... < AA < BA < CA < ... < AB < BB < CB ...)

    NB: this method will fail if chains do not conform to this mmCIF naming
    convention.

    Only to be used for third party metrics that rely on the chain order.
    Elsewhere chains should be identified by name and code should be agnostic to
    the order.
    """
    sorted_chains = sorted(self.chains, key=mmcif.str_id_to_int_id)
    return self.reorder_chains(new_order=sorted_chains)

  @functools.cached_property
  def atom_ids(self) -> Sequence[tuple[str, str, None, str]]:
    """Gets a list of atom ID tuples from Structure class arrays.

    Returns:
      A list of tuples of (chain_id, res_id, insertion_code, atom_name) where
      insertion code is always None. There is one element per atom, and the
      list is ordered according to the order of atoms in the input arrays.
    """
    # Convert to Numpy strings, then to Python strings (dtype=object).
    res_ids = self.residues_table.id.astype(str).astype(object)
    res_ids = res_ids[
        self.residues_table.index_by_key[self.atoms_table.res_key]
    ]
    ins_codes = [None] * self.num_atoms
    return list(
        zip(self.chain_id, res_ids, ins_codes, self.atom_name, strict=True)
    )

  def order_and_drop_atoms_to_match(
      self,
      other: 'Structure',
      *,
      allow_missing_atoms: bool = False,
  ) -> Self:
    """Returns a new structure with atoms ordered & dropped to match another's.

    This performs two operations simultaneously:
      * Ordering the atoms in this structure to match the order in the other.
      * Dropping atoms in this structure that do not appear in the other.

    Example:
    Consider a prediction and ground truth with the following atoms, described
    using tuples of `(chain_id, res_id, atom_name)`:
      * `prediction:   [(A, 1, CA), (A, 1, N), (A, 2, CA), (B, 1, CA)]`
      * `ground_truth: [(B, 1, CA), (A, 1, N), (A, 1, CA)]`
    Note how the ground truth is missing the `(A, 2, CA)` atom and also
    has the atoms in a different order. This method returns a modified
    prediction that has reordered atoms and without any atoms not in the ground
    truth so that its atom list looks the same as the ground truth atom list.
    This means `prediction.coords` and `ground_truth.coords` now have the
    same shape and can be compared across the atom dimension.

    Note that matching residues with no atoms and matching chains with no
    residues will also be kept. E.g. in the example above, if prediction and
    ground truth both had an unresolved residue (A, 3), the output structure
    will also have an unresolved residue (A, 3).

    Args:
      other: Another `Structure`. This provides the reference ordering that is
        used to sort this structure's atom arrays.
      allow_missing_atoms: Whether to skip atoms present in `other` but not this
        structure and return a structure containing a subset of the atoms in the
        other structure.

    Returns:
      A new `Structure`, based on this structure, which, if
      `allow_missing_atoms` is False, contains exactly the same atoms as in
      the `other` structure and which matches the `other` structure in terms
      of the order of the atoms in the field arrays. Otherwise, if missing
      atoms are allowed then the resulting structure contains a subset of
      those atoms in the other structure.

    Raises:
      MissingAtomError: If there are atoms present in the other structure that
        cannot be found in this structure.
    """
    atom_index_map = {atom_id: i for i, atom_id in enumerate(self.atom_ids)}
    try:
      if allow_missing_atoms:
        # Only include atoms that were found in the other structure.
        atom_indices = [
            atom_index
            for atom_id in other.atom_ids
            if (atom_index := atom_index_map.get(atom_id)) is not None
        ]
      else:
        atom_indices = [
            atom_index_map[atom_id]  # Hard fail on missing.
            for atom_id in other.atom_ids
        ]
    except KeyError as e:
      if len(e.args[0]) == 4:
        chain_id, res_id, ins_code, atom_name = e.args[0]
        raise MissingAtomError(
            f'No atom in this structure (name: {self._name}) matches atom in '
            f'other structure (name: {other.name}) with internal (label) chain '
            f'ID {chain_id}, residue ID {res_id}, insertion code {ins_code} '
            f'and atom name {atom_name}.'
        ) from e
      else:
        raise

    def _iter_residues(struc: Self) -> Iterable[tuple[str, str]]:
      yield from zip(
          struc.chains_table['id', struc.residues_table.chain_key],
          struc.residues_table.id,
          strict=True,
      )

    chain_index_map = {
        chain_id: i for i, chain_id in enumerate(self._chains.id)
    }
    chain_indices = [
        chain_index
        for chain_id in other.chains_table.id
        if (chain_index := chain_index_map.get(chain_id)) is not None
    ]
    residue_index_map = {
        res_id: i for i, res_id in enumerate(_iter_residues(self))
    }
    res_indices = [
        residue_index
        for res_id in _iter_residues(other)
        if (residue_index := residue_index_map.get(res_id)) is not None
    ]

    # Reorder all tables.
    chains = self._chains.apply_index(np.array(chain_indices, dtype=np.int64))
    residues = self._residues.apply_index(np.array(res_indices, dtype=np.int64))
    atoms = self._atoms.apply_index(np.array(atom_indices, dtype=np.int64))

    # Get chain keys in the order they appear in the atoms table.
    new_chain_boundaries = _get_change_indices(atoms.chain_key)
    new_chain_key_order = atoms.chain_key[new_chain_boundaries]
    if len(new_chain_key_order) != len(set(new_chain_key_order)):
      raise ValueError(
          f'Chain keys not contiguous after reordering: {new_chain_key_order}'
      )

    # Get residue keys in the order they appear in the atoms table.
    new_res_boundaries = _get_change_indices(atoms.res_key)
    new_res_key_order = atoms.res_key[new_res_boundaries]
    if len(new_res_key_order) != len(set(new_res_key_order)):
      raise ValueError(
          f'Residue keys not contiguous after reordering: {new_res_key_order}'
      )

    # If any atoms were deleted, propagate that into the bonds table.
    updated_tables = self._cascade_delete(
        chains=chains,
        residues=residues,
        atoms=atoms,
    )
    return self.copy_and_update(
        chains=chains,
        residues=residues,
        atoms=updated_tables.atoms,
        bonds=updated_tables.bonds,
    )

  def copy_and_update(
      self,
      *,
      name: str | _UnsetType = _UNSET,
      release_date: datetime.date | None | _UnsetType = _UNSET,
      resolution: float | None | _UnsetType = _UNSET,
      structure_method: str | None | _UnsetType = _UNSET,
      bioassembly_data: (
          bioassemblies.BioassemblyData | None | _UnsetType
      ) = _UNSET,
      chemical_components_data: (
          struc_chem_comps.ChemicalComponentsData | None | _UnsetType
      ) = _UNSET,
      chains: structure_tables.Chains | None | _UnsetType = _UNSET,
      residues: structure_tables.Residues | None | _UnsetType = _UNSET,
      atoms: structure_tables.Atoms | None | _UnsetType = _UNSET,
      bonds: structure_tables.Bonds | None | _UnsetType = _UNSET,
      skip_validation: bool = False,
  ) -> Self:
    """Performs a shallow copy but with specified fields updated."""

    def all_unset(fields):
      return all(field == _UNSET for field in fields)

    if all_unset((chains, residues, atoms, bonds)):
      if all_unset((
          name,
          release_date,
          resolution,
          structure_method,
          bioassembly_data,
          chemical_components_data,
      )):
        raise ValueError(
            'Unnecessary call to copy_and_update with no changes. As Structure'
            ' and its component tables are immutable, there is no need to copy'
            ' it. Any subsequent operation that modifies structure will return'
            ' a new object.'
        )
      else:
        raise ValueError(
            'When only changing global fields, prefer to use the specialised '
            'copy_and_update_globals.'
        )

    def select(field, default):
      return field if field != _UNSET else default

    return Structure(
        name=select(name, self.name),
        release_date=select(release_date, self.release_date),
        resolution=select(resolution, self.resolution),
        structure_method=select(structure_method, self.structure_method),
        bioassembly_data=select(bioassembly_data, self.bioassembly_data),
        chemical_components_data=select(
            chemical_components_data, self.chemical_components_data
        ),
        chains=select(chains, self._chains),
        residues=select(residues, self._residues),
        atoms=select(atoms, self._atoms),
        bonds=select(bonds, self._bonds),
        skip_validation=skip_validation,
    )

  def _copy_and_update(
      self, skip_validation: bool = False, **changes: Any
  ) -> Self:
    """Performs a shallow copy but with specified fields updated."""
    if not changes:
      raise ValueError(
          'Unnecessary call to copy_and_update with no changes. As Structure '
          'and its component tables are immutable, there is no need to copy '
          'it. Any subsequent operation that modifies structure will return a '
          'new object.'
      )

    if 'author_naming_scheme' in changes:
      raise ValueError(
          'Updating using author_naming_scheme is not supported. Update '
          'auth_asym_id, entity_id, entity_desc fields directly in the chains '
          'table and auth_seq_id, insertion_code in the residues table.'
      )

    if all(k in GLOBAL_FIELDS for k in changes):
      raise ValueError(
          'When only changing global fields, prefer to use the specialised '
          'copy_and_update_globals.'
      )

    if all(k in V2_FIELDS for k in changes):
      constructor_kwargs = {field: self[field] for field in V2_FIELDS}
      constructor_kwargs.update(changes)
    elif any(k in ('atoms', 'residues', 'chains') for k in changes):
      raise ValueError(
          'Cannot specify atoms/chains/residues table changes with non-v2'
          f' constructor params: {changes.keys()}'
      )
    elif all(k in ATOM_FIELDS for k in changes):
      if 'atom_key' not in changes:
        raise ValueError(
            'When only changing atom fields, prefer to use the specialised '
            'copy_and_update_atoms.'
        )
      # Only atom fields are being updated, do that directly on the atoms table.
      updated_atoms = self._atoms.copy_and_update(
          **{ATOM_FIELDS[k]: v for k, v in changes.items()}
      )
      constructor_kwargs = {
          field: self[field] for field in V2_FIELDS if field != 'atoms'
      }
      constructor_kwargs['atoms'] = updated_atoms
    else:
      constructor_kwargs = {field: self[field] for field in _UPDATEABLE_FIELDS}
      constructor_kwargs.update(changes)
    return Structure(skip_validation=skip_validation, **constructor_kwargs)

  def copy_and_update_coords(self, coords: np.ndarray) -> Self:
    """Performs a shallow copy but with coordinates updated."""
    if coords.shape[-2:] != (self.num_atoms, 3):
      raise ValueError(
          f'{coords.shape=} does not have last dimensions ({self.num_atoms}, 3)'
      )
    updated_atoms = self._atoms.copy_and_update_coords(coords)
    return self.copy_and_update(atoms=updated_atoms, skip_validation=True)

  def copy_and_update_from_res_arrays(
      self,
      *,
      include_unresolved: bool = False,
      **changes: np.ndarray,
  ) -> Self:
    """Like copy_and_update but changes are arrays of length num_residues.

    These changes are first scattered into arrays of length num_atoms such
    that each value is repeated across the residue at that index, then they
    are used as the new values of these fields.

    E.g.
    * This structure's res_id: 1, 1, 1, 2, 3, 3 (3 res, 6 atoms)
    * new atom_b_factor: 7, 8, 9
    * Returned structure's atom_b_factor: 7, 7, 7, 8, 9, 9

    Args:
      include_unresolved: Whether the provided list of new values per residue
        include values for all residues, or only those that are resolved.
      **changes: kwargs corresponding to atom array fields, e.g. atom_x or
        atom_b_factor, but with length num_residues rather than num_atoms. Note
        that changing atom_key this way is is not supported.

    Returns:
      A new `Structure` with all fields other than those specified as kwargs
      shallow copied from this structure. The values of the kwargs are
      scattered across the atom arrays and then used to overwrite these
      fields for the returned structure.
    """
    if not all(c in set(ATOM_FIELDS) - {'atom_key'} for c in changes):
      raise ValueError(
          'Changes must only be to atom fields, got changes to'
          f' {changes.keys()}'
      )

    num_residues = self.num_residues(count_unresolved=include_unresolved)

    for field_name, new_values in changes.items():
      if len(new_values) != num_residues:
        raise ValueError(
            f'{field_name} array of length {len(new_values)} does not match '
            f'{num_residues=} - is include_unresolved set correctly?'
        )

    # We cannot assume that atom_table.res_keys are the relevant indices of the
    # residue table.

    # Therefore we need to construct a map from res_key to the new values and
    # update the atoms_table with that.
    if include_unresolved:
      target_keys = self.residues_table.key
    else:
      target_keys = self.present_residues.key

    new_atom_columns = {}
    for field_name, new_values in changes.items():
      value_by_key = dict(zip(target_keys, new_values, strict=True))
      # pylint: disable=cell-var-from-loop
      new_atom_columns[field_name] = np.vectorize(lambda x: value_by_key[x])(
          self.atoms_table.res_key
      )
      # pylint: enable=cell-var-from-loop

    return self.copy_and_update_atoms(**new_atom_columns)

  def copy_and_update_globals(
      self,
      *,
      name: str | _UnsetType = _UNSET,
      release_date: datetime.date | _UnsetType | None = _UNSET,
      resolution: float | _UnsetType | None = _UNSET,
      structure_method: str | _UnsetType | None = _UNSET,
      bioassembly_data: (
          bioassemblies.BioassemblyData | _UnsetType | None
      ) = _UNSET,
      chemical_components_data: (
          struc_chem_comps.ChemicalComponentsData | _UnsetType | None
      ) = _UNSET,
  ) -> Self:
    """Returns a shallow copy with the global columns updated."""

    def select(field, default):
      return field if field != _UNSET else default

    name = select(name, self.name)
    release_date = select(release_date, self.release_date)
    resolution = select(resolution, self.resolution)
    structure_method = select(structure_method, self.structure_method)
    bioassembly_data = select(bioassembly_data, self.bioassembly_data)
    chem_data = select(chemical_components_data, self.chemical_components_data)

    return Structure(
        name=name,
        release_date=release_date,
        resolution=resolution,
        structure_method=structure_method,
        bioassembly_data=bioassembly_data,
        chemical_components_data=chem_data,
        atoms=self._atoms,
        residues=self._residues,
        chains=self._chains,
        bonds=self._bonds,
    )

  def copy_and_update_atoms(
      self,
      *,
      atom_name: np.ndarray | None = None,
      atom_element: np.ndarray | None = None,
      atom_x: np.ndarray | None = None,
      atom_y: np.ndarray | None = None,
      atom_z: np.ndarray | None = None,
      atom_b_factor: np.ndarray | None = None,
      atom_occupancy: np.ndarray | None = None,
  ) -> Self:
    """Returns a shallow copy with the atoms table updated."""
    new_atoms = structure_tables.Atoms(
        key=self._atoms.key,
        res_key=self._atoms.res_key,
        chain_key=self._atoms.chain_key,
        name=atom_name if atom_name is not None else self.atom_name,
        element=atom_element if atom_element is not None else self.atom_element,
        x=atom_x if atom_x is not None else self.atom_x,
        y=atom_y if atom_y is not None else self.atom_y,
        z=atom_z if atom_z is not None else self.atom_z,
        b_factor=(
            atom_b_factor if atom_b_factor is not None else self.atom_b_factor
        ),
        occupancy=(
            atom_occupancy
            if atom_occupancy is not None
            else self.atom_occupancy
        ),
    )
    return self.copy_and_update(atoms=new_atoms)

  def copy_and_update_residues(
      self,
      *,
      res_id: np.ndarray | None = None,
      res_name: np.ndarray | None = None,
      res_auth_seq_id: np.ndarray | None = None,
      res_insertion_code: np.ndarray | None = None,
  ) -> Self:
    """Returns a shallow copy with the residues table updated."""
    new_residues = structure_tables.Residues(
        key=self._residues.key,
        chain_key=self._residues.chain_key,
        id=res_id if res_id is not None else self._residues.id,
        name=res_name if res_name is not None else self._residues.name,
        auth_seq_id=res_auth_seq_id
        if res_auth_seq_id is not None
        else self._residues.auth_seq_id,
        insertion_code=res_insertion_code
        if res_insertion_code is not None
        else self._residues.insertion_code,
    )
    return self.copy_and_update(residues=new_residues)

  def _cascade_delete(
      self,
      *,
      chains: structure_tables.Chains | None = None,
      residues: structure_tables.Residues | None = None,
      atoms: structure_tables.Atoms | None = None,
      bonds: structure_tables.Bonds | None = None,
  ) -> StructureTables:
    """Performs a cascade delete operation on the structure's tables.

    Cascade delete ensures all the tables are consistent after any table fields
    are being updated by cascading any deletions down the hierarchy of tables:
    chains > residues > atoms > bonds.

    E.g.: if a row from residues table is removed then all the atoms in that
    residue will also be removed from the atoms table. In turn this cascades
    also to the bond table, by removing any bond row which involves any of those
    removed atoms. However the chains table will not be modified, even if
    that was the only residue in its chain, because the chains table is above
    the residues table in the hierarchy.

    Args:
      chains: An optional new chains table.
      residues: An optional new residues table.
      atoms: An optional new atoms table.
      bonds: An optional new bonds table.

    Returns:
      A StructureTables object with the updated tables.
    """
    if chains_unchanged := chains is None:
      chains = self._chains
    if residues_unchanged := residues is None:
      residues = self._residues
    if atoms_unchanged := atoms is None:
      atoms = self._atoms
    if bonds is None:
      bonds = self._bonds

    if not chains_unchanged:
      residues_mask = membership.isin(residues.chain_key, set(chains.key))  # pylint:disable=attribute-error
      if not np.all(residues_mask):  # Only apply if this is not a no-op.
        residues = residues[residues_mask]
        residues_unchanged = False
    if not residues_unchanged:
      atoms_mask = membership.isin(atoms.res_key, set(residues.key))  # pylint:disable=attribute-error
      if not np.all(atoms_mask):  # Only apply if this is not a no-op.
        atoms = atoms[atoms_mask]
        atoms_unchanged = False
    if not atoms_unchanged:
      bonds = bonds.restrict_to_atoms(atoms.key)
    return StructureTables(
        chains=chains, residues=residues, atoms=atoms, bonds=bonds
    )

  def filter(
      self,
      mask: np.ndarray | None = None,
      *,
      apply_per_element: bool = False,
      invert: bool = False,
      cascade_delete: CascadeDelete = CascadeDelete.CHAINS,
      **predicate_by_field_name: table.FilterPredicate,
  ) -> Self:
    """Filters the structure by field values and returns a new structure.

    Predicates are specified as keyword arguments, with names following the
    pattern: <table_name>_<col_name>, where table_name := (chain|res|atom).
    For instance the auth_seq_id column in the residues table can be filtered
    by passing `res_auth_seq_id=pred_value`. The full list of valid options
    are defined in the `col_by_field_name` fields on the different Table
    dataclasses.

    Predicate values can be either:
      1. A constant value, e.g. 'CA'. In this case then only rows that match
        this value for the given field are retained.
      2. A (non-string) iterable e.g. ('A', 'B'). In this
        case then rows are retained if they match any of the provided values for
        the given field.
      3. A boolean function e.g. lambda b_fac: b_fac < 100.0.
        In this case then only rows that evaluate to True are retained. By
        default this function's parameter is expected to be an array, unless
        apply_per_element=True.

    Example usage:
      # Filter to backbone atoms in residues up to 100 in chain B.
      filtered_struc = struc.filter(
          chain_id='B',
          atom_name=('N', 'CA', 'C'),
          res_id=lambda res_id: res_id < 100)

    Example usage where predicate must be applied per-element:
      # Filter to residues with IDs in either [1, 100) or [300, 400).
      ranges = ((1, 100), (300, 400))
      filtered_struc = struc.filter(
          res_id=lambda i: np.any([start <= i < end for start, end in ranges]),
          apply_per_element=True)

    Example usage of providing a raw mask:
      filtered_struc = struc.filter(struc.atom_b_factor < 10.0)

    Args:
      mask: An optional boolean NumPy array with length equal to num_atoms. If
        provided then this will be combined with the other predicates so that an
        atom is included if it is masked-in *and* matches all the predicates.
      apply_per_element: Whether apply predicates to each element individually,
        or to pass the whole column array to the predicate.
      invert: Whether to remove, rather than retain, the entities which match
        the specified predicates.
      cascade_delete: Whether to remove residues and chains which are left
        unresolved in a cascade. filter operates on the atoms table, removing
        atoms which match the predicate. If all atoms in a residue are removed,
        the residue is "unresolved". The value of this argument then determines
        whether such residues and their parent chains should be deleted. FULL
        implies that all unresolved residues should be deleted, and any chains
        which are left with no resolved residues should be deleted. CHAINS is
        the default behaviour - only chains with no resolved residues, and their
        child residues are deleted. Unresolved residues in partially resolved
        chains remain. NONE implies that no unresolved residues or chains should
        be deleted.
      **predicate_by_field_name: A mapping from field name to a predicate.
        Filtered columns must be 1D arrays. If multiple fields are provided as
        keyword arguments then each predicate is applied and the results are
        combined using a boolean AND operation, so an atom is only retained if
        it passes all predicates.

    Returns:
      A new structure representing a filtered version of the current structure.

    Raises:
      ValueError: If mask is provided and is not a bool array with shape
        (num_atoms,).
    """
    chain_predicates, res_predicates, atom_predicates = (
        _unpack_filter_predicates(predicate_by_field_name)
    )
    # Get boolean masks for each table. These are None if none of the filter
    # parameters affect the table in question.
    chain_mask = self._chains.make_filter_mask(
        **chain_predicates, apply_per_element=apply_per_element
    )
    res_mask = self._residues.make_filter_mask(
        **res_predicates, apply_per_element=apply_per_element
    )
    atom_mask = self._atoms.make_filter_mask(
        mask, **atom_predicates, apply_per_element=apply_per_element
    )
    if atom_mask is None:
      atom_mask = np.ones((self._atoms.size,), dtype=bool)

    # Remove atoms that belong to filtered out chains.
    if chain_mask is not None:
      atom_chain_mask = membership.isin(
          self._atoms.chain_key, set(self._chains.key[chain_mask])
      )
      np.logical_and(atom_mask, atom_chain_mask, out=atom_mask)

    # Remove atoms that belong to filtered out residues.
    if res_mask is not None:
      atom_res_mask = membership.isin(
          self._atoms.res_key, set(self._residues.key[res_mask])
      )
      np.logical_and(atom_mask, atom_res_mask, out=atom_mask)

    final_atom_mask = ~atom_mask if invert else atom_mask

    if cascade_delete == CascadeDelete.NONE and np.all(final_atom_mask):
      return self  # Shortcut: The filter is a no-op, so just return itself.

    filtered_atoms = typing.cast(
        structure_tables.Atoms, self._atoms[final_atom_mask]
    )

    match cascade_delete:
      case CascadeDelete.FULL:
        nonempty_residues_mask = np.isin(
            self._residues.key, filtered_atoms.res_key
        )
        filtered_residues = self._residues[nonempty_residues_mask]
        nonempty_chain_mask = np.isin(
            self._chains.key, filtered_atoms.chain_key
        )
        filtered_chains = self._chains[nonempty_chain_mask]
        updated_tables = self._cascade_delete(
            chains=filtered_chains,
            residues=filtered_residues,
            atoms=filtered_atoms,
        )
      case CascadeDelete.CHAINS:
        # To match v1 behavior we remove chains that have no atoms remaining,
        # and we remove residues in those chains.
        # NB we do not remove empty residues.
        nonempty_chain_mask = membership.isin(
            self._chains.key, set(filtered_atoms.chain_key)
        )
        filtered_chains = self._chains[nonempty_chain_mask]
        updated_tables = self._cascade_delete(
            chains=filtered_chains, atoms=filtered_atoms
        )
      case CascadeDelete.NONE:
        updated_tables = self._cascade_delete(atoms=filtered_atoms)
      case _:
        raise ValueError(f'Unknown cascade_delete behaviour: {cascade_delete}')
    return self.copy_and_update(
        chains=updated_tables.chains,
        residues=updated_tables.residues,
        atoms=updated_tables.atoms,
        bonds=updated_tables.bonds,
        skip_validation=True,
    )

  def filter_out(self, *args, **kwargs) -> Self:
    """Returns a new structure with the specified elements removed."""
    return self.filter(*args, invert=True, **kwargs)

  def filter_to_entity_type(
      self,
      *,
      protein: bool = False,
      rna: bool = False,
      dna: bool = False,
      dna_rna_hybrid: bool = False,
      ligand: bool = False,
      water: bool = False,
  ) -> Self:
    """Filters the structure to only include the selected entity types.

    This convenience method abstracts away the specifics of mmCIF entity
    type names which, especially for ligands, are non-trivial.

    Args:
      protein: Whether to include protein (polypeptide(L)) chains.
      rna: Whether to include RNA chains.
      dna: Whether to include DNA chains.
      dna_rna_hybrid: Whether to include DNA RNA hybrid chains.
      ligand: Whether to include ligand (i.e. not polymer) chains.
      water: Whether to include water chains.

    Returns:
      The filtered structure.
    """
    include_types = []
    if protein:
      include_types.append(mmcif_names.PROTEIN_CHAIN)
    if rna:
      include_types.append(mmcif_names.RNA_CHAIN)
    if dna:
      include_types.append(mmcif_names.DNA_CHAIN)
    if dna_rna_hybrid:
      include_types.append(mmcif_names.DNA_RNA_HYBRID_CHAIN)
    if ligand:
      include_types.extend(mmcif_names.LIGAND_CHAIN_TYPES)
    if water:
      include_types.append(mmcif_names.WATER)
    return self.filter(chain_type=include_types)

  def get_stoichiometry(
      self, *, fix_non_standard_polymer_res: bool = False
  ) -> Sequence[int]:
    """Returns the structure's stoichiometry using chain_res_name_sequence.

    Note that everything is considered (protein, RNA, DNA, ligands) except for
    water molecules. If you are interested only in a certain type of entities,
    filter them out before calling this method.

    Args:
      fix_non_standard_polymer_res: If True, maps non standard residues in
        protein / RNA / DNA chains to standard residues (e.g. MSE -> MET) or UNK
        / N if a match is not found.

    Returns:
      A list of integers, one for each unique chain in the structure,
      determining the number of that chain appearing in the structure. The
      numbers are sorted highest to lowest. E.g. for an A3B2 protein this method
      will return [3, 2].
    """
    filtered = self.filter_to_entity_type(
        protein=True,
        rna=True,
        dna=True,
        dna_rna_hybrid=True,
        ligand=True,
        water=False,
    )
    seqs = filtered.chain_res_name_sequence(
        include_missing_residues=True,
        fix_non_standard_polymer_res=fix_non_standard_polymer_res,
    )

    unique_seq_counts = collections.Counter(seqs.values())
    return sorted(unique_seq_counts.values(), reverse=True)

  def without_hydrogen(self) -> Self:
    """Returns the structure without hydrogen atoms."""
    return self.filter(
        np.logical_and(self._atoms.element != 'H', self._atoms.element != 'D')
    )

  def without_terminal_oxygens(self) -> Self:
    """Returns the structure without terminal oxygen atoms."""
    terminal_oxygen_filter = np.zeros(self.num_atoms, dtype=bool)
    for chain_type, atom_name in mmcif_names.TERMINAL_OXYGENS.items():
      chain_keys = self._chains.key[self._chains.type == chain_type]
      chain_atom_filter = np.logical_and(
          self._atoms.name == atom_name,
          np.isin(self._atoms.chain_key, chain_keys),
      )
      np.logical_or(
          terminal_oxygen_filter, chain_atom_filter, out=terminal_oxygen_filter
      )
    return self.filter_out(terminal_oxygen_filter)

  def reset_author_naming_scheme(self) -> Self:
    """Remove author chain/residue ids, entity info and use internal ids."""
    new_chains = structure_tables.Chains(
        key=self._chains.key,
        id=self._chains.id,
        type=self._chains.type,
        auth_asym_id=self._chains.id,
        entity_id=np.arange(1, self.num_chains + 1).astype(str).astype(object),
        entity_desc=np.full(self.num_chains, '.', dtype=object),
    )
    new_residues = structure_tables.Residues(
        key=self._residues.key,
        chain_key=self._residues.chain_key,
        id=self._residues.id,
        name=self._residues.name,
        auth_seq_id=self._residues.id.astype(str).astype(object),
        insertion_code=np.full(
            self.num_residues(count_unresolved=True), '?', dtype=object
        ),
    )
    return self.copy_and_update(
        chains=new_chains, residues=new_residues, skip_validation=True
    )

  def filter_residues(self, res_mask: np.ndarray) -> Self:
    """Filter resolved residues using a boolean mask."""
    required_shape = (self.num_residues(count_unresolved=False),)
    if res_mask.shape != required_shape:
      raise ValueError(
          f'res_mask must have shape {required_shape}. Got: {res_mask.shape}.'
      )
    if res_mask.dtype != bool:
      raise ValueError(f'res_mask must have dtype bool. Got: {res_mask.dtype}.')

    filtered_residues = self.present_residues.filter(res_mask)
    atom_mask = np.isin(self._atoms.res_key, filtered_residues.key)
    return self.filter(atom_mask)

  def filter_coords(
      self, coord_predicate: Callable[[np.ndarray], bool]
  ) -> Self:
    """Filter a structure's atoms by a function of their coordinates.

    Args:
      coord_predicate: A boolean function of coordinate vectors (shape (3,)).

    Returns:
      A Structure filtered so that only atoms with coords passing the predicate
      function are present.

    Raises:
      ValueError: If the coords are not shaped (num_atom, 3).
    """
    coords = self.coords
    if coords.ndim != 2 or coords.shape[-1] != 3:
      raise ValueError(
          f'coords should have shape (num_atom, 3). Got {coords.shape}.'
      )
    mask = np.vectorize(coord_predicate, signature='(n)->()')(coords)
    # This use of _apply_atom_index_array is safe because a boolean mask is
    # used, which means the chain/residue/atom ordering will stay unchanged.
    return self._apply_atom_index_array(mask, skip_validation=True)

  def filter_polymers_to_single_atom_per_res(
      self,
      representative_atom_by_chain_type: Mapping[
          str, str
      ] = mmcif_names.RESIDUE_REPRESENTATIVE_ATOMS,
  ) -> Self:
    """Filter to one representative atom per polymer residue, ligands unchanged.

    Args:
      representative_atom_by_chain_type: Chain type str to atom name, only atoms
        with this name will be kept for this chain type. Chains types from the
        structure not found in this mapping will keep all their atoms.

    Returns:
      A Structure filtered so that per chain types, only specified atoms are
      present.
    """
    polymer_chain_keys = self._chains.key[
        string_array.isin(
            self._chains.type, set(representative_atom_by_chain_type)
        )
    ]
    polymer_atoms_mask = np.isin(self._atoms.chain_key, polymer_chain_keys)

    wanted_atom_by_chain_key = {
        chain_key: representative_atom_by_chain_type.get(chain_type, None)
        for chain_key, chain_type in zip(self._chains.key, self._chains.type)
    }
    wanted_atoms = string_array.remap(
        self._atoms.chain_key.astype(object), mapping=wanted_atom_by_chain_key
    )

    representative_polymer_atoms_mask = polymer_atoms_mask & (
        wanted_atoms == self._atoms.name
    )

    return self.filter(representative_polymer_atoms_mask | ~polymer_atoms_mask)

  def drop_non_standard_protein_atoms(self, *, drop_oxt: bool = True) -> Self:
    """Drops non-standard atom names from protein chains.

    Args:
      drop_oxt: If True, also drop terminal oxygens (OXT).

    Returns:
      A new Structure object where the protein chains have been filtered to
        only contain atoms with names listed in `atom_types`
        (including OXT unless `drop_oxt` is `True`). Non-protein chains are
        unaltered.
    """
    allowed_names = set(atom_types.ATOM37)
    if drop_oxt:
      allowed_names = {n for n in allowed_names if n != atom_types.OXT}

    return self.filter_out(
        chain_type=mmcif_names.PROTEIN_CHAIN,
        atom_name=lambda n: string_array.isin(n, allowed_names, invert=True),
    )

  def drop_non_standard_atoms(
      self,
      *,
      ccd: chemical_components.Ccd,
      drop_unk: bool,
      drop_non_ccd: bool,
      drop_terminal_oxygens: bool = False,
  ) -> Self:
    """Drops atoms that are not in the CCD for the given residue type."""

    # We don't remove any atoms in UNL, as it has no standard atoms.
    def _keep(atom_index: int) -> bool:
      atom_name = self._atoms.name[atom_index]
      res_name = self._residues.name[
          self._residues.index_by_key[self._atoms.res_key[atom_index]]
      ]
      if drop_unk and res_name in residue_names.UNKNOWN_TYPES:
        return False
      else:
        return (
            (not drop_non_ccd and not ccd.get(res_name))
            or atom_name in struc_chem_comps.get_res_atom_names(ccd, res_name)
            or res_name == residue_names.UNL
        )

    standard_atom_mask = np.array(
        [_keep(atom_i) for atom_i in range(self.num_atoms)], dtype=bool
    )
    standard_atoms = self.filter(mask=standard_atom_mask)
    if drop_terminal_oxygens:
      standard_atoms = standard_atoms.without_terminal_oxygens()
    return standard_atoms

  def find_chains_with_unknown_sequence(self) -> Sequence[str]:
    """Returns a sequence of chain IDs that contain only unknown residues."""
    unknown_sequences = []
    for start, end in self.iter_chain_ranges():
      try:
        unknown_id = residue_names.UNKNOWN_TYPES.index(self.res_name[start])
        if start + 1 == end or np.all(
            self.res_name[start + 1 : end]
            == residue_names.UNKNOWN_TYPES[unknown_id]
        ):
          unknown_sequences.append(self.chain_id[start])
      except ValueError:
        pass
    return unknown_sequences

  def add_bonds(
      self,
      bonded_atom_pairs: Sequence[
          tuple[tuple[str, int, str], tuple[str, int, str]],
      ],
      bond_type: str | None = None,
  ) -> Self:
    """Returns a structure with new bonds added.

    Args:
      bonded_atom_pairs: A sequence of pairs of atoms, with one pair per bond.
        Each element of the pair is a tuple of (chain_id, res_id, atom_name),
        matching values from the respective fields of this structure. The first
        element is the start atom, and the second atom is the end atom of the
        bond.
      bond_type: This type will be used for all bonds in the structure, where
        type follows PDB scheme, e.g. unknown (?), hydrog, metalc, covale,
        disulf.

    Returns:
      A copy of this structure with the new bonds added. If this structure has
      bonds already then the new bonds are concatenated onto the end of the
      old bonds. NB: bonds are not deduplicated.
    """
    atom_key_lookup: dict[tuple[str, str, None, str], int] = dict(
        zip(self.atom_ids, self._atoms.key, strict=True)
    )

    # iter_atoms returns a 4-tuple (chain_id, res_id, ins_code, atom_name) but
    # the insertion code is always None. It also uses string residue IDs.
    def _to_internal_res_id(
        bonded_atom_id: tuple[str, int, str],
    ) -> tuple[str, str, None, str]:
      return bonded_atom_id[0], str(bonded_atom_id[1]), None, bonded_atom_id[2]

    from_atom_key = []
    dest_atom_key = []
    for from_atom, dest_atom in bonded_atom_pairs:
      from_atom_key.append(atom_key_lookup[_to_internal_res_id(from_atom)])
      dest_atom_key.append(atom_key_lookup[_to_internal_res_id(dest_atom)])
    num_bonds = len(bonded_atom_pairs)
    bonds_key = np.arange(num_bonds, dtype=np.int64)
    from_atom_key = np.array(from_atom_key, dtype=np.int64)
    dest_atom_key = np.array(dest_atom_key, dtype=np.int64)
    all_unk_col = np.array(['?'] * num_bonds, dtype=object)
    if bond_type is None:
      bond_type_col = all_unk_col
    else:
      bond_type_col = np.full((num_bonds,), bond_type, dtype=object)

    max_key = -1 if not self._bonds.size else np.max(self._bonds.key)
    new_bonds = structure_tables.Bonds(
        key=np.concatenate([self._bonds.key, bonds_key + max_key + 1]),
        from_atom_key=np.concatenate(
            [self._bonds.from_atom_key, from_atom_key]
        ),
        dest_atom_key=np.concatenate(
            [self._bonds.dest_atom_key, dest_atom_key]
        ),
        type=np.concatenate([self._bonds.type, bond_type_col]),
        role=np.concatenate([self._bonds.role, all_unk_col]),
    )
    return self.copy_and_update(bonds=new_bonds)

  @property
  def coords(self) -> np.ndarray:
    """A [..., num_atom, 3] shaped array of atom coordinates."""
    return np.stack([self._atoms.x, self._atoms.y, self._atoms.z], axis=-1)

  def chain_single_letter_sequence(
      self, include_missing_residues: bool = True
  ) -> Mapping[str, str]:
    """Returns a mapping from chain ID to a single letter residue sequence.

    Args:
      include_missing_residues: Whether to include residues that have no atoms.
    """
    res_table = (
        self._residues if include_missing_residues else self.present_residues
    )
    residue_chain_boundaries = _get_change_indices(res_table.chain_key)
    boundaries = self._iter_residue_ranges(
        residue_chain_boundaries,
        count_unresolved=include_missing_residues,
    )
    chain_keys = res_table.chain_key[residue_chain_boundaries]
    chain_ids = self._chains.apply_array_to_column('id', chain_keys)
    chain_types = self._chains.apply_array_to_column('type', chain_keys)
    chain_seqs = {}
    for idx, (start, end) in enumerate(boundaries):
      chain_id = chain_ids[idx]
      chain_type = chain_types[idx]
      chain_res = res_table.name[start:end]
      if chain_type in mmcif_names.PEPTIDE_CHAIN_TYPES:
        unknown_default = 'X'
      elif chain_type in mmcif_names.NUCLEIC_ACID_CHAIN_TYPES:
        unknown_default = 'N'
      else:
        chain_seqs[chain_id] = 'X' * chain_res.size
        continue

      chain_res = string_array.remap(
          chain_res,
          mapping=residue_names.CCD_NAME_TO_ONE_LETTER,
          inplace=False,
          default_value=unknown_default,
      )
      chain_seqs[chain_id] = ''.join(chain_res.tolist())

    return chain_seqs

  def polymer_auth_asym_id_to_label_asym_id(
      self,
      *,
      protein: bool = True,
      rna: bool = True,
      dna: bool = True,
      other: bool = True,
  ) -> Mapping[str, str]:
    """Mapping from author chain ID to internal chain ID, polymers only.

    This mapping is well defined only for polymers (protein, DNA, RNA), but not
    for ligands or water.

    E.g. if a structure had the following internal chain IDs (label_asym_id):
      A (protein), B (DNA), C (ligand bound to A), D (ligand bound to A),
      E (ligand bound to B).

    Such structure would have this internal chain ID (label_asym_id) -> author
    chain ID (auth_asym_id) mapping:
      A -> A, B -> B, C -> A, D -> A, E -> B

    This is a bijection only for polymers (A, B), but not for ligands.

    Args:
      protein: Whether to include protein (polypeptide(L)) chains.
      rna: Whether to include RNA chains.
      dna: Whether to include DNA chains.
      other: Whether to include other polymer chains, e.g. RNA/DNA hybrid or
        polypeptide(D). Note that include_other=True must be set in from_mmcif.

    Returns:
      A mapping from author chain ID to the internal (label) chain ID for the
      given polymer types in the Structure, ligands/water are ignored.

    Raises:
      ValueError: If the mapping from internal chain IDs to author chain IDs is
        not a bijection for polymer chains.
    """
    allowed_types = set()
    if protein:
      allowed_types.add(mmcif_names.PROTEIN_CHAIN)
    if rna:
      allowed_types.add(mmcif_names.RNA_CHAIN)
    if dna:
      allowed_types.add(mmcif_names.DNA_CHAIN)
    if other:
      non_standard_chain_types = (
          mmcif_names.POLYMER_CHAIN_TYPES
          - mmcif_names.STANDARD_POLYMER_CHAIN_TYPES
      )
      allowed_types |= non_standard_chain_types

    auth_asym_id_to_label_asym_id = {}
    for chain in self.iter_chains():
      if chain['chain_type'] not in allowed_types:
        continue
      label_asym_id = chain['chain_id']
      auth_asym_id = chain['chain_auth_asym_id']
      # The mapping from author chain id to label chain id is only one-to-one if
      # we restrict our attention to polymers. But check nevertheless.
      if auth_asym_id in auth_asym_id_to_label_asym_id:
        raise ValueError(
            f'Author chain ID "{auth_asym_id}" does not have a unique mapping '
            f'to internal chain ID "{label_asym_id}", it is already mapped to '
            f'"{auth_asym_id_to_label_asym_id[auth_asym_id]}".'
        )
      auth_asym_id_to_label_asym_id[auth_asym_id] = label_asym_id

    return auth_asym_id_to_label_asym_id

  def polymer_author_chain_single_letter_sequence(
      self,
      *,
      include_missing_residues: bool = True,
      protein: bool = True,
      rna: bool = True,
      dna: bool = True,
      other: bool = True,
  ) -> Mapping[str, str]:
    """Mapping from author chain ID to single letter aa sequence, polymers only.

    This mapping is well defined only for polymers (protein, DNA, RNA), but not
    for ligands or water.

    Args:
      include_missing_residues: If True then all residues will be returned for
        each polymer chain present in the structure. This uses the all_residues
        field and will include residues missing due to filtering operations as
        well as e.g. unresolved residues specified in an mmCIF header.
      protein: Whether to include protein (polypeptide(L)) chains.
      rna: Whether to include RNA chains.
      dna: Whether to include DNA chains.
      other: Whether to include other polymer chains, e.g. RNA/DNA hybrid or
        polypeptide(D). Note that include_other=True must be set in from_mmcif.

    Returns:
      A mapping from (author) chain IDs to their single-letter sequences for all
      polymers in the Structure, ligands/water are ignored.

    Raises:
      ValueError: If the mapping from internal chain IDs to author chain IDs is
        not a bijection for polymer chains.
    """
    label_chain_id_to_seq = self.chain_single_letter_sequence(
        include_missing_residues=include_missing_residues
    )
    auth_to_label = self.polymer_auth_asym_id_to_label_asym_id(
        protein=protein, rna=rna, dna=dna, other=other
    )
    return {
        auth: label_chain_id_to_seq[label]
        for auth, label in auth_to_label.items()
    }

  def chain_res_name_sequence(
      self,
      *,
      include_missing_residues: bool = True,
      fix_non_standard_polymer_res: bool = False,
  ) -> Mapping[str, Sequence[str]]:
    """A mapping from internal chain ID to a sequence of residue names.

    The residue names are the full residue names rather than single letter
    codes. For instance, for proteins these are the 3 letter CCD codes.

    Args:
      include_missing_residues: Whether to include residues with no atoms in the
        returned sequences.
      fix_non_standard_polymer_res: Whether to map non standard residues in
        protein / RNA / DNA chains to standard residues (e.g. MSE -> MET) or UNK
        / N if a match is not found.

    Returns:
      A mapping from (internal) chain IDs to a sequence of residue names.
    """
    res_table = (
        self._residues if include_missing_residues else self.present_residues
    )
    residue_chain_boundaries = _get_change_indices(res_table.chain_key)
    boundaries = self._iter_residue_ranges(
        residue_chain_boundaries, count_unresolved=include_missing_residues
    )
    chain_keys = res_table.chain_key[residue_chain_boundaries]
    chain_ids = self._chains.apply_array_to_column('id', chain_keys)
    chain_types = self._chains.apply_array_to_column('type', chain_keys)
    chain_seqs = {}
    for idx, (start, end) in enumerate(boundaries):
      chain_id = chain_ids[idx]
      chain_type = chain_types[idx]
      chain_res = res_table.name[start:end]
      if (
          fix_non_standard_polymer_res
          and chain_type in mmcif_names.POLYMER_CHAIN_TYPES
      ):
        chain_seqs[chain_id] = tuple(
            fix_non_standard_polymer_residues(
                res_names=chain_res, chain_type=chain_type
            )
        )
      else:
        chain_seqs[chain_id] = tuple(chain_res)

    return chain_seqs

  def fix_non_standard_polymer_res(
      self,
      res_mapper: Callable[
          [np.ndarray, str], np.ndarray
      ] = fix_non_standard_polymer_residues,
  ) -> Self:
    """Replaces non-standard polymer residues with standard alternatives or UNK.

    e.g. maps 'ACE' -> 'UNK', 'MSE' -> 'MET'.

    NB: Only fixes the residue names, but does not fix the atom names.
    E.g., 'MSE' will be renamed to 'MET' but its 'SE' atom will not be renamed
    to 'S'. Fixing MSE should be done during conversion from mmcif with the
    `fix_mse_residues` flag.

    Args:
      res_mapper: An optional function that accepts a numpy array of residue
        names and chain_type, and returns an array with fixed res_names. This
        defaults to fix_non_standard_polymer_residues.

    Returns:
      A Structure containing only standard residue types (or 'UNK') in its
      polymer chains.
    """
    fixed_res_name = self._residues.name.copy()
    chain_change_indices = _get_change_indices(self._residues.chain_key)
    for start, end in self._iter_atom_ranges(chain_change_indices):
      chain_key = self._residues.chain_key[start]
      chain_type = self._chains.type[self._chains.index_by_key[chain_key]]
      if chain_type not in mmcif_names.POLYMER_CHAIN_TYPES:
        continue  # We don't need to change anything for non-polymers.
      fixed_res_name[start:end] = res_mapper(
          fixed_res_name[start:end], chain_type
      )
    fixed_residues = self._residues.copy_and_update(name=fixed_res_name)
    return self.copy_and_update(residues=fixed_residues, skip_validation=True)

  @property
  def slice_leading_dims(self) -> '_LeadingDimSlice':
    """Used to create a new Structure by slicing into the leading dimensions.

    Example usage 1:

    ```
    final_state = multi_state_struc.slice_leading_dims[-1]
    ```

    Example usage 2:

    ```
    # Structure has leading batch and time dimensions.
    # Get final 3 time frames from first two batch elements.
    sliced_strucs = batched_trajectories.slice_leading_dims[:2, -3:]
    ```
    """
    return _LeadingDimSlice(self)

  def unstack(self, axis: int = 0) -> Sequence[Self]:
    """Unstacks a multi-model structure into a list of Structures.

    This method is the inverse of `stack`.

    Example usage:
    ```
    strucs = multi_dim_struc.unstack(axis=0)
    ```

    Args:
      axis: The axis to unstack over. The structures in the returned list won't
        have this axis in their coordinate of b-factor fields.

    Returns:
      A list of `Structure`s with length equal to the size of the specified
      axis in the coorinate field arrays.

    Raises:
      IndexError: If axis does not refer to one of the leading dimensions of
        `self.atoms_table.size`.
    """
    ndim = self._atoms.ndim
    if not (-ndim <= axis < ndim):
      raise IndexError(
          f'{axis=} is out of range for atom coordinate fields with {ndim=}.'
      )
    elif axis < 0:
      axis += ndim
    if axis == ndim - 1:
      raise IndexError(
          'axis must refer to one of the leading dimensions, not the final '
          f'dimension. The atom fields have {ndim=} and {axis=} was specified.'
      )
    unstacked = []
    leading_dim_slice = self.slice_leading_dims  # Compute once here.
    for i in range(self._atoms.shape[axis]):
      slice_i = (slice(None),) * axis + (i,)
      unstacked.append(leading_dim_slice[slice_i])
    return unstacked

  def split_by_chain(self) -> Sequence[Self]:
    """Splits a Structure into single-chain Structures, one for each chain.

    The obtained structures can be merged back together into the original
    structure using the `concat` function.

    Returns:
      A list of `Structure`s, one for each chain. The order is the same as the
      chain order in the original Structure.
    """
    return [self.filter(chain_id=chain_id) for chain_id in self.chains]

  def transform_states_to_chains(self) -> Self:
    """Transforms states to chains.

    A multi-state protein structure will be transformed to a multi-chain
    single-state protein structure. Useful for visualising multiples states to
    examine diversity. This structure's coordinate fields must have shape
    `(num_states, num_atoms)`.

    Returns:
      A new `Structure`, based on this structure, but with the multiple states
      now represented as `num_states * num_chains` chains in a
      single-state protein.

    Raises:
      ValueError: If this structure's array fields don't have shape
        `(num_states, num_atoms)`.
    """
    if self._atoms.ndim != 2:
      raise ValueError(
          'Coordinate field tensor must have 2 dimensions: '
          f'(num_states, num_atoms), got {self._atoms.ndim}.'
      )
    return concat(self.unstack(axis=0))

  def merge_chains(
      self,
      *,
      chain_groups: Sequence[Sequence[str]],
      chain_group_ids: Sequence[str] | None = None,
      chain_group_types: Sequence[str] | None = None,
      chain_group_entity_ids: Sequence[str] | None = None,
  ) -> Self:
    """Merges chains in each group into a single chain.

    If a Structure has chains A, B, C, D, E, and
    `merge_chains([[A, C], [B, D], [E]])` is called, the new Structure will have
    3 chains A, B, C, the first being concatenation of A+C, the second B+D, the
    third just the original chain E.

    Args:
      chain_groups: Each group defines what chains should be merged into a
        single chain. The output structure will therefore have len(chain_groups)
        chains. Residue IDs (label and author) are renumbered to preserve
        uniqueness within new chains. Order of chain groups and within each
        group matters.
      chain_group_ids: Optional sequence of new chain IDs for each group. If not
        given, the new internal chain IDs (label_asym_id) are assigned in the
        standard mmCIF order (i.e. A, B, ..., Z, AA, BA, CA, ...). Author chain
        names (auth_asym_id) are set to be equal to the new internal chain IDs.
      chain_group_types: Optional sequence of new chain types for each group. If
        not given, only chains with the same type can be merged.
      chain_group_entity_ids: Optional sequence of new entity IDs for each
        group. If not given, the new internal entity IDs (entity_id) are
        assigned in the standard mmCIF order (i.e. 1, 2, 3, ...). Entity
        descriptions (entity_desc) are set to '.' for each chain.

    Returns:
      A new `Structure` with chains merged together into a single chain within
      each chain group.

    Raises:
      ValueError: If chain_group_ids or chain_group_types are given but don't
        match the length of chain_groups.
      ValueError: If the chain IDs in the flattened chain_groups don't match the
        chain IDs in the Structure.
      ValueError: If chains in any of the groups don't have the same chain type.
    """
    if chain_group_ids and len(chain_group_ids) != len(chain_groups):
      raise ValueError(
          'chain_group_ids must the same length as chain_groups: '
          f'{len(chain_group_ids)=} != {len(chain_groups)=}'
      )
    if chain_group_types and len(chain_group_types) != len(chain_groups):
      raise ValueError(
          'chain_group_types must the same length as chain_groups: '
          f'{len(chain_group_types)=} != {len(chain_groups)=}'
      )
    if chain_group_entity_ids and len(chain_group_entity_ids) != len(
        chain_groups
    ):
      raise ValueError(
          'chain_group_entity_ids must the same length as chain_groups: '
          f'{len(chain_group_entity_ids)=} != {len(chain_groups)=}'
      )

    flattened = sorted(itertools.chain.from_iterable(chain_groups))
    if flattened != sorted(self.chains):
      raise ValueError(
          'IDs in chain groups do not match Structure chain IDs: '
          f'{chain_groups=}, chains={self.chains}'
      )

    new_chain_key_by_chain_id = {}
    for new_chain_key, group_chain_ids in enumerate(chain_groups):
      for chain_id in group_chain_ids:
        new_chain_key_by_chain_id[chain_id] = new_chain_key

    chain_key_remap = {}
    new_chain_type_by_chain_key = {}
    for old_chain_key, old_chain_id, old_chain_type in zip(
        self._chains.key, self._chains.id, self._chains.type
    ):
      new_chain_key = new_chain_key_by_chain_id[old_chain_id]
      chain_key_remap[old_chain_key] = new_chain_key

      if new_chain_key not in new_chain_type_by_chain_key:
        new_chain_type_by_chain_key[new_chain_key] = old_chain_type
      elif not chain_group_types:
        if new_chain_type_by_chain_key[new_chain_key] != old_chain_type:
          bad_types = [
              f'{cid}: {self._chains.type[np.where(self._chains.id == cid)][0]}'
              for cid in chain_groups[new_chain_key]
          ]
          raise ValueError(
              'Inconsistent chain types within group:\n' + '\n'.join(bad_types)
          )

    new_chain_key = np.arange(len(chain_groups), dtype=np.int64)
    if chain_group_ids:
      new_chain_id = np.array(chain_group_ids, dtype=object)
    else:
      new_chain_id = np.array(
          [mmcif.int_id_to_str_id(k) for k in new_chain_key + 1], dtype=object
      )
    if chain_group_types:
      new_chain_type = np.array(chain_group_types, dtype=object)
    else:
      new_chain_type = np.array(
          [new_chain_type_by_chain_key[k] for k in new_chain_key], dtype=object
      )
    if chain_group_entity_ids:
      new_chain_entity_id = np.array(chain_group_entity_ids, dtype=object)
    else:
      new_chain_entity_id = np.char.mod('%d', new_chain_key + 1).astype(object)
    new_chains = structure_tables.Chains(
        key=new_chain_key,
        id=new_chain_id,
        type=new_chain_type,
        auth_asym_id=new_chain_id,
        entity_id=new_chain_entity_id,
        entity_desc=np.full(len(chain_groups), fill_value='.', dtype=object),
    )

    # Remap chain keys and sort residues to match the chain table order.
    new_residues = self._residues.copy_and_remap(chain_key=chain_key_remap)
    new_residues = new_residues.apply_index(
        np.argsort(new_residues.chain_key, kind='stable')
    )
    # Renumber uniquely residues in each chain.
    indices = np.arange(new_residues.chain_key.size, dtype=np.int32)
    new_res_ids = (indices + 1) - np.maximum.accumulate(
        indices * (new_residues.chain_key != np.roll(new_residues.chain_key, 1))
    )
    new_residues = new_residues.copy_and_update(
        id=new_res_ids,
        auth_seq_id=np.char.mod('%d', new_res_ids).astype(object),
    )

    # Remap chain keys and sort atoms to match the chain table order.
    new_atoms = self._atoms.copy_and_remap(chain_key=chain_key_remap)
    new_atoms = new_atoms.apply_index(
        np.argsort(new_atoms.chain_key, kind='stable')
    )

    return self.copy_and_update(
        chains=new_chains,
        residues=new_residues,
        atoms=new_atoms,
        bonds=self._bonds,
    )

  def to_res_arrays(
      self,
      *,
      include_missing_residues: bool,
      atom_order: Mapping[str, int] = atom_types.ATOM37_ORDER,
  ) -> ResArrays:
    """Returns atom-level information in arrays containing a num_res dimension.

    NB: All residues in the structure will appear in the residue dimension but
    atoms will only have a True (1.0) mask value if the residue + atom
    combination is defined in `atom_order`. E.g. for the standard ATOM37_ORDER,
    atoms are guaranteed to be present only for standard protein residues.

    Args:
      include_missing_residues: If True then the res arrays will include rows
        for missing residues where all atoms will be masked out. Otherwise these
        will simply be skipped.
      atom_order: Atom order mapping atom names to their index in the atom
        dimension of the returned arrays. Default is atom_order for proteins,
        choose atom_types.ATOM29_ORDER for nucleics.

    Returns:
      A ResArrays object.
    """
    num_res = self.num_residues(count_unresolved=include_missing_residues)
    atom_type_num = len(atom_order)
    atom_positions = np.zeros((num_res, atom_type_num, 3), dtype=np.float32)
    atom_mask = np.zeros((num_res, atom_type_num), dtype=np.float32)
    atom_b_factor = np.zeros((num_res, atom_type_num), dtype=np.float32)
    atom_occupancy = np.zeros((num_res, atom_type_num), dtype=np.float32)

    all_residues = None if not include_missing_residues else self.all_residues
    for i, atom in enumerate_residues(self.iter_atoms(), all_residues):
      atom_idx = atom_order.get(atom['atom_name'])
      if atom_idx is not None:
        atom_positions[i, atom_idx, 0] = atom['atom_x']
        atom_positions[i, atom_idx, 1] = atom['atom_y']
        atom_positions[i, atom_idx, 2] = atom['atom_z']
        atom_mask[i, atom_idx] = 1.0
        atom_b_factor[i, atom_idx] = atom['atom_b_factor']
        atom_occupancy[i, atom_idx] = atom['atom_occupancy']

    return ResArrays(
        atom_positions=atom_positions,
        atom_mask=atom_mask,
        atom_b_factor=atom_b_factor,
        atom_occupancy=atom_occupancy,
    )

  def to_res_atom_lists(
      self, *, include_missing_residues: bool
  ) -> Sequence[Sequence[Mapping[str, Any]]]:
    """Returns list of atom dictionaries grouped by residue.

    If this is a multi-model structure, each atom will store its fields
    atom_x, atom_y, atom_z, and atom_b_factor as Numpy arrays of shape of the
    leading dimension(s). If this is a single-mode structure, these fields will
    just be scalars.

    Args:
      include_missing_residues: If True, then the output list will contain an
        empty list of atoms for missing residues. Otherwise missing residues
        will simply be skipped.

    Returns:
      A list of size `num_res`. Each element in the list represents atoms of one
      residue. If a residue is present is present, the list will contain an atom
      dictionary for every atom present in that residue. If a residue is missing
      and `include_missing_residues=True`, the list for that missing residue
      will be empty.
    """
    num_res = self.num_residues(count_unresolved=include_missing_residues)
    residue_atoms = [[] for _ in range(num_res)]
    all_residues = None if not include_missing_residues else self.all_residues

    # We could yield directly in this loop but the code would be more complex.
    # Let's optimise if memory usage is an issue.
    for res_index, atom in enumerate_residues(self.iter_atoms(), all_residues):
      residue_atoms[res_index].append(atom)

    return residue_atoms

  def reorder_chains(self, new_order: Sequence[str]) -> Self:
    """Reorders tables so that the label_asym_ids are in the given order.

    This method changes the order of the chains, residues, and atoms tables so
    that they are all consistent with each other. Moreover, it remaps chain keys
    so that they stay monotonically increasing in chains/residues/atoms tables.

    Args:
      new_order: The order in which the chain IDs (label_asym_id) should be.
        This must be a permutation of the current chain IDs.

    Returns:
      A structure with chains reorded.
    """
    if len(new_order) != len(self.chains):
      raise ValueError(
          f'The new number of chains ({len(new_order)}) does not match the '
          f'current number of chains ({len(self.chains)}).'
      )
    new_chain_set = set(new_order)
    if len(new_chain_set) != len(new_order):
      raise ValueError(f'The new order {new_order} contains non-unique IDs.')
    if new_chain_set.symmetric_difference(set(self.chains)):
      raise ValueError(
          f'New chain IDs {new_order} do not match the old {set(self.chains)}'
      )

    if self.chains == tuple(new_order):
      return self  # Shortcut: the new order is the same as the current one.

    desired_chain_id_pos = {chain_id: i for i, chain_id in enumerate(new_order)}

    current_chain_index_order = np.empty(self.num_chains, dtype=np.int64)
    for index, old_chain_id in enumerate(self._chains.id):
      current_chain_index_order[index] = desired_chain_id_pos[old_chain_id]
    chain_reorder = np.argsort(current_chain_index_order, kind='stable')
    chain_key_map = dict(
        zip(self._chains.key[chain_reorder], range(self.num_chains))
    )
    chains = self._chains.apply_index(chain_reorder)
    chains = chains.copy_and_remap(key=chain_key_map)

    # The stable sort keeps the original residue ordering within each chain.
    residues = self._residues.copy_and_remap(chain_key=chain_key_map)
    residue_reorder = np.argsort(residues.chain_key, kind='stable')
    residues = residues.apply_index(residue_reorder)

    # The stable sort keeps the original atom ordering within each chain.
    atoms = self._atoms.copy_and_remap(chain_key=chain_key_map)
    atoms_reorder = np.argsort(atoms.chain_key, kind='stable')
    atoms = atoms.apply_index(atoms_reorder)

    # Bonds unchanged - each references 2 atom keys, hence ordering not defined.
    return self.copy_and_update(chains=chains, residues=residues, atoms=atoms)

  def rename_auth_asym_ids(self, new_id_by_old_id: Mapping[str, str]) -> Self:
    """Returns a new structure with renamed author chain IDs (auth_asym_ids).

    Args:
      new_id_by_old_id: A mapping from original auth_asym_ids to their new
        values. Any auth_asym_ids in this structure that are not in the mapping
        will remain unchanged.

    Raises:
      ValueError: If any two previously distinct polymer chains do not have
        unique names anymore after the rename.
    """
    mapped_chains = self._chains.copy_and_remap(auth_asym_id=new_id_by_old_id)
    mapped_polymer_ids = mapped_chains.filter(
        type=mmcif_names.POLYMER_CHAIN_TYPES
    ).auth_asym_id
    if len(mapped_polymer_ids) != len(set(mapped_polymer_ids)):
      raise ValueError(
          'The new polymer auth_asym_ids are not unique:'
          f' {sorted(mapped_polymer_ids)}.'
      )
    return self.copy_and_update(chains=mapped_chains, skip_validation=True)

  def rename_chain_ids(self, new_id_by_old_id: Mapping[str, str]) -> Self:
    """Returns a new structure with renamed chain IDs (label_asym_ids).

    The chains' auth_asym_ids will be updated to be identical to the chain ID
    since there isn't one unambiguous way to maintain the auth_asym_ids after
    renaming the chain IDs (depending on whether you view the auth_asym_id as
    more strongly associated with a given physical chain, or with a given
    chain ID).

    The residues' auth_seq_id will be updated to be identical to the residue ID
    since they are strongly tied to the original author chain naming and keeping
    them would be misleading.

    Args:
      new_id_by_old_id: A mapping from original chain ID to their new values.
        Any chain IDs in this structure that are not in this mapping will remain
        unchanged.

    Returns:
      A new structure with renamed chains (and bioassembly data if it is
      present).

    Raises:
      ValueError: If any two previously distinct chains do not have unique names
        anymore after the rename.
    """
    new_chain_id = string_array.remap(self._chains.id, new_id_by_old_id)
    if len(new_chain_id) != len(set(new_chain_id)):
      raise ValueError(f"New chain names aren't unique: {sorted(new_chain_id)}")

    # Map label_asym_ids in the bioassembly data.
    if self._bioassembly_data is None:
      new_bioassembly_data = None
    else:
      new_bioassembly_data = self._bioassembly_data.rename_label_asym_ids(
          new_id_by_old_id, present_chains=set(self.present_chains.id)
      )

    # Set author residue IDs to be the string version of internal residue IDs.
    new_residues = self._residues.copy_and_update(
        auth_seq_id=self._residues.id.astype(str).astype(object)
    )

    new_chains = self._chains.copy_and_update(
        id=new_chain_id, auth_asym_id=new_chain_id
    )

    return self.copy_and_update(
        bioassembly_data=new_bioassembly_data,
        chains=new_chains,
        residues=new_residues,
        skip_validation=True,
    )

  @functools.cached_property
  def chains(self) -> tuple[str, ...]:
    """Ordered internal chain IDs (label_asym_id) present in the Structure."""
    return tuple(self._chains.id)

  def rename_res_name(
      self,
      res_name_map: Mapping[str, str],
      fail_if_not_found: bool = True,
  ) -> Self:
    """Returns a copy of this structure with residues renamed.

    Residue names in chemical components data will also be renamed.

    Args:
      res_name_map: A mapping from old residue names to new residue names. Any
        residues that are not in this mapping will be left unchanged.
      fail_if_not_found: Whether to fail if keys in the res_name_map mapping are
        not found in this structure's residues' `name` column.

    Raises:
      ValueError: If `fail_if_not_found=True` and a residue name isn't found in
        the residues table's `name` field.
    """
    res_name_set = set(self._residues.name)
    if fail_if_not_found:
      for res_name in res_name_map:
        if res_name not in res_name_set:
          raise ValueError(f'"{res_name}" not found in this structure.')
    new_residues = self._residues.copy_and_remap(name=res_name_map)

    if self._chemical_components_data is not None:
      chem_comp = {
          res_name_map.get(res_name, res_name): data
          for res_name, data in self._chemical_components_data.chem_comp.items()
      }
      new_chem_comp = struc_chem_comps.ChemicalComponentsData(chem_comp)
    else:
      new_chem_comp = None

    return self.copy_and_update(
        residues=new_residues,
        chemical_components_data=new_chem_comp,
        skip_validation=True,
    )

  def remap_res_id(self, res_id_map: Mapping[str, Mapping[int, int]]) -> Self:
    """Returns a copy of this structure with residue IDs remapped.

    Example structure with 2 chains:
    Chain A: residues 1, 2, 3; chain B: residues 6, 7, 8
    res_id_map: {'A': {1: 1, 2: 5, 3: 6}, 'B': {6: 1, 7: 2, 8: 8}}
    Will result in:
    Chain A: residues 1, 5, 6; chain B: residues 1, 2, 8

    Args:
      res_id_map: A mapping from internal chain ID to a mapping from old residue
        ID to new residue ID. A mapping must be provided for each residue in
        each chain.

    Raises:
      KeyError: If residue ID in a given chain is not found in the mapping for
        that chain.
      ValueError: If residue IDs are not unique in each chain after remapping.
    """
    chain_ids = self._chains.apply_array_to_column(
        column_name='id', arr=self._residues.chain_key
    )
    flat_res_id_map = {}
    for chain_id, chain_res_id_map in res_id_map.items():
      flat_res_id_map.update({
          (chain_id, old_res_id): new_res_id
          for old_res_id, new_res_id in chain_res_id_map.items()
      })

    try:
      new_res_id = string_array.remap_multiple(
          (chain_ids, self._residues.id), flat_res_id_map
      )
    except KeyError as e:
      raise KeyError(
          f'Could not find new residue ID for residue {e} in {res_id_map=}'
      ) from e

    residue_chain_boundaries = _get_change_indices(self._residues.chain_key)
    res_boundaries = self._iter_residue_ranges(
        residue_chain_boundaries, count_unresolved=True
    )
    for idx, (start, end) in enumerate(res_boundaries):
      chain_id = chain_ids[idx]
      chain_res_ids = new_res_id[start:end]
      if len(chain_res_ids) != len(set(chain_res_ids)):
        raise ValueError(
            f'New residue IDs not unique in chain {chain_id}: {chain_res_ids}'
        )

    return self.copy_and_update(
        residues=self._residues.copy_and_update(id=new_res_id.astype(np.int32)),
        skip_validation=True,
    )

  def rename_chains_to_match(
      self,
      other: 'Structure',
      *,
      fuzzy_match_non_standard_res: bool = True,
  ) -> Self:
    """Returns a new structure with renamed chains to match another's.

    Example:
    This structure has chains: {'A': 'DEEP', 'B': 'MIND', 'C': 'MIND'}
    Other structure has chains: {'X': 'DEEP', 'Z': 'MIND', 'Y': 'MIND'}

    After calling this method, you will get a structure that has chains named:
    {'X': 'DEEP', 'Z': 'MIND', Y: 'MIND'}

    Args:
      other: Another `Structure`. This provides the reference chain names that
        is used to rename this structure's chains.
      fuzzy_match_non_standard_res: If True, protein/RNA/DNA chains with the
        same one letter sequence will be matched. e.g. "MET-MET-UNK1" will match
        "MET-MSE-UNK2", since both will be mapped to "MMX". If False, we require
        the full res_names to match.

    Returns:
      A new `Structure`, based on this structure, which has chains renamed to
      match the other structure.
    """
    sequences = self.chain_res_name_sequence(
        include_missing_residues=True,
        fix_non_standard_polymer_res=fuzzy_match_non_standard_res,
    )

    other_sequences = other.chain_res_name_sequence(
        include_missing_residues=True,
        fix_non_standard_polymer_res=fuzzy_match_non_standard_res,
    )

    # Check that the sequences are the same.
    sequence_counts = collections.Counter(sequences.values())
    other_sequence_counts = collections.Counter(other_sequences.values())
    if other_sequence_counts != sequence_counts:
      raise ValueError(
          'The other structure does not have the same sequences\n'
          f' other: {other_sequence_counts}\n self: {sequence_counts}'
      )

    new_decoy_id_by_old_id = {}
    used_chain_ids = set()
    # Sort self keys and take min over other to make matching deterministic.
    # The matching is arbitrary but this helps debugging.
    for self_chain_id, self_seq in sorted(sequences.items()):
      # Find corresponding chains in the other structure.
      other_chain_id = min(
          k
          for k, v in other_sequences.items()
          if v == self_seq and k not in used_chain_ids
      )

      new_decoy_id_by_old_id[self_chain_id] = other_chain_id
      used_chain_ids.add(other_chain_id)

    return self.rename_chain_ids(new_decoy_id_by_old_id)

  def _apply_bioassembly_transform(
      self, transform: bioassemblies.Transform
  ) -> Self:
    """Applies a bioassembly transform to this structure."""
    base_struc = self.filter(chain_id=transform.chain_ids)
    transformed_atoms = base_struc.atoms_table.copy_and_update_coords(
        transform.apply_to_coords(base_struc.coords)
    )
    transformed_chains = base_struc.chains_table.copy_and_remap(
        id=transform.chain_id_rename_map
    )
    # Set the transformed author chain ID to match the label chain ID.
    transformed_chains = transformed_chains.copy_and_update(
        auth_asym_id=transformed_chains.id
    )
    return base_struc.copy_and_update(
        chains=transformed_chains,
        atoms=transformed_atoms,
        skip_validation=True,
    )

  def generate_bioassembly(self, assembly_id: str | None = None) -> Self:
    """Generates a biological assembly as a new `Structure`.

    When no assembly ID is provided this method produces a default assembly.
    If this structure has no `bioassembly_data` then this returns itself
    unchanged. Otherwise a default assembly ID is picked with
    `BioassemblyData.get_default_assembly_id()`.

    Args:
      assembly_id: The assembly ID to generate, or None to generate a default
        bioassembly.

    Returns:
      A new `Structure`, based on this one, representing the specified
      bioassembly. Note that if the bioassembly contains copies of chains
      in the original structure then they will be given new unique chain IDs.

    Raises:
      ValueError: If this structure's `bioassembly_data` is `None` and
      `assembly_id` is not `None`.
    """
    if self._bioassembly_data is None:
      if assembly_id is None:
        return self
      else:
        raise ValueError(
            f'Unset bioassembly_data, cannot generate assembly {assembly_id}'
        )

    if assembly_id is None:
      assembly_id = self._bioassembly_data.get_default_assembly_id()

    transformed_strucs = [
        self._apply_bioassembly_transform(transform)
        for transform in self._bioassembly_data.get_transforms(assembly_id)
    ]

    # We don't need to assign unique chain IDs because the bioassembly
    # transform takes care of remapping chain IDs to be unique.
    concatenated = concat(
        transformed_strucs,
        assign_unique_chain_ids=False,
        assign_unique_entity_ids=False,
    )

    # Copy over all scalar fields (e.g. name, release date, etc.) other than
    # bioassembly_data because it relates only to the pre-transformed structure.
    return concatenated.copy_and_update_globals(
        name=self.name,
        release_date=self.release_date,
        resolution=self.resolution,
        structure_method=self.structure_method,
        bioassembly_data=None,
        chemical_components_data=self.chemical_components_data,
    )

  def _to_mmcif_header(self) -> Mapping[str, Sequence[str]]:
    raw_mmcif = collections.defaultdict(list)
    raw_mmcif['data_'] = [self._name.replace(' ', '-')]
    raw_mmcif['_entry.id'] = [self._name]

    if self._release_date is not None:
      date = [datetime.datetime.strftime(self._release_date, '%Y-%m-%d')]
      raw_mmcif['_pdbx_audit_revision_history.revision_date'] = date
      raw_mmcif['_pdbx_database_status.recvd_initial_deposition_date'] = date

    if self._resolution is not None:
      raw_mmcif['_refine.ls_d_res_high'] = ['%.2f' % self._resolution]

    if self._structure_method is not None:
      for method in self._structure_method.split(','):
        raw_mmcif['_exptl.method'].append(method)

    if self._bioassembly_data is not None:
      raw_mmcif.update(self._bioassembly_data.to_mmcif_dict())

    # Populate chemical components data for all residues of this Structure.
    if self._chemical_components_data:
      raw_mmcif.update(self._chemical_components_data.to_mmcif_dict())

    # Add _software table to store version number used to generate mmCIF.
    # Only required data items are used (+ _software.version).
    raw_mmcif['_software.pdbx_ordinal'] = ['1']
    raw_mmcif['_software.name'] = ['DeepMind Structure Class']
    raw_mmcif['_software.version'] = [self._VERSION]
    raw_mmcif['_software.classification'] = ['other']  # Required.

    return raw_mmcif

  def to_mmcif_dict(
      self,
      *,
      coords_decimal_places: int = _COORDS_DECIMAL_PLACES,
  ) -> mmcif.Mmcif:
    """Returns an Mmcif representing the structure."""
    header = self._to_mmcif_header()
    sequence_tables = structure_tables.to_mmcif_sequence_and_entity_tables(
        self._chains, self._residues, self._atoms.res_key
    )
    atom_and_bond_tables = structure_tables.to_mmcif_atom_site_and_bonds_table(
        chains=self._chains,
        residues=self._residues,
        atoms=self._atoms,
        bonds=self._bonds,
        coords_decimal_places=coords_decimal_places,
    )
    return mmcif.Mmcif({**header, **sequence_tables, **atom_and_bond_tables})

  def to_mmcif(
      self, *, coords_decimal_places: int = _COORDS_DECIMAL_PLACES
  ) -> str:
    """Returns an mmCIF string representing the structure.

    Args:
      coords_decimal_places: The number of decimal places to keep for atom
        coordinates, including trailing zeros.
    """
    return self.to_mmcif_dict(
        coords_decimal_places=coords_decimal_places
    ).to_string()


class _LeadingDimSlice:
  """Helper class for slicing the leading dimensions of a `Structure`.

  Wraps a `Structure` instance and applies a slice operation to the coordinate
  fields and other fields that may have leading dimensions (e.g. b_factor).

  Example usage:
    t0_struc = multi_state_struc.slice_leading_dims[0]
  """

  def __init__(self, struc: Structure):
    self._struc = struc

  def __getitem__(self, *args, **kwargs) -> Structure:
    sliced_atom_cols = {}
    for col_name in structure_tables.Atoms.multimodel_cols:
      if (col := self._struc.atoms_table.get_column(col_name)).ndim > 1:
        sliced_col = col.__getitem__(*args, **kwargs)
        if (
            not sliced_col.shape
            or sliced_col.shape[-1] != self._struc.num_atoms
        ):
          raise ValueError(
              'Coordinate slice cannot change final (atom) dimension.'
          )
        sliced_atom_cols[col_name] = sliced_col
    sliced_atoms = self._struc.atoms_table.copy_and_update(**sliced_atom_cols)
    return self._struc.copy_and_update(atoms=sliced_atoms, skip_validation=True)


def stack(strucs: Sequence[Structure], axis: int = 0) -> Structure:
  """Stacks multiple structures into a single multi-model Structure.

  This function is the inverse of `Structure.unstack()`.

  NB: this function assumes that every structure in `strucs` is identical
  other than the coordinates and b-factors. Under this assumption we can safely
  copy all these identical fields from the first element of strucs w.l.o.g.
  However this is not checked in full detail as full comparison is expensive.
  Instead this only checks that the `atom_name` field is identical, and that
  the coordinates have the same shape.

  Usage example:
  ```
  multi_model_struc = structure.stack(strucs, axis=0)
  ```

  Args:
    strucs: A sequence of structures, each with the same atoms, but they may
      have different coordinates and b-factors. If any b-factors are not None
      then they must have the same shape as each of the coordinate fields.
    axis: The axis in the returned structure that represents the different
      structures in `strucs` and will have size `len(strucs)`. This cannot be
      the final dimension as this is reserved for `num_atoms`.

  Returns:
    A `Structure` with the same atoms as the structures in `strucs` but with
    all of their coordinates stacked into a new leading axis.

  Raises:
    ValueError: If `strucs` is empty.
    ValueError: If `strucs` do not all have the same `atom_name` field.
  """
  if not strucs:
    raise ValueError('Need at least one Structure to stack.')
  struc_0, *other_strucs = strucs
  for i, struc in enumerate(other_strucs, start=1):
    # Check that every structure has the same atom name column.
    # This check is intended to catch cases where the input structures might
    # contain the same atoms, but in different orders. This won't catch every
    # such case, e.g. if these are carbon-alpha-only structures, but should
    # catch most cases.
    if np.any(struc.atoms_table.name != struc_0.atoms_table.name):
      raise ValueError(
          f'strucs[0] and strucs[{i}] have mismatching atom name columns.'
      )

  stacked_atoms = struc_0.atoms_table.copy_and_update(
      x=np.stack([s.atoms_table.x for s in strucs], axis=axis),
      y=np.stack([s.atoms_table.y for s in strucs], axis=axis),
      z=np.stack([s.atoms_table.z for s in strucs], axis=axis),
      b_factor=np.stack([s.atoms_table.b_factor for s in strucs], axis=axis),
      occupancy=np.stack([s.atoms_table.occupancy for s in strucs], axis=axis),
  )
  return struc_0.copy_and_update(atoms=stacked_atoms, skip_validation=True)


def _assign_unique_chain_ids(
    strucs: Iterable[Structure],
) -> Sequence[Structure]:
  """Creates a sequence of `Structure` objects with unique chain IDs.

  Let e.g. [A, B] denote a structure of two chains A and B, then this function
  performs the following kind of renaming operation:

  e.g.: [Z], [C], [B, C] -> [A], [B], [C, D]

  NB: This function uses Structure.rename_chain_ids which will define each
  structure's chains.auth_asym_id to be identical to its chains.id columns.

  Args:
    strucs: Structures whose chains ids are to be uniquified.

  Returns:
    A sequence with the same number of elements as `strucs` but where each
    element has had its chains renamed so that they aren't shared with any
    other `Structure` in the sequence.
  """
  # Start counting at 1 because mmcif.int_id_to_str_id expects integers >= 1.
  chain_counter = 1
  strucs_with_new_chain_ids = []
  for struc in strucs:
    rename_map = {}
    for chain_id in struc.chains:
      rename_map[chain_id] = mmcif.int_id_to_str_id(chain_counter)
      chain_counter += 1
    renamed = struc.rename_chain_ids(rename_map)
    strucs_with_new_chain_ids.append(renamed)
  return strucs_with_new_chain_ids


def concat(
    strucs: Sequence[Structure],
    *,
    name: str | None = None,
    assign_unique_chain_ids: bool = True,
    assign_unique_entity_ids: bool = True,
) -> Structure:
  """Concatenates structures along the atom dimension.

  NB: By default this function will first assign unique chain IDs to all chains
  in `strucs` so that the resulting structure does not contain duplicate chain
  IDs. This will also fix entity IDs and author chain IDs. If this is disabled
  via `assign_unique_chain_ids=False` the user must ensure that there are no
  duplicate chains (label_asym_id). However, duplicate entity IDs and author
  chain IDs are allowed as that might be the desired behavior.

  If `assign_unique_chain_ids=True`, note also that the chain_ids may be
  overwritten even if they are already unique.

  Let e.g. [A, B] denote a structure of two chains A and B, then this function
  performs the following kind of concatenation operation:

  assign_unique_chain_ids=True:
    label chain IDS : [Z], [C], [B, C] -> [A, B, C, D]
    author chain IDS: [U], [V], [V, C] -> [A, B, C, D]
    entity IDs      : [1], [1], [3, 3] -> [1, 2, 3, 4]
  assign_unique_chain_ids=False:
    label chain IDS : [D], [B], [C, A] -> [D, B, C, A]  (inputs must be unique)
    author chain IDS: [U], [V], [V, A] -> [U, V, V, A]
    entity IDs      : [1], [1], [3, 3] -> [1, 1, 3, 3]

  NB: This operation loses some information from the elements of `strucs`,
  namely the `name`, `resolution`, `release_date` and `bioassembly_data` fields.

  Args:
    strucs: The `Structure` instances to concatenate. These should all have the
      same number and shape of leading dimensions (i.e. if any are multi-model
      structures then they should all have the same number of models).
    name: Optional name to give to the concatenated structure. If None, the name
      will be concatenation of names of all concatenated structures.
    assign_unique_chain_ids: If True, this function first assigns new unique
      chain IDs, entity IDs and author chain IDs to every chain in `strucs`. If
      False, you must ensure chain IDs are already unique, otherwise an
      exception is raised. See `_assign_unique_chain_ids` for more information
      on how this is performed.
    assign_unique_entity_ids: If True, this function first assigns new unique
      entity IDs to every chain in `strucs`. If False, you must ensure entity
      IDs are already set in a way so that same entity ID implies for two chains
      in `strucs` that they have the same residues. This option applies only if
      `assign_unique_chain_ids == False`, otherwise it must be set to True.

  Returns:
    A new concatenated `Structure` with all of the chains in `strucs` combined
    into one new structure. The new structure will be named by joining the
    names of `strucs` with underscores.

  Raises:
    ValueError: If `strucs` is empty.
    ValueError: If `assign_unique_chain_ids=False` and not all chains in
      `strucs` have unique chain IDs.
  """
  if not strucs:
    raise ValueError('Need at least one Structure to concatenate.')

  if assign_unique_chain_ids and not assign_unique_entity_ids:
    raise ValueError(
        'If assign_unique_chain_ids is True, assign_unique_entity_ids must be '
        'True as well.'
    )

  if assign_unique_chain_ids:
    strucs = _assign_unique_chain_ids(strucs)

  chemical_components_data = {}
  seen_label_chain_ids = set()
  for i, struc in enumerate(strucs):
    if not assign_unique_chain_ids:
      if seen_cid := seen_label_chain_ids.intersection(struc.chains):
        raise ValueError(
            f'Chain IDs {seen_cid} from strucs[{i}] also exist in other'
            ' members of strucs. All given structures must have unique chain'
            ' IDs. Consider setting assign_unique_chain_ids=True.'
        )
      seen_label_chain_ids.update(struc.chains)

    if struc.chemical_components_data is not None:
      chemical_components_data.update(struc.chemical_components_data.chem_comp)  # pytype: disable=attribute-error  # always-use-property-annotation

  concatted_struc = table.concat_databases(strucs)
  name = name if name is not None else '_'.join(s.name for s in strucs)
  # Chain IDs (label and author) are fixed at this point, fix also entity IDs.
  if assign_unique_chain_ids or assign_unique_entity_ids:
    numeric_ids = np.arange(1, concatted_struc.num_chains + 1)
    entity_id = np.char.mod('%d', numeric_ids).astype(object)
    chains = concatted_struc.chains_table.copy_and_update(entity_id=entity_id)
  else:
    chains = concatted_struc.chains_table

  return concatted_struc.copy_and_update(
      name=name,
      release_date=None,
      resolution=None,
      structure_method=None,
      bioassembly_data=None,
      chemical_components_data=(
          struc_chem_comps.ChemicalComponentsData(chemical_components_data)
          if chemical_components_data
          else None
      ),
      chains=chains,
      skip_validation=True,  # Already validated by table.concat_databases.
  )


def multichain_residue_index(
    struc: Structure, chain_offset: int = 9000, between_chain_buffer: int = 1000
) -> np.ndarray:
  """Compute a residue index array that is monotonic across all chains.

  Lots of metrics (lddt, l1_long, etc) require computing a
  distance-along-chain between two residues.  For multimers we want to ensure
  that any residues on different chains have a high along-chain distance
  (i.e. they should always count as long-range contacts for example).  To
  do this we add 10000 to the residue indices of each chain, and enforce that
  the residue index is monotonically increasing across the whole complex.

  Note: This returns the same as struc.res_id for monomers.

  Args:
    struc: The structure to make a multichain residue index for.
    chain_offset: The start of each chain is offset by at least this amount.
      This must be larger than the absolute range of standard residue IDs.
    between_chain_buffer: The final residue in one chain will have at least this
      much of a buffer before the first residue in the next chain.

  Returns:
    A monotonically increasing residue index, with at least
    `between_chain_buffer` residues in between each chain.
  """
  if struc.num_atoms:
    res_id_range = np.max(struc.res_id) - np.min(struc.res_id)
    assert res_id_range < chain_offset
  chain_id_int = struc.chain_id
  monotonic_chain_id_int = np.concatenate(
      ([0], np.cumsum(chain_id_int[1:] != chain_id_int[:-1]))
  )
  return struc.res_id + monotonic_chain_id_int * (
      chain_offset + between_chain_buffer
  )


def make_empty_structure() -> Structure:
  """Returns a new structure consisting of empty array fields."""
  return Structure(
      chains=structure_tables.Chains.make_empty(),
      residues=structure_tables.Residues.make_empty(),
      atoms=structure_tables.Atoms.make_empty(),
      bonds=structure_tables.Bonds.make_empty(),
  )


def enumerate_residues(
    atom_iter: Iterable[Mapping[str, Any]],
    all_residues: AllResidues | None = None,
) -> Iterator[tuple[int, Mapping[str, Any]]]:
  """Provides a zero-indexed enumeration of residues in an atom iterable.

  Args:
    atom_iter: An iterable of atom dicts as returned by Structure.iter_atoms().
    all_residues: (Optional) A structure's all_residues field. If present then
      this will be used to count missing residues by adding appropriate gaps in
      the residue enumeration.

  Yields:
    (res_i, atom) pairs where atom is the unmodified atom dict and res_i is a
    zero-based index for the residue that the atom belongs to.
  """
  if all_residues is None:
    prev_res = None
    res_i = -1
    for atom in atom_iter:
      res = (atom['chain_id'], atom['res_id'])
      if res != prev_res:
        prev_res = res
        res_i += 1
      yield res_i, atom
  else:
    all_res_seq = []  # Sequence of (chain_id, res_id) for all chains.
    prev_chain = None
    res_i = 0
    for atom in atom_iter:
      chain_id = atom['chain_id']
      if chain_id not in all_residues:
        raise ValueError(
            f'Atom {atom} does not belong to any residue in all_residues.'
        )
      if chain_id != prev_chain:
        prev_chain = chain_id
        all_res_seq.extend(
            (chain_id, res_id) for (_, res_id) in all_residues[chain_id]
        )
      res = (chain_id, atom['res_id'])
      while res_i < len(all_res_seq) and res != all_res_seq[res_i]:
        res_i += 1
      if res_i == len(all_res_seq):
        raise ValueError(
            f'Atom {atom} does not belong to a residue in all_residues.'
        )
      yield res_i, atom
