# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Table implementations for the Structure class."""

import collections
from collections.abc import Mapping, Sequence
import dataclasses
import functools
import itertools
import typing
from typing import Any, ClassVar, Self

from alphafold3.constants import mmcif_names
from alphafold3.constants import residue_names
from alphafold3.cpp import aggregation
from alphafold3.cpp import string_array
from alphafold3.structure import bonds as bonds_module
from alphafold3.structure import mmcif
from alphafold3.structure import table
import numpy as np


Bonds = bonds_module.Bonds


def _residue_name_to_record_name(
    residue_name: np.ndarray,
    polymer_mask: np.ndarray,
) -> np.ndarray:
  """Returns record names (ATOM/HETATM) given residue names and polymer mask."""
  record_name = np.array(['HETATM'] * len(residue_name), dtype=object)
  record_name[polymer_mask] = string_array.remap(
      residue_name[polymer_mask],
      mapping={r: 'ATOM' for r in residue_names.STANDARD_POLYMER_TYPES},
      default_value='HETATM',
  )
  return record_name


@dataclasses.dataclass(frozen=True, slots=True, kw_only=True)
class AuthorNamingScheme:
  """A mapping from internal values to author values in a mmCIF.

  Fields:
    auth_asym_id: A mapping from label_asym_id to auth_asym_id.
    auth_seq_id: A mapping from label_asym_id to a mapping from
      label_seq_id to auth_seq_id.
    insertion_code: A mapping from label_asym_id to a mapping from
      label_seq_id to insertion codes.
    entity_id: A mapping from label_asym_id to _entity.id.
    entity_desc: A mapping from _entity.id to _entity.pdbx_description.
  """

  auth_asym_id: Mapping[str, str]
  auth_seq_id: Mapping[str, Mapping[int, str]]
  insertion_code: Mapping[str, Mapping[int, str | None]]
  entity_id: Mapping[str, str]
  entity_desc: Mapping[str, str]


def _default(
    candidate_value: np.ndarray | None, default_value: Sequence[Any], dtype: Any
) -> np.ndarray:
  if candidate_value is None:
    return np.array(default_value, dtype=dtype)
  return np.array(candidate_value, dtype=dtype)


@dataclasses.dataclass(frozen=True, slots=True, kw_only=True)
class Atoms(table.Table):
  """Table of atoms in a Structure."""

  chain_key: np.ndarray
  res_key: np.ndarray
  name: np.ndarray
  element: np.ndarray
  x: np.ndarray
  y: np.ndarray
  z: np.ndarray
  b_factor: np.ndarray
  occupancy: np.ndarray
  multimodel_cols: ClassVar[tuple[str, ...]] = (
      'x',
      'y',
      'z',
      'b_factor',
      'occupancy',
  )

  def __post_init__(self):
    # Validates that the atom coordinates, b-factors and occupancies are finite.
    for column_name in ('x', 'y', 'z', 'b_factor', 'occupancy'):
      column = self.get_column(column_name)
      if not np.isfinite(column).all():
        raise ValueError(
            f'Column {column_name} must not contain NaN/inf values.'
        )
    # super().__post_init__() can't be used as that causes the following error:
    # TypeError: super(type, obj): obj must be an instance or subtype of type
    super(Atoms, self).__post_init__()

  @classmethod
  def make_empty(cls) -> Self:
    return cls(
        key=np.array([], dtype=np.int64),
        chain_key=np.array([], dtype=np.int64),
        res_key=np.array([], dtype=np.int64),
        name=np.array([], dtype=object),
        element=np.array([], dtype=object),
        x=np.array([], dtype=np.float32),
        y=np.array([], dtype=np.float32),
        z=np.array([], dtype=np.float32),
        b_factor=np.array([], dtype=np.float32),
        occupancy=np.array([], dtype=np.float32),
    )

  @classmethod
  def from_defaults(
      cls,
      *,
      chain_key: np.ndarray,
      res_key: np.ndarray,
      key: np.ndarray | None = None,
      name: np.ndarray | None = None,
      element: np.ndarray | None = None,
      x: np.ndarray | None = None,
      y: np.ndarray | None = None,
      z: np.ndarray | None = None,
      b_factor: np.ndarray | None = None,
      occupancy: np.ndarray | None = None,
  ) -> Self:
    """Create an Atoms table with minimal user inputs."""
    num_atoms = len(chain_key)
    if not num_atoms:
      return cls.make_empty()
    return Atoms(
        chain_key=chain_key,
        res_key=res_key,
        key=_default(key, np.arange(num_atoms), np.int64),
        name=_default(name, ['?'] * num_atoms, object),
        element=_default(element, ['?'] * num_atoms, object),
        x=_default(x, [0.0] * num_atoms, np.float32),
        y=_default(y, [0.0] * num_atoms, np.float32),
        z=_default(z, [0.0] * num_atoms, np.float32),
        b_factor=_default(b_factor, [0.0] * num_atoms, np.float32),
        occupancy=_default(occupancy, [1.0] * num_atoms, np.float32),
    )

  def get_value_by_index(
      self, column_name: str, index: int
  ) -> table.TableEntry | np.ndarray:
    if column_name in self.multimodel_cols:
      return self.get_column(column_name)[..., index]
    else:
      return self.get_column(column_name)[index]

  def copy_and_update_coords(self, coords: np.ndarray) -> Self:
    """Returns a copy with the x, y and z columns updated."""
    if coords.shape[-1] != 3:
      raise ValueError(
          f'Expecting 3-dimensional coordinates, got {coords.shape}'
      )
    return typing.cast(
        Atoms,
        self.copy_and_update(
            x=coords[..., 0], y=coords[..., 1], z=coords[..., 2]
        ),
    )

  @property
  def shape(self) -> tuple[int, ...]:
    return self.x.shape

  @property
  def ndim(self) -> int:
    return len(self.shape)

  @functools.cached_property
  def num_models(self) -> int:
    """The number of models of this Structure."""
    leading_dims = self.shape[:-1]
    match leading_dims:
      case ():
        return 1
      case (single_leading_dim_size,):
        return single_leading_dim_size
      case _:
        raise ValueError(
            'num_models not defined for atom tables with more than one '
            'leading dimension.'
        )


@dataclasses.dataclass(frozen=True, slots=True, kw_only=True)
class Residues(table.Table):
  """Table of residues in a Structure."""

  chain_key: np.ndarray
  id: np.ndarray
  name: np.ndarray
  auth_seq_id: np.ndarray
  insertion_code: np.ndarray

  @classmethod
  def make_empty(cls) -> Self:
    return cls(
        key=np.array([], dtype=np.int64),
        chain_key=np.array([], dtype=np.int64),
        id=np.array([], dtype=np.int32),
        name=np.array([], dtype=object),
        auth_seq_id=np.array([], dtype=object),
        insertion_code=np.array([], dtype=object),
    )

  @classmethod
  def from_defaults(
      cls,
      *,
      id: np.ndarray,  # pylint:disable=redefined-builtin
      chain_key: np.ndarray,
      key: np.ndarray | None = None,
      name: np.ndarray | None = None,
      auth_seq_id: np.ndarray | None = None,
      insertion_code: np.ndarray | None = None,
  ) -> Self:
    """Create a Residues table with minimal user inputs."""
    num_res = len(id)
    if not num_res:
      return cls.make_empty()
    return Residues(
        key=_default(key, np.arange(num_res), np.int64),
        id=id,
        chain_key=chain_key,
        name=_default(name, ['UNK'] * num_res, object),
        auth_seq_id=_default(auth_seq_id, id.astype(str), object),
        insertion_code=_default(insertion_code, ['?'] * num_res, object),
    )


@dataclasses.dataclass(frozen=True, slots=True, kw_only=True)
class Chains(table.Table):
  """Table of chains in a Structure."""

  id: np.ndarray
  type: np.ndarray
  auth_asym_id: np.ndarray
  entity_id: np.ndarray
  entity_desc: np.ndarray

  @classmethod
  def make_empty(cls) -> Self:
    return cls(
        key=np.array([], dtype=np.int64),
        id=np.array([], dtype=object),
        type=np.array([], dtype=object),
        auth_asym_id=np.array([], dtype=object),
        entity_id=np.array([], dtype=object),
        entity_desc=np.array([], dtype=object),
    )

  @classmethod
  def from_defaults(
      cls,
      *,
      id: np.ndarray,  # pylint:disable=redefined-builtin
      key: np.ndarray | None = None,
      type: np.ndarray | None = None,  # pylint:disable=redefined-builtin
      auth_asym_id: np.ndarray | None = None,
      entity_id: np.ndarray | None = None,
      entity_desc: np.ndarray | None = None,
  ) -> Self:
    """Create a Chains table with minimal user inputs."""
    num_chains = len(id)
    if not num_chains:
      return cls.make_empty()

    return Chains(
        key=_default(key, np.arange(num_chains), np.int64),
        id=id,
        type=_default(type, [mmcif_names.PROTEIN_CHAIN] * num_chains, object),
        auth_asym_id=_default(auth_asym_id, id, object),
        entity_id=_default(
            entity_id, np.arange(1, num_chains + 1).astype(str), object
        ),
        entity_desc=_default(entity_desc, ['.'] * num_chains, object),
    )


def to_mmcif_sequence_and_entity_tables(
    chains: Chains,
    residues: Residues,
    atom_res_key: np.ndarray,
) -> Mapping[str, Sequence[str]]:
  """Returns raw sequence and entity mmCIF tables."""
  raw_mmcif = collections.defaultdict(list)
  chains_by_entity_id = {}
  written_entity_poly_seq_ids = set()
  present_res_keys = set(atom_res_key)

  # Performance optimisation: Find residue indices for each chain in advance, so
  # that we don't have to do redunant masking work for each chain.
  res_indices_for_chain = aggregation.indices_grouped_by_value(
      residues.chain_key
  )

  for chain in chains.iterrows():
    # Add all chain information to the _struct_asym table.
    chain_id = chain['id']  # Saves multiple dict lookups.
    auth_asym_id = chain['auth_asym_id']
    entity_id = chain['entity_id']
    chains_by_entity_id.setdefault(entity_id, []).append(chain)
    raw_mmcif['_struct_asym.id'].append(chain_id)
    raw_mmcif['_struct_asym.entity_id'].append(entity_id)

    res_chain_indices = res_indices_for_chain[chain['key']]
    chain_type = chain['type']
    is_polymer = chain_type in mmcif_names.POLYMER_CHAIN_TYPES
    is_water = chain_type == mmcif_names.WATER
    is_branched = len(res_chain_indices) > 1 and not is_polymer and not is_water
    write_entity_poly_seq = entity_id not in written_entity_poly_seq_ids

    # Iterate over the individual masked residue table columns, as that doesn't
    # create a copy (only a view), while residues[res_chain_indices] does.
    for res_key, res_name, res_id, pdb_seq_num, res_ins_code in zip(
        residues.key[res_chain_indices],
        residues.name[res_chain_indices],
        residues.id[res_chain_indices],
        residues.auth_seq_id[res_chain_indices],
        residues.insertion_code[res_chain_indices],
        strict=True,
    ):
      is_missing = res_key not in present_res_keys
      str_res_id = str(res_id)
      # While atom_site uses "?" for insertion codes, scheme tables use ".".
      ins_code = (res_ins_code or '.').replace('?', '.')
      auth_seq_num = '?' if is_missing else pdb_seq_num

      if is_polymer:
        raw_mmcif['_pdbx_poly_seq_scheme.asym_id'].append(chain_id)
        raw_mmcif['_pdbx_poly_seq_scheme.entity_id'].append(entity_id)
        raw_mmcif['_pdbx_poly_seq_scheme.seq_id'].append(str_res_id)
        raw_mmcif['_pdbx_poly_seq_scheme.mon_id'].append(res_name)
        raw_mmcif['_pdbx_poly_seq_scheme.pdb_seq_num'].append(pdb_seq_num)
        raw_mmcif['_pdbx_poly_seq_scheme.auth_seq_num'].append(auth_seq_num)
        raw_mmcif['_pdbx_poly_seq_scheme.pdb_strand_id'].append(auth_asym_id)
        raw_mmcif['_pdbx_poly_seq_scheme.pdb_ins_code'].append(ins_code)
        # Structure doesn't support heterogeneous sequences.
        raw_mmcif['_pdbx_poly_seq_scheme.hetero'].append('n')
        if write_entity_poly_seq:
          raw_mmcif['_entity_poly_seq.entity_id'].append(entity_id)
          raw_mmcif['_entity_poly_seq.num'].append(str_res_id)
          raw_mmcif['_entity_poly_seq.mon_id'].append(res_name)
          # Structure doesn't support heterogeneous sequences.
          raw_mmcif['_entity_poly_seq.hetero'].append('n')
          written_entity_poly_seq_ids.add(entity_id)
      elif is_branched:
        raw_mmcif['_pdbx_branch_scheme.asym_id'].append(chain_id)
        raw_mmcif['_pdbx_branch_scheme.entity_id'].append(entity_id)
        raw_mmcif['_pdbx_branch_scheme.mon_id'].append(res_name)
        raw_mmcif['_pdbx_branch_scheme.num'].append(str_res_id)
        raw_mmcif['_pdbx_branch_scheme.pdb_asym_id'].append(auth_asym_id)
        raw_mmcif['_pdbx_branch_scheme.pdb_seq_num'].append(pdb_seq_num)
        raw_mmcif['_pdbx_branch_scheme.auth_asym_id'].append(auth_asym_id)
        raw_mmcif['_pdbx_branch_scheme.auth_seq_num'].append(auth_seq_num)
        raw_mmcif['_pdbx_branch_scheme.pdb_ins_code'].append(ins_code)
        # Structure doesn't support heterogeneous sequences.
        raw_mmcif['_pdbx_branch_scheme.hetero'].append('n')
      else:
        raw_mmcif['_pdbx_nonpoly_scheme.asym_id'].append(chain_id)
        raw_mmcif['_pdbx_nonpoly_scheme.entity_id'].append(entity_id)
        raw_mmcif['_pdbx_nonpoly_scheme.mon_id'].append(res_name)
        raw_mmcif['_pdbx_nonpoly_scheme.pdb_seq_num'].append(pdb_seq_num)
        raw_mmcif['_pdbx_nonpoly_scheme.auth_seq_num'].append(auth_seq_num)
        raw_mmcif['_pdbx_nonpoly_scheme.pdb_strand_id'].append(auth_asym_id)
        raw_mmcif['_pdbx_nonpoly_scheme.pdb_ins_code'].append(ins_code)

  # Add _entity and _entity_poly tables.
  for entity_id, chains in chains_by_entity_id.items():
    # chains should always be a non-empty list because of how we constructed
    # chains_by_entity_id.
    assert chains
    # All chains for a given entity should have the same type and sequence
    # so we can pick the first one without losing information.
    key_chain = chains[0]
    raw_mmcif['_entity.id'].append(entity_id)
    raw_mmcif['_entity.pdbx_description'].append(key_chain['entity_desc'])
    entity_type = key_chain['type']
    if entity_type not in mmcif_names.POLYMER_CHAIN_TYPES:
      raw_mmcif['_entity.type'].append(entity_type)
    else:
      raw_mmcif['_entity.type'].append('polymer')
      raw_mmcif['_entity_poly.entity_id'].append(entity_id)
      raw_mmcif['_entity_poly.type'].append(entity_type)

      # _entity_poly.pdbx_strand_id is a comma-separated list of
      # auth_asym_ids that are part of the entity.
      raw_mmcif['_entity_poly.pdbx_strand_id'].append(
          ','.join(chain['auth_asym_id'] for chain in chains)
      )
  return raw_mmcif


def to_mmcif_atom_site_and_bonds_table(
    *,
    chains: Chains,
    residues: Residues,
    atoms: Atoms,
    bonds: Bonds,
    coords_decimal_places: int,
) -> Mapping[str, Sequence[str]]:
  """Returns raw _atom_site and _struct_conn mmCIF tables."""
  raw_mmcif: dict[str, Sequence[str]] = {}
  # Use [value] * num wherever possible since it is about 10x faster than list
  # comprehension in such cases. Also use f-strings instead of str() - faster.
  total_atoms = atoms.size * atoms.num_models
  raw_mmcif['_atom_site.id'] = [f'{i}' for i in range(1, total_atoms + 1)]
  raw_mmcif['_atom_site.label_alt_id'] = ['.'] * total_atoms
  # Use format_float_array instead of list comprehension for performance.
  raw_mmcif['_atom_site.Cartn_x'] = mmcif.format_float_array(
      values=atoms.x.ravel(), num_decimal_places=coords_decimal_places
  )
  raw_mmcif['_atom_site.Cartn_y'] = mmcif.format_float_array(
      values=atoms.y.ravel(), num_decimal_places=coords_decimal_places
  )
  raw_mmcif['_atom_site.Cartn_z'] = mmcif.format_float_array(
      values=atoms.z.ravel(), num_decimal_places=coords_decimal_places
  )

  # atoms.b_factor or atoms.occupancy can be flat even when the coordinates have
  # leading dimensions. In this case we tile it to match.
  if atoms.b_factor.ndim == 1:
    atom_b_factor = np.tile(atoms.b_factor, atoms.num_models)
  else:
    atom_b_factor = atoms.b_factor.ravel()
  raw_mmcif['_atom_site.B_iso_or_equiv'] = mmcif.format_float_array(
      values=atom_b_factor, num_decimal_places=2
  )

  if atoms.occupancy.ndim == 1:
    atom_occupancy = np.tile(atoms.occupancy, atoms.num_models)
  else:
    atom_occupancy = atoms.occupancy.ravel()
  raw_mmcif['_atom_site.occupancy'] = mmcif.format_float_array(
      values=atom_occupancy.ravel(), num_decimal_places=2
  )

  label_atom_id = atoms.name
  type_symbol = atoms.element
  label_comp_id = residues.apply_array_to_column('name', atoms.res_key)
  label_asym_id = chains.apply_array_to_column('id', atoms.chain_key)
  label_entity_id = chains.apply_array_to_column('entity_id', atoms.chain_key)
  # Performance optimisation: Do the int->str conversion on num_residue-sized,
  # array, then select instead of selecting and then converting.
  label_seq_id = residues.id.astype('str').astype(object)[
      ..., residues.index_by_key[atoms.res_key]
  ]

  # _atom_site.label_seq_id is '.' for non-polymers.
  non_polymer_chain_mask = string_array.isin(
      chains.type, mmcif_names.POLYMER_CHAIN_TYPES, invert=True
  )
  non_polymer_chain_keys = chains.key[non_polymer_chain_mask]
  non_polymer_atom_mask = np.isin(atoms.chain_key, non_polymer_chain_keys)
  label_seq_id[non_polymer_atom_mask] = '.'

  auth_asym_id = chains.apply_array_to_column('auth_asym_id', atoms.chain_key)
  auth_seq_id = residues.apply_array_to_column('auth_seq_id', atoms.res_key)
  pdbx_pdb_ins_code = residues.apply_array_to_column(
      'insertion_code', atoms.res_key
  )
  string_array.remap(pdbx_pdb_ins_code, mapping={None: '?'}, inplace=True)

  group_pdb = _residue_name_to_record_name(
      residue_name=label_comp_id, polymer_mask=~non_polymer_atom_mask
  )

  def tile_for_models(arr: np.ndarray) -> list[str]:
    if atoms.num_models == 1:
      return arr.tolist()  # Memory optimisation: np.tile(arr, 1) does a copy.
    return np.tile(arr, atoms.num_models).tolist()

  raw_mmcif['_atom_site.group_PDB'] = tile_for_models(group_pdb)
  raw_mmcif['_atom_site.label_atom_id'] = tile_for_models(label_atom_id)
  raw_mmcif['_atom_site.type_symbol'] = tile_for_models(type_symbol)
  raw_mmcif['_atom_site.label_comp_id'] = tile_for_models(label_comp_id)
  raw_mmcif['_atom_site.label_asym_id'] = tile_for_models(label_asym_id)
  raw_mmcif['_atom_site.label_entity_id'] = tile_for_models(label_entity_id)
  raw_mmcif['_atom_site.label_seq_id'] = tile_for_models(label_seq_id)
  raw_mmcif['_atom_site.auth_asym_id'] = tile_for_models(auth_asym_id)
  raw_mmcif['_atom_site.auth_seq_id'] = tile_for_models(auth_seq_id)
  raw_mmcif['_atom_site.pdbx_PDB_ins_code'] = tile_for_models(pdbx_pdb_ins_code)
  model_id = np.array(
      [str(i + 1) for i in range(atoms.num_models)], dtype=object
  )
  raw_mmcif['_atom_site.pdbx_PDB_model_num'] = np.repeat(
      model_id, [atoms.size] * atoms.num_models
  ).tolist()

  if bonds.key.size > 0:
    raw_mmcif.update(
        bonds.to_mmcif_dict_from_atom_arrays(
            atom_key=atoms.key,
            chain_id=label_asym_id,
            res_id=label_seq_id,
            res_name=label_comp_id,
            atom_name=label_atom_id,
            auth_asym_id=auth_asym_id,
            auth_seq_id=auth_seq_id,
            insertion_code=np.array(pdbx_pdb_ins_code),
        )
    )
  return raw_mmcif


def _flatten_author_naming_scheme_table(
    res_table: Mapping[str, Mapping[int, str]],
    chain_ids: np.ndarray,
    res_chain_ids: np.ndarray,
    res_ids: np.ndarray,
    default_if_missing: str,
    table_name: str,
) -> np.ndarray:
  """Flattens an author naming scheme table consistently with res_ids."""
  if not set(chain_ids).issubset(res_table):
    raise ValueError(
        f'Chain IDs in the chain_id array must be a subset of {table_name} in '
        'author naming scheme:\n'
        f'chain_ids: {sorted(chain_ids)}\n'
        f'{table_name} keys: {sorted(res_table.keys())}'
    )

  chain_change_mask = res_chain_ids[1:] != res_chain_ids[:-1]
  res_chain_boundaries = np.concatenate(
      ([0], np.where(chain_change_mask)[0] + 1, [len(res_chain_ids)])
  )

  flat_vals = np.empty(len(res_ids), dtype=object)
  for chain_start, chain_end in itertools.pairwise(res_chain_boundaries):
    chain_id = res_chain_ids[chain_start]
    chain_res_ids = res_ids[chain_start:chain_end]
    chain_mapping = res_table[chain_id]
    flat_vals[chain_start:chain_end] = [
        chain_mapping.get(r, default_if_missing) for r in chain_res_ids
    ]

  return flat_vals


def tables_from_atom_arrays(
    *,
    res_id: np.ndarray,
    author_naming_scheme: AuthorNamingScheme | None = None,
    all_residues: Mapping[str, Sequence[tuple[str, int]]] | None = None,
    chain_id: np.ndarray | None = None,
    chain_type: np.ndarray | None = None,
    res_name: np.ndarray | None = None,
    atom_key: np.ndarray | None = None,
    atom_name: np.ndarray | None = None,
    atom_element: np.ndarray | None = None,
    atom_x: np.ndarray | None = None,
    atom_y: np.ndarray | None = None,
    atom_z: np.ndarray | None = None,
    atom_b_factor: np.ndarray | None = None,
    atom_occupancy: np.ndarray | None = None,
) -> tuple[Atoms, Residues, Chains]:
  """Returns Structure tables constructed from atom array level data.

  All fields except name and, res_id are optional, all array fields consist of a
  value for each atom in the structure - so residue and chain values should hold
  the same value for each atom in the chain or residue. Fields which are not
  defined are filled with default values.

  Validation is performed by the Structure constructor where possible - but
  author_naming scheme and all_residues must be checked in this function.

  It is not possible to construct structures with chains that do not contain
  any resolved residues using this function. If this is necessary, use the
  structure.Structure constructor directly.

  Args:
    res_id: Integer array of shape [num_atom]. The unique residue identifier for
      each residue. mmCIF field - _atom_site.label_seq_id.
    author_naming_scheme: An optional instance of AuthorNamingScheme to use when
      converting this structure to mmCIF.
    all_residues: An optional mapping from each chain ID (i.e. label_asym_id) to
      a sequence of (label_comp_id, label_seq_id) tuples, one per residue. This
      can contain residues that aren't present in the atom arrays. This is
      common in experimental data where some residues are not resolved but are
      known to be present.
    chain_id: String array of shape [num_atom] of unique chain identifiers.
      mmCIF field - _atom_site.label_asym_id.
    chain_type: String array of shape [num_atom]. The molecular type of the
      current chain (e.g. polyribonucleotide). mmCIF field - _entity_poly.type
      OR _entity.type (for non-polymers).
    res_name: String array of shape [num_atom].. The name of each residue,
      typically a 3 letter string for polypeptides or 1-2 letter strings for
      polynucleotides. mmCIF field - _atom_site.label_comp_id.
    atom_key: A unique sorted integer array, used only by the bonds table to
      identify the atoms participating in each bond. If the bonds table is
      specified then this column must be non-None.
    atom_name: String array of shape [num_atom]. The name of each atom (e.g CA,
      O2', etc.). mmCIF field - _atom_site.label_atom_id.
    atom_element: String array of shape [num_atom]. The element type of each
      atom (e.g. C, O, N, etc.). mmCIF field - _atom_site.type_symbol.
    atom_x: Float array of shape [..., num_atom] of atom x coordinates. May have
      arbitrary leading dimensions, provided that these are consistent across
      all coordinate fields.
    atom_y: Float array of shape [..., num_atom] of atom y coordinates. May have
      arbitrary leading dimensions, provided that these are consistent across
      all coordinate fields.
    atom_z: Float array of shape [..., num_atom] of atom z coordinates. May have
      arbitrary leading dimensions, provided that these are consistent across
      all coordinate fields.
    atom_b_factor: Float array of shape [..., num_atom] or [num_atom] of atom
      b-factors or equivalent. If there are no extra leading dimensions then
      these values are assumed to apply to all coordinates for a given atom. If
      there are leading dimensions then these must match those used by the
      coordinate fields.
    atom_occupancy: Float array of shape [..., num_atom] or [num_atom] of atom
      occupancies or equivalent. If there are no extra leading dimensions then
      these values are assumed to apply to all coordinates for a given atom. If
      there are leading dimensions then these must match those used by the
      coordinate fields.
  """
  num_atoms = len(res_id)

  for arr_name, array, dtype in (
      ('chain_id', chain_id, object),
      ('chain_type', chain_type, object),
      ('res_id', res_id, np.int32),
      ('res_name', res_name, object),
      ('atom_key', atom_key, np.int64),
      ('atom_name', atom_name, object),
      ('atom_element', atom_element, object),
  ):
    if array is not None and array.shape != (num_atoms,):
      raise ValueError(f'{arr_name} shape {array.shape} != ({num_atoms},)')
    if array is not None and array.dtype != dtype:
      raise ValueError(f'{arr_name} dtype {array.dtype} != {dtype}')

  for arr_name, array in (
      ('atom_x', atom_x),
      ('atom_y', atom_y),
      ('atom_z', atom_z),
      ('atom_b_factor', atom_b_factor),
      ('atom_occupancy', atom_occupancy),
  ):
    if array is not None and array.shape[-1] != num_atoms:
      raise ValueError(f'{arr_name} last dim {array.shape[-1]} != {num_atoms=}')
    if (
        array is not None
        and array.dtype != np.float32
        and array.dtype != np.float64
    ):
      raise ValueError(
          f'{arr_name} must be np.float32 or np.float64, got {array.dtype=}'
      )

  if all_residues is not None and (res_name is None or res_id is None):
    raise ValueError(
        'If all_residues != None, res_name and res_id must not be None either.'
    )

  if num_atoms == 0:
    return Atoms.make_empty(), Residues.make_empty(), Chains.make_empty()

  if chain_id is None:
    chain_id = np.full(shape=num_atoms, fill_value='A', dtype=object)
  if res_name is None:
    res_name = np.full(shape=num_atoms, fill_value='UNK', dtype=object)

  chain_change_mask = chain_id[1:] != chain_id[:-1]
  chain_start = np.concatenate(([0], np.where(chain_change_mask)[0] + 1))
  res_start = np.concatenate(
      ([0], np.where((res_id[1:] != res_id[:-1]) | chain_change_mask)[0] + 1)
  )

  if len(set(chain_id)) != len(chain_start):
    raise ValueError(f'Chain IDs must be contiguous, but got {chain_id}')

  # We do not support chains with unresolved residues-only in this function.
  chain_ids = chain_id[chain_start]
  if all_residues and set(all_residues.keys()) != set(chain_ids):
    raise ValueError(
        'all_residues must contain the same set of chain IDs as the chain_id '
        f'array:\nall_residues keys: {sorted(all_residues.keys())}\n'
        f'chain_ids: {sorted(chain_ids)}.'
    )
  # Make sure all_residue ordering is consistent with chain_id.
  if all_residues and np.any(list(all_residues.keys()) != chain_ids):
    all_residues = {cid: all_residues[cid] for cid in chain_ids}

  # Create the chains table.
  num_chains = len(chain_ids)
  chain_keys = np.arange(num_chains, dtype=np.int64)
  chain_key_by_chain_id = dict(zip(chain_ids, chain_keys, strict=True))

  if chain_type is not None:
    chain_types = chain_type[chain_start]
  else:
    chain_types = np.full(num_chains, mmcif_names.PROTEIN_CHAIN, dtype=object)

  if author_naming_scheme is not None:
    auth_asym_id = string_array.remap(
        chain_ids, author_naming_scheme.auth_asym_id
    )
    entity_id = string_array.remap(
        chain_ids, author_naming_scheme.entity_id, default_value='.'
    )
    entity_desc = string_array.remap(
        entity_id, author_naming_scheme.entity_desc, default_value='.'
    )
  else:
    auth_asym_id = chain_ids
    entity_id = (chain_keys + 1).astype(str).astype(object)
    entity_desc = np.full(num_chains, '.', dtype=object)

  chains = Chains(
      key=chain_keys,
      id=chain_ids,
      type=chain_types,
      auth_asym_id=auth_asym_id,
      entity_id=entity_id,
      entity_desc=entity_desc,
  )

  # Create the residues table.
  if all_residues is not None:
    residue_order = []
    for cid, residues in all_residues.items():
      residue_order.extend((cid, rname, int(rid)) for (rname, rid) in residues)
    res_chain_ids, res_names, res_ids = zip(*residue_order)
    res_chain_ids = np.array(res_chain_ids, dtype=object)
    res_ids = np.array(res_ids, dtype=np.int32)
    res_names = np.array(res_names, dtype=object)
  else:
    res_chain_ids = chain_id[res_start]
    res_ids = res_id[res_start]
    res_names = res_name[res_start]
    residue_order = list(zip(res_chain_ids, res_names, res_ids))

  if author_naming_scheme is not None and author_naming_scheme.auth_seq_id:
    auth_seq_id = _flatten_author_naming_scheme_table(
        author_naming_scheme.auth_seq_id,
        chain_ids=chain_ids,
        res_chain_ids=res_chain_ids,
        res_ids=res_ids,
        default_if_missing='.',
        table_name='auth_seq_id',
    )
  else:
    auth_seq_id = res_ids.astype(str).astype(object)

  if author_naming_scheme is not None and author_naming_scheme.insertion_code:
    insertion_code = _flatten_author_naming_scheme_table(
        author_naming_scheme.insertion_code,
        chain_ids=chain_ids,
        res_chain_ids=res_chain_ids,
        res_ids=res_ids,
        default_if_missing='?',
        table_name='insertion_code',
    )
    # Make sure insertion code of None is mapped to '.'.
    insertion_code = string_array.remap(insertion_code, {None: '?'})
  else:
    insertion_code = np.full(shape=len(res_ids), fill_value='?', dtype=object)

  res_key_by_res = {res: i for i, res in enumerate(residue_order)}
  res_keys = np.arange(len(residue_order), dtype=np.int64)
  res_chain_keys = string_array.remap(
      res_chain_ids, chain_key_by_chain_id
  ).astype(np.int64)
  residues = Residues(
      chain_key=res_chain_keys,
      key=res_keys,
      id=res_ids,
      name=res_names,
      auth_seq_id=auth_seq_id,
      insertion_code=insertion_code,
  )

  if atom_key is None:
    atom_key = np.arange(num_atoms, dtype=np.int64)

  atom_chain_keys = string_array.remap(chain_id, chain_key_by_chain_id).astype(
      np.int64
  )

  try:
    atom_res_keys = [res_key_by_res[r] for r in zip(chain_id, res_name, res_id)]
  except KeyError as e:
    missing_chain_id, missing_res_name, missing_res_id = e.args[0]
    raise ValueError(
        'Inconsistent res_name, res_id and all_residues. Could not find '
        f'residue with chain_id={missing_chain_id}, '
        f'res_name={missing_res_name}, res_id={missing_res_id} in all_residues.'
    ) from e

  atoms = Atoms(
      key=atom_key,
      chain_key=atom_chain_keys,
      res_key=np.array(atom_res_keys, dtype=np.int64),
      name=_default(atom_name, ['?'] * num_atoms, object),
      element=_default(atom_element, ['?'] * num_atoms, object),
      x=_default(atom_x, [0.0] * num_atoms, np.float32),
      y=_default(atom_y, [0.0] * num_atoms, np.float32),
      z=_default(atom_z, [0.0] * num_atoms, np.float32),
      b_factor=_default(atom_b_factor, [0.0] * num_atoms, np.float32),
      occupancy=_default(atom_occupancy, [1.0] * num_atoms, np.float32),
  )
  return atoms, residues, chains
