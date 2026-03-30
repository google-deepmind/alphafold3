# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Bond representation for structure module."""

import collections
from collections.abc import Mapping, Sequence
import dataclasses
import typing
from typing import Self

from alphafold3.structure import table
import numpy as np


@dataclasses.dataclass(frozen=True, kw_only=True)
class Bonds(table.Table):
  """Table of atomic bonds."""

  # mmCIF column: _struct_conn.conn_type_id
  # mmCIF desc: This data item is a pointer to _struct_conn_type.id in the
  #             STRUCT_CONN_TYPE category.
  # E.g.: "covale", "disulf", "hydrog", "metalc".
  type: np.ndarray

  # mmCIF column: _struct_conn.pdbx_role
  # mmCIF desc: The chemical or structural role of the interaction.
  # E.g.: "N-Glycosylation", "O-Glycosylation".
  role: np.ndarray

  # mmCIF columns: _struct_conn.ptnr1_*
  from_atom_key: np.ndarray

  # mmCIF columns: _struct_conn.ptnr2_*
  dest_atom_key: np.ndarray

  @classmethod
  def make_empty(cls) -> Self:
    return cls(
        key=np.empty((0,), dtype=np.int64),
        from_atom_key=np.empty((0,), dtype=np.int64),
        dest_atom_key=np.empty((0,), dtype=np.int64),
        type=np.empty((0,), dtype=object),
        role=np.empty((0,), dtype=object),
    )

  def get_atom_indices(
      self,
      atom_key: np.ndarray,
  ) -> tuple[np.ndarray, np.ndarray]:
    """Returns the indices of the from/dest atoms in the atom_key array."""
    from_atom_missing = ~np.isin(self.from_atom_key, atom_key)
    dest_atom_missing = ~np.isin(self.dest_atom_key, atom_key)
    if np.any(from_atom_missing):
      raise ValueError(
          f'No atoms for from_atom_key {self.from_atom_key[from_atom_missing]}'
      )
    if np.any(dest_atom_missing):
      raise ValueError(
          f'No atoms for dest_atom_key {self.dest_atom_key[dest_atom_missing]}'
      )
    sort_indices = np.argsort(atom_key)
    from_indices_sorted = np.searchsorted(
        atom_key, self.from_atom_key, sorter=sort_indices
    )
    dest_indices_sorted = np.searchsorted(
        atom_key, self.dest_atom_key, sorter=sort_indices
    )
    from_indices = sort_indices[from_indices_sorted]
    dest_indices = sort_indices[dest_indices_sorted]
    return from_indices, dest_indices

  def restrict_to_atoms(self, atom_key: np.ndarray) -> Self:
    if not self.size:  # Early-out for empty table.
      return self
    from_atom_mask = np.isin(self.from_atom_key, atom_key)
    dest_atom_mask = np.isin(self.dest_atom_key, atom_key)
    mask = np.logical_and(from_atom_mask, dest_atom_mask)
    return typing.cast(Bonds, self.filter(mask=mask))

  def to_mmcif_dict_from_atom_arrays(
      self,
      atom_key: np.ndarray,
      chain_id: np.ndarray,
      res_id: np.ndarray,
      res_name: np.ndarray,
      atom_name: np.ndarray,
      auth_asym_id: np.ndarray,
      auth_seq_id: np.ndarray,
      insertion_code: np.ndarray,
  ) -> Mapping[str, Sequence[str] | np.ndarray]:
    """Returns a dict suitable for building a CifDict, representing bonds.

    Args:
      atom_key: A (num_atom,) integer array of atom_keys.
      chain_id: A (num_atom,) array of label_asym_id strings.
      res_id: A (num_atom,) array of label_seq_id strings.
      res_name: A (num_atom,) array of label_comp_id strings.
      atom_name: A (num_atom,) array of label_atom_id strings.
      auth_asym_id: A (num_atom,) array of auth_asym_id strings.
      auth_seq_id: A (num_atom,) array of auth_seq_id strings.
      insertion_code: A (num_atom,) array of insertion code strings.
    """
    cif = collections.defaultdict(list)
    p1_indices, p2_indices = self.get_atom_indices(atom_key)

    cif['_struct_conn.ptnr1_label_asym_id'] = chain_id[p1_indices].tolist()
    cif['_struct_conn.ptnr2_label_asym_id'] = chain_id[p2_indices].tolist()
    cif['_struct_conn.ptnr1_label_comp_id'] = res_name[p1_indices].tolist()
    cif['_struct_conn.ptnr2_label_comp_id'] = res_name[p2_indices].tolist()
    cif['_struct_conn.ptnr1_label_seq_id'] = res_id[p1_indices].tolist()
    cif['_struct_conn.ptnr2_label_seq_id'] = res_id[p2_indices].tolist()
    cif['_struct_conn.ptnr1_label_atom_id'] = atom_name[p1_indices].tolist()
    cif['_struct_conn.ptnr2_label_atom_id'] = atom_name[p2_indices].tolist()

    cif['_struct_conn.ptnr1_auth_asym_id'] = auth_asym_id[p1_indices].tolist()
    cif['_struct_conn.ptnr2_auth_asym_id'] = auth_asym_id[p2_indices].tolist()
    cif['_struct_conn.ptnr1_auth_seq_id'] = auth_seq_id[p1_indices].tolist()
    cif['_struct_conn.ptnr2_auth_seq_id'] = auth_seq_id[p2_indices].tolist()
    cif['_struct_conn.pdbx_ptnr1_PDB_ins_code'] = insertion_code[
        p1_indices
    ].tolist()
    cif['_struct_conn.pdbx_ptnr2_PDB_ins_code'] = insertion_code[
        p2_indices
    ].tolist()

    label_alt_id = ['?'] * self.size
    cif['_struct_conn.pdbx_ptnr1_label_alt_id'] = label_alt_id
    cif['_struct_conn.pdbx_ptnr2_label_alt_id'] = label_alt_id

    # We need to set this to make visualisation work in NGL/PyMOL.
    cif['_struct_conn.pdbx_value_order'] = ['?'] * self.size

    # We use a symmetry of 1_555 which is the no-op transformation. Other
    # values are used when bonds involve atoms that only exist after expanding
    # the bioassembly, but we don't support this kind of bond at the moment.
    symmetry = ['1_555'] * self.size
    cif['_struct_conn.ptnr1_symmetry'] = symmetry
    cif['_struct_conn.ptnr2_symmetry'] = symmetry
    bond_type_counter = collections.Counter()
    for bond_row in self.iterrows():
      bond_type = bond_row['type']
      bond_type_counter[bond_type] += 1
      cif['_struct_conn.id'].append(
          f'{bond_type}{bond_type_counter[bond_type]}'
      )
      cif['_struct_conn.pdbx_role'].append(bond_row['role'])
      cif['_struct_conn.conn_type_id'].append(bond_type)

    bond_types = np.unique(self.type)
    cif['_struct_conn_type.id'] = bond_types
    unknown = ['?'] * len(bond_types)
    cif['_struct_conn_type.criteria'] = unknown
    cif['_struct_conn_type.reference'] = unknown

    return dict(cif)


def concat_with_atom_keys(
    bonds_tables: Sequence[Bonds | None],
    atom_key_arrays: Sequence[np.ndarray],
) -> tuple[Bonds | None, np.ndarray]:
  """Concatenates bonds tables and atom keys simultaneously.

  Args:
    bonds_tables: A sequence of `Bonds` instances to concatenate. If any are
      None then these are skipped.
    atom_key_arrays: A sequence of integer `atom_key` arrays, where the n-th
      bonds_table referrs to the atoms in the n-th atom_key array. These must
      all be non-None.

  Returns:
    A pair of (bonds, atom_key) where atom_key is a unique atom_key array with
    length equal to the sum of the input atom array sizes, and the bonds table
    contains all the bonds from the individual bonds table inputs.
  """
  if not bonds_tables or not atom_key_arrays:
    if bonds_tables or atom_key_arrays:
      raise ValueError(
          'bonds_tables and atom_keys must have same length but got'
          f' {len(bonds_tables)=} and {len(atom_key_arrays)=}'
      )
    return None, np.array([], dtype=np.int64)
  max_key = -1
  atom_keys_to_concat = []
  types_to_concat = []
  roles_to_concat = []
  from_atom_keys_to_concat = []
  dest_atom_keys_to_concat = []
  for bonds, atom_key in zip(bonds_tables, atom_key_arrays, strict=True):
    if not atom_key.size:
      assert bonds is None or bonds.size == 0
      continue
    assert np.min(atom_key, initial=0) >= 0  # Should always be non-negative!
    offset = max_key + 1
    offset_atom_key = atom_key + offset
    atom_keys_to_concat.append(offset_atom_key)
    max_key = np.max(offset_atom_key)
    if bonds is not None:
      types_to_concat.append(bonds.type)
      roles_to_concat.append(bonds.role)
      from_atom_keys_to_concat.append(bonds.from_atom_key + offset)
      dest_atom_keys_to_concat.append(bonds.dest_atom_key + offset)

  if atom_keys_to_concat:
    concatted_atom_keys = np.concatenate(atom_keys_to_concat, axis=0)
  else:
    concatted_atom_keys = np.array([], dtype=np.int64)

  if types_to_concat:
    assert (
        len(types_to_concat)
        == len(roles_to_concat)
        == len(from_atom_keys_to_concat)
        == len(dest_atom_keys_to_concat)
    )
    num_bonds = sum(b.size for b in bonds_tables if b is not None)
    concatted_bonds = Bonds(
        key=np.arange(num_bonds, dtype=np.int64),
        type=np.concatenate(types_to_concat, axis=0),
        role=np.concatenate(roles_to_concat, axis=0),
        from_atom_key=np.concatenate(from_atom_keys_to_concat, axis=0),
        dest_atom_key=np.concatenate(dest_atom_keys_to_concat, axis=0),
    )
  else:
    concatted_bonds = None

  return concatted_bonds, concatted_atom_keys
