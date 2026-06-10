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

"""Utilities for parsing and manipulating bioassembly data."""

from collections.abc import Mapping, Sequence
import copy
import dataclasses
from typing import Self

from alphafold3.structure import mmcif
import numpy as np


@dataclasses.dataclass(frozen=True)
class Operation:
  """A rigid transformation operation."""

  trans: np.ndarray  # shape: (3,)
  rot: np.ndarray  # shape: (3, 3)

  def apply_to_coords(self, coords: np.ndarray) -> np.ndarray:
    """Applies the rotation followed by the translation to `coords`."""
    return np.dot(coords, self.rot.T) + self.trans[np.newaxis, :]


@dataclasses.dataclass(frozen=True)
class Transform:
  """A rigid transformation composed of a sequence of `Operation`s."""

  # The sequence of operations that form the transform. These will be applied
  # right-to-left (last-to-first).
  operations: Sequence[Operation]

  # The chain IDs that this transform should be applied to. These are
  # label_asym_ids in the mmCIF spec.
  chain_ids: Sequence[str]

  # A mapping from chain IDs (of chains that participate in this transform)
  # to their new values in the bioassembly.
  chain_id_rename_map: Mapping[str, str]

  def apply_to_coords(self, coords: np.ndarray) -> np.ndarray:
    """Applies the `operations` in right-to-left order."""
    for operation in reversed(self.operations):
      coords = operation.apply_to_coords(coords)
    return coords


def _get_operation(oper_data: Mapping[str, str]) -> Operation:
  """Parses an `Operation` from a mmCIF _pdbx_struct_oper_list row."""
  trans = np.zeros((3,), dtype=np.float32)
  rot = np.zeros((3, 3), dtype=np.float32)
  for i in range(3):
    trans[i] = float(oper_data[f'_pdbx_struct_oper_list.vector[{i + 1}]'])
  for i in range(3):
    for j in range(3):
      rot[i][j] = float(
          oper_data[f'_pdbx_struct_oper_list.matrix[{i + 1}][{j + 1}]']
      )
  return Operation(trans=trans, rot=rot)


class MissingBioassemblyDataError(Exception):
  """Raised when bioassembly data is missing from an mmCIF."""


class BioassemblyData:
  """Stores and processes bioassembly data from mmCIF tables."""

  # Not all of these columns are required for internal operations, but all
  # should be present whenever bioassemblies are defined in an mmCIF to stay
  # consistent with external mmCIFs.
  _REQUIRED_COLUMNS = (
      '_pdbx_struct_assembly.id',
      '_pdbx_struct_assembly.details',
      '_pdbx_struct_assembly.method_details',
      '_pdbx_struct_assembly.oligomeric_details',
      '_pdbx_struct_assembly.oligomeric_count',
      '_pdbx_struct_assembly_gen.assembly_id',
      '_pdbx_struct_assembly_gen.oper_expression',
      '_pdbx_struct_assembly_gen.asym_id_list',
      '_pdbx_struct_oper_list.id',
      '_pdbx_struct_oper_list.type',
      '_pdbx_struct_oper_list.name',
      '_pdbx_struct_oper_list.symmetry_operation',
      '_pdbx_struct_oper_list.matrix[1][1]',
      '_pdbx_struct_oper_list.matrix[1][2]',
      '_pdbx_struct_oper_list.matrix[1][3]',
      '_pdbx_struct_oper_list.vector[1]',
      '_pdbx_struct_oper_list.matrix[2][1]',
      '_pdbx_struct_oper_list.matrix[2][2]',
      '_pdbx_struct_oper_list.matrix[2][3]',
      '_pdbx_struct_oper_list.vector[2]',
      '_pdbx_struct_oper_list.matrix[3][1]',
      '_pdbx_struct_oper_list.matrix[3][2]',
      '_pdbx_struct_oper_list.matrix[3][3]',
      '_pdbx_struct_oper_list.vector[3]',
  )

  def __init__(
      self,
      *,
      pdbx_struct_assembly: Mapping[str, Mapping[str, str]],
      pdbx_struct_assembly_gen: Mapping[str, Sequence[Mapping[str, str]]],
      pdbx_struct_oper_list: Mapping[str, Mapping[str, str]],
      assembly_ids: Sequence[str],
      oper_ids: Sequence[str],
  ):
    for assembly_id in assembly_ids:
      for table, table_name in (
          (pdbx_struct_assembly, '_pdbx_struct_assembly'),
          (pdbx_struct_assembly_gen, '_pdbx_struct_assembly_gen'),
      ):
        if assembly_id not in table:
          raise ValueError(
              f'Assembly ID "{assembly_id}" missing from {table_name} '
              f'with keys: {table.keys()}'
          )
    for oper_id in oper_ids:
      if oper_id not in pdbx_struct_oper_list:
        raise ValueError(
            f'Oper ID "{oper_id}" missing from _pdbx_struct_oper_list '
            f'with keys: {pdbx_struct_oper_list.keys()}'
        )

    self._pdbx_struct_assembly = pdbx_struct_assembly
    self._pdbx_struct_assembly_gen = pdbx_struct_assembly_gen
    self._pdbx_struct_oper_list = pdbx_struct_oper_list
    self._operations = {
        oper_id: _get_operation(oper_data)
        for oper_id, oper_data in self._pdbx_struct_oper_list.items()
    }
    self._assembly_ids = assembly_ids
    self._oper_ids = oper_ids

  @classmethod
  def from_mmcif(cls, cif: mmcif.Mmcif) -> Self:
    """Constructs an instance of `BioassemblyData` from an `Mmcif` object."""
    for col in cls._REQUIRED_COLUMNS:
      if col not in cif:
        raise MissingBioassemblyDataError(col)

    pdbx_struct_assembly = cif.extract_loop_as_dict(
        prefix='_pdbx_struct_assembly.', index='_pdbx_struct_assembly.id'
    )
    pdbx_struct_oper_list = cif.extract_loop_as_dict(
        prefix='_pdbx_struct_oper_list.', index='_pdbx_struct_oper_list.id'
    )

    # _pdbx_struct_assembly_gen is unlike the other two tables because it can
    # have multiple rows share the same assembly ID. This can happen when an
    # assembly is constructed by applying different sets of transforms to
    # different sets of chain IDs. Each of these would have its own row.
    # Here we group rows by their assembly_id.
    pdbx_struct_assembly_gen = {}
    for assembly_id, oper_expression, asym_id_list in zip(
        cif['_pdbx_struct_assembly_gen.assembly_id'],
        cif['_pdbx_struct_assembly_gen.oper_expression'],
        cif['_pdbx_struct_assembly_gen.asym_id_list'],
    ):
      pdbx_struct_assembly_gen.setdefault(assembly_id, []).append({
          '_pdbx_struct_assembly_gen.assembly_id': assembly_id,
          '_pdbx_struct_assembly_gen.oper_expression': oper_expression,
          '_pdbx_struct_assembly_gen.asym_id_list': asym_id_list,
      })

    # We provide these separately to keep track of the original order that they
    # appear in the mmCIF.
    assembly_ids = cif['_pdbx_struct_assembly.id']
    oper_ids = cif['_pdbx_struct_oper_list.id']
    return cls(
        pdbx_struct_assembly=pdbx_struct_assembly,
        pdbx_struct_assembly_gen=pdbx_struct_assembly_gen,
        pdbx_struct_oper_list=pdbx_struct_oper_list,
        assembly_ids=assembly_ids,
        oper_ids=oper_ids,
    )

  @property
  def assembly_ids(self) -> Sequence[str]:
    return self._assembly_ids

  def asym_id_by_assembly_chain_id(self, assembly_id: str) -> Mapping[str, str]:
    asym_id_by_assembly_chain_id = {}
    for transform in self.get_transforms(assembly_id):
      for asym_id, assembly_chain_id in transform.chain_id_rename_map.items():
        asym_id_by_assembly_chain_id[assembly_chain_id] = asym_id
    return asym_id_by_assembly_chain_id

  def assembly_chain_ids_by_asym_id(
      self, assembly_id: str
  ) -> Mapping[str, set[str]]:
    assembly_chain_ids_by_asym_id = {}
    for transform in self.get_transforms(assembly_id):
      for asym_id, assembly_chain_id in transform.chain_id_rename_map.items():
        assembly_chain_ids_by_asym_id.setdefault(asym_id, set()).add(
            assembly_chain_id
        )
    return assembly_chain_ids_by_asym_id

  def get_default_assembly_id(self) -> str:
    """Gets a default assembly ID."""
    # The first assembly is usually (though not always) the best choice.
    # If we find a better heuristic for picking bioassemblies then this
    # method should be updated.
    return min(self._assembly_ids)

  def get_assembly_info(self, assembly_id: str) -> Mapping[str, str]:
    return {
        k.replace('_pdbx_struct_assembly.', ''): v
        for k, v in self._pdbx_struct_assembly[assembly_id].items()
    }

  def get_transforms(self, assembly_id: str) -> Sequence[Transform]:
    """Returns the transforms required to generate the given assembly."""
    partial_transforms = []
    all_chain_ids = set()
    for row in self._pdbx_struct_assembly_gen[assembly_id]:
      oper_expression = row['_pdbx_struct_assembly_gen.oper_expression']
      parsed_oper_id_seqs = mmcif.parse_oper_expr(oper_expression)
      label_asym_ids = row['_pdbx_struct_assembly_gen.asym_id_list'].split(',')
      all_chain_ids |= set(label_asym_ids)
      for parsed_oper_id_seq in parsed_oper_id_seqs:
        partial_transforms.append((parsed_oper_id_seq, label_asym_ids))

    # We start assigning new chain IDs by finding the largest chain ID in
    # the original structure that is involved in this bioassembly, and then
    # starting from the next one.
    max_int_chain_id = max(mmcif.str_id_to_int_id(c) for c in all_chain_ids)
    next_int_chain_id = max_int_chain_id + 1

    transforms = []
    has_been_renamed = set()
    for parsed_oper_id_seq, label_asym_ids in partial_transforms:
      chain_id_rename_map = {}
      for label_asym_id in label_asym_ids:
        if label_asym_id not in has_been_renamed:
          # The first time we see a label_asym_id we don't need to rename it.
          # This isn't strictly necessary since we don't provide any
          # guarantees about chain naming after bioassembly extraction but
          # can make it a bit easier to inspect and compare structures
          # pre and post bioassembly extraction.
          chain_id_rename_map[label_asym_id] = label_asym_id
          has_been_renamed.add(label_asym_id)
        else:
          chain_id_rename_map[label_asym_id] = mmcif.int_id_to_str_id(
              next_int_chain_id
          )
          next_int_chain_id += 1
      transforms.append(
          Transform(
              operations=[
                  self._operations[oper_id] for oper_id in parsed_oper_id_seq
              ],
              chain_ids=label_asym_ids,
              chain_id_rename_map=chain_id_rename_map,
          )
      )
    return transforms

  def to_mmcif_dict(self) -> Mapping[str, Sequence[str]]:
    """Returns the bioassembly data as a dict suitable for `mmcif.Mmcif`."""
    mmcif_dict = {}
    for assembly_id in self._assembly_ids:
      for column, val in self._pdbx_struct_assembly[assembly_id].items():
        mmcif_dict.setdefault(column, []).append(val)
      for row in self._pdbx_struct_assembly_gen[assembly_id]:
        for column, val in row.items():
          mmcif_dict.setdefault(column, []).append(val)
    for oper_id in self._oper_ids:
      for column, val in self._pdbx_struct_oper_list[oper_id].items():
        mmcif_dict.setdefault(column, []).append(val)
    return mmcif_dict

  def rename_label_asym_ids(
      self,
      mapping: Mapping[str, str],
      present_chains: set[str],
  ) -> Self:
    """Returns a new BioassemblyData with renamed label_asym_ids.

    Args:
      mapping: A mapping from original label_asym_ids to their new values. Any
        label_asym_ids in this BioassemblyData that are not in this mapping will
        remain unchanged.
      present_chains: A set of label_asym_ids that are actually present in the
        atom site list. All label_asym_ids that are in the BioassemblyData but
        not in present_chains won't be included in the output BioassemblyData.

    Returns:
      A new BioassemblyData with renamed label_asym_ids.

    Raises:
      ValueError: If any two previously distinct chains do not have unique names
          anymore after the rename.
    """
    new_pdbx_struct_assembly_gen = copy.deepcopy(self._pdbx_struct_assembly_gen)
    for rows in new_pdbx_struct_assembly_gen.values():
      for row in rows:
        old_asym_ids = row['_pdbx_struct_assembly_gen.asym_id_list'].split(',')
        new_asym_ids = [
            mapping.get(label_asym_id, label_asym_id)
            for label_asym_id in old_asym_ids
            if label_asym_id in present_chains
        ]
        if len(set(old_asym_ids) & present_chains) != len(set(new_asym_ids)):
          raise ValueError(
              'Can not rename chains, the new names are not unique: '
              f'{sorted(new_asym_ids)}.'
          )
        row['_pdbx_struct_assembly_gen.asym_id_list'] = ','.join(new_asym_ids)  # pytype: disable=unsupported-operands

    return BioassemblyData(
        pdbx_struct_assembly=copy.deepcopy(self._pdbx_struct_assembly),
        pdbx_struct_assembly_gen=new_pdbx_struct_assembly_gen,
        pdbx_struct_oper_list=copy.deepcopy(self._pdbx_struct_oper_list),
        assembly_ids=copy.deepcopy(self._assembly_ids),
        oper_ids=copy.deepcopy(self._oper_ids),
    )
