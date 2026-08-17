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

"""Low level mmCIF parsing operations and wrappers for nicer C++/Py errors.

Note that the cif_dict.CifDict class has many useful methods to help with data
extraction which are not shown in this file. You can find them in cif_dict.clif
together with docstrings. The cif_dict.CifDict class behaves like an immutable
Python dictionary (some methods are not implemented though).
"""
from collections.abc import Callable, Mapping, Sequence
import functools
import itertools
import re
from typing import ParamSpec, TypeAlias, TypeVar

from alphafold3.constants import chemical_components
from alphafold3.cpp import cif_dict
from alphafold3.cpp import mmcif_atom_site
from alphafold3.cpp import mmcif_struct_conn
from alphafold3.cpp import string_array
import numpy as np

Mmcif = cif_dict.CifDict


_P = ParamSpec('_P')
_T = TypeVar('_T')
_WappedFn: TypeAlias = Callable[_P, _T]


@functools.lru_cache(maxsize=256)
def int_id_to_str_id(num: int) -> str:
  """Encodes a number as a string, using reverse spreadsheet style naming.

  Args:
    num: A positive integer.

  Returns:
    A string that encodes the positive integer using reverse spreadsheet style,
    naming e.g. 1 = A, 2 = B, ..., 27 = AA, 28 = BA, 29 = CA, ... This is the
    usual way to encode chain IDs in mmCIF files.
  """
  if num <= 0:
    raise ValueError(f'Only positive integers allowed, got {num}.')

  num = num - 1  # 1-based indexing.
  output = []
  while num >= 0:
    output.append(chr(num % 26 + ord('A')))
    num = num // 26 - 1
  return ''.join(output)


@functools.lru_cache(maxsize=256)
def str_id_to_int_id(str_id: str) -> int:
  """Encodes an mmCIF-style string chain ID as an integer.

  The integer IDs are one based so this function is the inverse of
  int_id_to_str_id.

  Args:
    str_id: A string chain ID consisting only of upper case letters A-Z.

  Returns:
    An integer that can be used to order mmCIF chain IDs in the standard
    (reverse spreadsheet style) ordering.
  """
  if not re.match('^[A-Z]+$', str_id):
    raise ValueError(f'String ID must be upper case letters, got {str_id}.')

  offset = ord('A') - 1
  output = 0
  for i, c in enumerate(str_id):
    output += (ord(c) - offset) * int(26**i)
  return output


def from_string(mmcif_string: str | bytes) -> Mmcif:
  return cif_dict.from_string(mmcif_string)


def parse_multi_data_cif(cif_string: str) -> dict[str, Mmcif]:
  """Parses a CIF string with multiple data records.

  For instance, the CIF string:

  ```
  data_001
  _foo bar
  #
  data_002
  _foo baz
  ```

  is parsed as:

  ```
  {'001': Mmcif({'_foo': ['bar']}), '002': Mmcif({'_foo': ['baz']})}
  ```

  Args:
    cif_string: The multi-data CIF string to be parsed.

  Returns:
    A dictionary mapping record names to Mmcif objects with data.
  """
  return cif_dict.parse_multi_data_cif(cif_string)


class BondParsingError(Exception):
  """Exception raised by errors when getting bond atom indices."""


def get_bond_atom_indices(
    mmcif: Mmcif,
    model_id: str = '1',
) -> tuple[Sequence[int], Sequence[int]]:
  """Extracts the indices of the atoms that participate in bonds.

  Args:
    mmcif: The mmCIF object to process.
    model_id: The ID of the model that the returned atoms will belong to. This
      should be a value in the mmCIF's _atom_site.pdbx_PDB_model_num column.

  Returns:
    Two lists of atom indices, `from_atoms` and `to_atoms`, each one having
    length num_bonds (as defined by _struct_conn, the bonds table). The bond
    i, defined by the i'th row in _struct_conn, is a bond from atom at index
    from_atoms[i], to the atom at index to_atoms[i]. The indices are simple
    0-based indexes into the columns of the _atom_site table in the input
    mmCIF, and do not necessarily correspond to the values in _atom_site.id,
    or any other column.

  Raises:
    BondParsingError: If any of the required tables or columns are not present
    in
      the mmCIF, or if the _struct_conn table refers to atoms that cannot
      be found in the _atom_site table.
  """
  try:
    return mmcif_struct_conn.get_bond_atom_indices(mmcif, model_id)
  except ValueError as e:
    raise BondParsingError(str(e)) from e


def get_or_infer_type_symbol(
    mmcif: Mmcif, ccd: chemical_components.Ccd | None = None
) -> Sequence[str]:
  """Returns the type symbol (element) for all of the atoms.

  Args:
    mmcif: A parsed mmCIF file in the Mmcif format.
    ccd: The chemical component dictionary. If not provided, defaults to the
      cached CCD.

  If present, returns the _atom_site.type_symbol. If not, infers it using
  _atom_site.label_comp_id (residue name), _atom_site.label_atom_id (atom name)
  and the CCD.
  """
  ccd = ccd or chemical_components.Ccd()
  type_symbol_fn = lambda res_name, atom_name: chemical_components.type_symbol(
      ccd, res_name, atom_name
  )
  return mmcif_atom_site.get_or_infer_type_symbol(mmcif, type_symbol_fn)


def get_chain_type_by_entity_id(mmcif: Mmcif) -> Mapping[str, str]:
  """Returns mapping from entity ID to its type or polymer type if available.

  If the entity is in the _entity_poly table, returns its polymer chain type.
  If not, returns the type as specified in the _entity table.

  Args:
    mmcif: CifDict holding the mmCIF.
  """
  poly_entity_id = mmcif.get('_entity_poly.entity_id', [])
  poly_type = mmcif.get('_entity_poly.type', [])
  poly_type_by_entity_id = dict(zip(poly_entity_id, poly_type, strict=True))

  chain_type_by_entity_id = {}
  for entity_id, entity_type in zip(
      mmcif.get('_entity.id', []), mmcif.get('_entity.type', []), strict=True
  ):
    chain_type = poly_type_by_entity_id.get(entity_id) or entity_type
    chain_type_by_entity_id[entity_id] = chain_type

  return chain_type_by_entity_id


def get_internal_to_author_chain_id_map(mmcif: Mmcif) -> Mapping[str, str]:
  """Returns a mapping from internal chain ID to the author chain ID.

  Note that this is not a bijection. One author chain ID can map to multiple
  internal chain IDs. For example, a protein chain and a ligand bound to it will
  share the same author chain ID, but they will each have a unique internal
  chain ID).

  Args:
    mmcif: CifDict holding the mmCIF.
  """
  return mmcif_atom_site.get_internal_to_author_chain_id_map(mmcif)


def get_experimental_method(mmcif: Mmcif) -> str | None:
  field = '_exptl.method'
  return ','.join(mmcif[field]).lower() if field in mmcif else None


def get_release_date(mmcif: Mmcif) -> str | None:
  """Returns the oldest revision date."""
  if '_pdbx_audit_revision_history.revision_date' not in mmcif:
    return None

  # Release dates are ISO-8601, hence sort well.
  return min(mmcif['_pdbx_audit_revision_history.revision_date'])


def get_resolution(mmcif: Mmcif) -> float | None:
  """Returns the resolution of the structure.

  More than one resolution can be reported in an mmCIF. This function returns
  the first one (in the order _refine.ls_d_res_high,
  _em_3d_reconstruction.resolution, _reflns.d_resolution_high) that appears
  in the mmCIF as is parseable as a float.

  Args:
    mmcif: An `Mmcif` object.

  Returns:
    The resolution as reported in the mmCIF.
  """
  for res_key in ('_refine.ls_d_res_high',
                  '_em_3d_reconstruction.resolution',
                  '_reflns.d_resolution_high'):
    if res_key in mmcif:
      try:
        raw_resolution = mmcif[res_key][0]
        return float(raw_resolution)
      except ValueError:
        continue
  return None


def parse_oper_expr(oper_expression: str) -> list[tuple[str, ...]]:
  """Determines which transforms to apply based on an MMCIF oper_expression str.

  Args:
    oper_expression: the field oper_expression from MMCIF format data.
      Transform ids may be either numbers or single letters. Hyphens are used to
      denote a numeric range of transforms to apply, and commas are used to
      delimit a sequence of transforms. Where two sets of parentheses are
      adjacent without a comma, the two sets of transforms should be combined as
      a cartesian product, i.e. all possible pairs.
      example 1,2,3 -> generate 3 copies of each chain by applying 1, 2 or 3.
      example (1-3) -> generate 3 copies of each chain by applying 1, 2 or 3.
      example (1-3)(4-6) -> generate 9 copies of each chain by applying one of
        [(1,4), (1,5), (1,6),
         (2,4), (2,5), (2,6),
         (3,4), (3,5), (3,6)]
      example (P) -> apply transform with id P.

  Raises:
    ValueError: Failure to parse oper_expression.

  Returns:
    A list with one element for each chain copy that should be generated.
    Each element is a list of transform ids to apply.
  """
  # Expand ranges, e.g. 1-4 -> 1,2,3,4.
  def range_expander(match):
    return ','.join(
        [str(i) for i in range(int(match.group(1)),
                               int(match.group(2)) + 1)])

  ranges_expanded = re.sub(r'\b(\d+)-(\d+)', range_expander, oper_expression)

  if re.fullmatch(r'(\w+,)*\w+', ranges_expanded):
    # No brackets, just a single range, e.g. "1,2,3".
    return [(t,) for t in ranges_expanded.split(',')]
  elif re.fullmatch(r'\((\w+,)*\w+\)', ranges_expanded):
    # Single range in brackets, e.g. "(1,2,3)".
    return [(t,) for t in ranges_expanded[1:-1].split(',')]
  elif re.fullmatch(r'\((\w+,)*\w+\)\((\w+,)*\w+\)', ranges_expanded):
    # Cartesian product of two ranges, e.g. "(1,2,3)(4,5)".
    part1, part2 = ranges_expanded[1:-1].split(')(')
    return list(itertools.product(part1.split(','), part2.split(',')))
  else:
    raise ValueError(f'Unsupported oper_expression format: {oper_expression}')


def format_float_array(
    values: np.ndarray, num_decimal_places: int) -> Sequence[str]:
  """Converts 1D array to a list of strings with the given number of decimals.

  This function is faster than converting via Python list comprehension, e.g.:
  atoms_x = ['%.3f' % x for x in atoms_x]

  Args:
    values: A numpy array with values to convert. This array is casted to
      float32 before doing the conversion.
    num_decimal_places: The number of decimal points to keep, including trailing
      zeros. E.g. for 1.07 and num_decimal_places=1: 1.1,
      num_decimal_places=2: 1.07, num_decimal_places=3: 1.070.

  Returns:
    A list of formatted strings.
  """
  if values.ndim != 1:
    raise ValueError(f'The given array must be 1D, got {values.ndim}D')

  return string_array.format_float_array(
      values=values.astype(np.float32), num_decimal_places=num_decimal_places
  )
