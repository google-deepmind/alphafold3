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

"""Table module for atom/residue/chain tables in Structure.

Tables are intended to be lightweight collections of columns, loosely based
on a pandas dataframe, for use in the Structure class.
"""

import abc
from collections.abc import Callable, Collection, Iterable, Iterator, Mapping, Sequence
import dataclasses
import functools
import graphlib
import typing
from typing import Any, Protocol, Self, TypeAlias, TypeVar, overload

from alphafold3.cpp import indexing
from alphafold3.cpp import membership
from alphafold3.cpp import string_array
import numpy as np


TableEntry: TypeAlias = str | int | float | None
FilterPredicate: TypeAlias = (
    TableEntry
    | Iterable[Any]  # Workaround for b/326384670. Tighten once fixed.
    | Callable[[Any], bool]  # Workaround for b/326384670. Tighten once fixed.
    | Callable[[np.ndarray], bool | np.ndarray]
)


class RowLookup(Protocol):

  def get_row_by_key(
      self,
      key: int,
      column_name_map: Mapping[str, str] | None = None,
  ) -> Mapping[str, Any]:
    ...


def _as_bool_mask(mask_or_scalar: Any, *, size: int) -> np.ndarray:
  """Coerces a scalar or 1D predicate result into a `(size,)` bool mask."""
  mask = np.asarray(mask_or_scalar, dtype=bool)
  return mask if mask.shape == (size,) else np.broadcast_to(mask, (size,))


def _eval_callable_predicate(
    col_arr: np.ndarray,
    predicate: (
        Callable[[Any], bool] | Callable[[np.ndarray], bool | np.ndarray]
    ),
    *,
    apply_per_element: bool,
) -> np.ndarray:
  """Evaluates a per-element or vectorized callable predicate into a mask."""
  size = col_arr.shape[0]
  if apply_per_element:
    # Cast via `bool(...)` so non-bool truthy/sequence returns match `np.array`.
    return np.fromiter(
        (bool(predicate(elem)) for elem in col_arr), dtype=bool, count=size
    )
  pred_fn = typing.cast(Callable[[np.ndarray], bool | np.ndarray], predicate)
  return _as_bool_mask(pred_fn(col_arr), size=size)


def _validate_filter_values(
    col_dtype: np.dtype, target_vals: Sequence[Any], *, col_name: str
) -> None:
  """Raises ValueError if target values are incompatible with `col_dtype`."""
  if col_dtype == object or isinstance(col_dtype, np.dtypes.StringDType):
    expected_dtype = None
  elif np.issubdtype(col_dtype, np.number):
    expected_dtype = np.number
  elif np.issubdtype(col_dtype, np.character):
    expected_dtype = np.character
  else:
    return
  for target_val in target_vals:
    if not (
        isinstance(target_val, str) or target_val is None
        if expected_dtype is None
        else np.issubdtype(np.asarray(target_val).dtype, expected_dtype)
    ):
      raise ValueError(
          f'Cannot filter column {col_name} (with type {col_dtype}) by value'
          f' {target_val} (with type {type(target_val)}).'
      )


def _eval_value_predicate(
    col_arr: np.ndarray, predicate: Any, *, col_name: str
) -> np.ndarray:
  """Evaluates a scalar or multi-value equality predicate into a bool mask."""
  # Strings are Iterable, so treat `str` predicates as a single scalar value.
  if isinstance(predicate, Iterable) and not isinstance(predicate, str):
    target_vals = tuple(predicate)
  else:
    target_vals = (predicate,)
  col_dtype = col_arr.dtype
  _validate_filter_values(col_dtype, target_vals, col_name=col_name)

  size = col_arr.shape[0]
  if not target_vals or size == 0:
    return np.zeros((size,), dtype=bool)
  if col_dtype == object:
    return membership.isin(col_arr, set(target_vals))
  if len(target_vals) == 1:
    return _as_bool_mask(col_arr == target_vals[0], size=size)
  if col_dtype in (np.int64, np.int32) and all(
      isinstance(v, (int, np.integer)) for v in target_vals
  ):
    return membership.isin(col_arr, set(target_vals))

  field_mask = np.zeros((size,), dtype=bool)
  for target_val in target_vals:
    np.logical_or(field_mask, col_arr == target_val, out=field_mask)
  return field_mask


def _take_rows(col_arr: np.ndarray, selector: np.ndarray | slice) -> np.ndarray:
  """Selects rows along the last axis; 1D columns skip NumPy's N-D index path."""
  return col_arr[selector] if col_arr.ndim == 1 else col_arr[..., selector]


@dataclasses.dataclass(frozen=True, kw_only=True)
class Table:
  """Parent class for structure tables.

  A table is a collection of columns of equal length, where one column is the
  key. The key uniquely identifies each row in the table.

  A table can refer to other tables by including a foreign key column, whose
  values are key values from the other table's key column. These column can have
  arbitrary names and are treated like any other integer-valued column.

  See the `Database` class in this module for utilities for handing sets of
  tables that are related via foreign keys.

  NB: This does not correspond to an mmCIF table.
  """

  key: np.ndarray

  def __post_init__(self):
    for col_name in self.columns:
      if (col_len := self.get_column(col_name).shape[-1]) != self.size:
        raise ValueError(
            f'All columns should have length {self.size} but got "{col_name}"'
            f' with length {col_len}.'
        )
      self.get_column(col_name).flags.writeable = False  # Make col immutable.
    if self.key.size and self.key.min() < 0:
      raise ValueError(
          'Key values must be non-negative. Got negative values:'
          f' {set(self.key[self.key < 0])}'
      )
    self.key.flags.writeable = False  # Make key immutable.

  @classmethod
  def _from_validated_columns(
      cls, columns: Iterable[tuple[str, np.ndarray]]
  ) -> Self:
    """Constructs a Table from row slices of an already-validated Table.

    Row-slicing an already-validated Table preserves column names, equal row
    counts, non-negative keys, and finite values, so only immutability needs to
    be re-applied on the newly allocated column arrays.

    Args:
      columns: Iterable of `(column_name, sliced_array)` pairs.

    Returns:
      A new immutable Table instance of `cls`.
    """
    instance = cls.__new__(cls)
    for col_name, col in columns:
      col.flags.writeable = False
      object.__setattr__(instance, col_name, col)
    return instance

  def _take(self, selector: np.ndarray) -> Self:
    """Returns the rows picked by an index array or boolean mask."""
    return self._from_validated_columns(
        (col, _take_rows(col_arr, selector))
        for col, col_arr in self.items.items()
    )

  def __getstate__(self) -> dict[str, Any]:
    """Returns members with cached properties removed for pickling."""
    cached_props = {
        k
        for k, v in self.__class__.__dict__.items()
        if isinstance(v, functools.cached_property)
    }
    return {k: v for k, v in self.__dict__.items() if k not in cached_props}

  @functools.cached_property
  def index_by_key(self) -> np.ndarray:
    """Mapping from key values to their index in the column arrays.

    i.e.: self.key[index_by_key[k]] == k
    """
    if not self.key.size:
      return np.array([], dtype=np.int64)
    else:
      index_by_key = np.zeros(np.max(self.key) + 1, dtype=np.int64)
      index_by_key[self.key] = np.arange(self.size)
      return index_by_key

  @functools.cached_property
  def columns(self) -> tuple[str, ...]:
    """The names of the columns in the table, including the key column."""
    return tuple(field.name for field in dataclasses.fields(self))

  @functools.cached_property
  def items(self) -> Mapping[str, np.ndarray]:
    """Returns the mapping from column names to column values."""
    return {col: getattr(self, col) for col in self.columns}

  @functools.cached_property
  def size(self) -> int:
    """The number of rows in the table."""
    return self.key.shape[-1]

  def __len__(self) -> int:
    return self.size

  def get_column(self, column_name: str) -> np.ndarray:
    """Gets a column by name."""
    # Performance optimisation: use the cached columns, instead of getattr.
    return self.items[column_name]

  def apply_array(self, arr: np.ndarray) -> Self:
    """Returns a sliced table using a key (!= index) array or a boolean mask."""
    if arr.dtype != bool:
      return self._take(indexing.index_of(arr, self.key))
    if arr.shape != (self.size,):
      return self._take(arr)
    if (indices := indexing.nonzero_indices(arr)) is None:
      return self  # Shortcut: No-op, so just return.
    return self._take(indices)

  def apply_index(self, index_arr: np.ndarray) -> Self:
    """Returns a sliced table using an index (!= key) array."""
    if index_arr.dtype == bool:
      raise ValueError('The index array must not be a boolean mask.')

    if self.size == 0 and index_arr.size == 0:
      return self  # Shortcut: No-op, so just return.

    if (
        # Perform cheap checks first.
        len(index_arr) == self.size
        and index_arr[0] == 0
        and index_arr[-1] == self.size - 1
        # Then do the full expensive check.
        and np.array_equal(index_arr, np.arange(self.size, dtype=np.int64))
    ):
      return self  # Shortcut: No-op, so just return.

    return self._take(index_arr)

  def slice_range(self, *, start: int, end: int) -> Self:
    """Returns a Table sliced to the contiguous `[start, end)` row range.

    Proper slices copy column arrays so small slices never retain NumPy views
    pinning the parent Table's memory.

    Args:
      start: Inclusive start row index (`0 <= start <= end`).
      end: Exclusive end row index (`start <= end <= self.size`).

    Returns:
      `self` when `[start, end)` covers the entire table; otherwise a new
      immutable `Table` containing rows `[start, end)`.

    Raises:
      ValueError: If `start` or `end` is out of `[0, self.size]` or
        `start > end`.
    """
    if not 0 <= start <= end <= self.size:
      raise ValueError(
          f'Invalid slice range [{start}, {end}) for table of size {self.size}.'
      )
    if start == 0 and end == self.size:
      return self
    return self._from_validated_columns(
        (col, (arr[start:end] if arr.ndim == 1 else arr[..., start:end]).copy())
        for col, arr in self.items.items()
    )

  def slice_spans(self, starts: np.ndarray, ends: np.ndarray) -> Self:
    """Returns a Table sliced to non-overlapping sorted `[start, end)` spans."""
    spans: list[slice] = []
    for start, end in zip(starts.tolist(), ends.tolist(), strict=True):
      if start >= end:
        continue
      if spans and spans[-1].stop == start:
        spans[-1] = slice(spans[-1].start, end)
      else:
        spans.append(slice(start, end))

    if not spans:
      return self.slice_range(start=0, end=0)
    if len(spans) == 1:
      return self.slice_range(start=spans[0].start, end=spans[0].stop)
    return self._from_validated_columns(
        (col, np.concatenate([_take_rows(arr, s) for s in spans], axis=-1))
        for col, arr in self.items.items()
    )

  def apply_array_to_column(
      self,
      column_name: str,
      arr: np.ndarray,
  ) -> np.ndarray:
    """Returns a sliced column array using a key array or a boolean mask."""
    col_arr = self.get_column(column_name)
    if arr.dtype == bool:
      return _take_rows(col_arr, arr)
    # Fast path: gather 1D native-endian aligned C-contiguous POD/object columns
    # directly in C++ without allocating an intermediate row-index array.
    if (
        col_arr.ndim == 1
        and col_arr.flags.c_contiguous
        and col_arr.flags.aligned
        and col_arr.dtype.isnative
        and col_arr.dtype.kind in 'biufO'
        and col_arr.dtype.itemsize <= 8
    ):
      return indexing.take_by_key(col_arr, arr, self.key)
    return _take_rows(col_arr, indexing.index_of(arr, self.key))

  def get_value_by_index(self, column_name: str, index: int) -> Any:
    return self.get_column(column_name)[index]

  def get_value_by_key(
      self,
      column_name: str,
      key: int | np.integer,
  ) -> TableEntry:
    """Gets the value of a column at the row with specified key value."""
    return self.get_value_by_index(column_name, self.index_by_key[key])

  @overload
  def __getitem__(self, key: str) -> np.ndarray:
    ...

  @overload
  def __getitem__(self, key: np.ndarray) -> 'Table':
    ...

  @overload
  def __getitem__(self, key: tuple[str, int | np.integer]) -> TableEntry:
    ...

  @overload
  def __getitem__(self, key: tuple[str, np.ndarray]) -> np.ndarray:
    ...

  def __getitem__(self, key):
    match key:
      case str():
        return self.get_column(key)
      case np.ndarray() as key_arr_or_mask:
        return self.apply_array(key_arr_or_mask)
      case str() as col, int() | np.integer() as key_val:
        return self.get_value_by_key(col, key_val)
      case str() as col, np.ndarray() as key_arr_or_mask:
        return self.apply_array_to_column(col, key_arr_or_mask)
      case _:
        if isinstance(key, tuple):
          err_msg = f'{key}, type: tuple({[type(v) for v in key]})'
        else:
          err_msg = f'{key}, type: {type(key)}'
        raise KeyError(err_msg)

  def get_row_by_key(
      self,
      key: int,
      column_name_map: Mapping[str, str] | None = None,
  ) -> dict[str, Any]:
    """Gets the row with specified key value."""
    return self.get_row_by_index(
        self.index_by_key[key], column_name_map=column_name_map
    )

  def get_row_by_index(
      self,
      index: int,
      column_name_map: Mapping[str, str] | None = None,
  ) -> dict[str, Any]:
    """Gets the row at the specified index."""
    if column_name_map is not None:
      return {
          renamed_col: self.get_value_by_index(col, index)
          for renamed_col, col in column_name_map.items()
      }
    else:
      return {col: self.get_value_by_index(col, index) for col in self.columns}

  def iterrows(
      self,
      *,
      row_keys: np.ndarray | None = None,
      column_name_map: Mapping[str, str] | None = None,
      **table_by_foreign_key_col: RowLookup,
  ) -> Iterator[Mapping[str, Any]]:
    """Yields rows from the table.

    This can be used to easily convert a table to a Pandas dataframe:

    ```py
    df = pd.DataFrame(table.iterrows())
    ```

    Args:
      row_keys: An optional array of keys of rows to yield. If None, all rows
        will be yielded.
      column_name_map: An optional mapping from desired keys in the row dicts to
        the names of the columns they correspond to.
      **table_by_foreign_key_col: An optional mapping from column names in this
        table, which are expected to be columns of foreign keys, to the table
        that the foreign keys point into. If provided, then the yielded rows
        will include data from the foreign tables at the appropriate key.
    """
    for i in (
        range(self.size)
        if row_keys is None
        else indexing.index_of(row_keys, self.key)
    ):
      row = self.get_row_by_index(i, column_name_map=column_name_map)
      for key_col, table in table_by_foreign_key_col.items():
        row.update(table.get_row_by_key(self[key_col][i]))
      yield row

  def with_column_names(
      self, column_name_map: Mapping[str, str]
  ) -> 'RenamedTableView':
    """Returns a view of this table with mapped column names."""
    return RenamedTableView(self, column_name_map=column_name_map)

  def make_filter_mask(
      self,
      mask: np.ndarray | None = None,
      *,
      apply_per_element: bool = False,
      **predicate_by_col: FilterPredicate,
  ) -> np.ndarray | None:
    """Returns a boolean array of rows to keep, or None if all can be kept.

    Args:
      mask: See `Table.filter`.
      apply_per_element: See `Table.filter`.
      **predicate_by_col: See `Table.filter`.

    Returns:
      Either a boolean NumPy array of length `(self.size,)` denoting which rows
      should be kept according to the input mask and predicates, or None. None
      implies there is no filtering required, and is used where possible
      instead of an all-True array to save time and space.
    """
    if mask is not None:
      if mask.shape != (self.size,):
        raise ValueError(
            f'mask must have shape ({self.size},). Got: {mask.shape}.'
        )
      if mask.dtype != bool:
        raise ValueError(f'mask must have dtype bool. Got: {mask.dtype}.')
    if not predicate_by_col:
      return mask

    result_mask: np.ndarray | None = mask
    for col, predicate in predicate_by_col.items():
      col_arr = self.get_column(col)
      if col_arr.ndim > 1:
        raise ValueError(
            f'Cannot filter by column {col} with more than 1 dimension.'
        )
      if callable(predicate):
        field_mask = _eval_callable_predicate(
            col_arr, predicate, apply_per_element=apply_per_element
        )
      else:
        field_mask = _eval_value_predicate(col_arr, predicate, col_name=col)
      # Allocate `&` rather than mutating caller `mask` or callable views.
      result_mask = (
          field_mask if result_mask is None else result_mask & field_mask
      )
    return result_mask

  def filter(
      self,
      mask: np.ndarray | None = None,
      *,
      apply_per_element: bool = False,
      invert: bool = False,
      **predicate_by_col: FilterPredicate,
  ) -> Self:
    """Filters the table using mask and/or predicates and returns a new table.

    Predicates can be either:
      1. A constant value, e.g. `'CA'`. In this case then only rows that match
        this value for the given column are retained.
      2. A (non-string) iterable e.g. `('A', 'B')`. In this
        case then rows are retained if they match any of the provided values for
        the given column.
      3. A boolean function e.g. `lambda b_fac: b_fac < 100.0`.
        In this case then only rows that evaluate to `True` are retained. By
        default this function's parameter is expected to be an array, unless
        `apply_per_element=True`.

    Args:
      mask: An optional boolean NumPy array with length equal to the table size.
        If provided then this will be combined with the other predicates so that
        a row is included if it is masked-in *and* matches all the predicates.
      apply_per_element: Whether apply predicates to each element in the column
        individually, or to pass the whole column array to the predicate.
      invert: If True then the returned table will contain exactly those rows
        that would be removed if this was `False`.
      **predicate_by_col: A mapping from column name to a predicate. Filtered
        columns must be 1D arrays. If multiple columns are provided as keyword
        arguments then each predicate is applied and the results are combined
        using a boolean AND operation, so an atom is only retained if it passes
        all predicates.

    Returns:
      A new table with the desired rows retained (or filtered out if
      `invert=True`).

    Raises:
      ValueError: If mask is provided and is not a bool array with shape
        `(num_atoms,)`.
    """
    filter_mask = self.make_filter_mask(
        mask, apply_per_element=apply_per_element, **predicate_by_col
    )
    if filter_mask is None:
      # No mask or predicate was specified, so we can return early.
      if not invert:
        return self
      else:
        return self[np.array((), dtype=np.int64)]  # pyrefly: ignore[bad-return]
    else:
      return self[~filter_mask if invert else filter_mask]  # pyrefly: ignore[bad-return]

  def _validate_keys_are_column_names(self, keys: Collection[str]) -> None:
    """Raises an error if any of the keys are not column names."""
    if mismatches := set(keys) - set(self.columns):
      raise ValueError(f'Invalid column names: {sorted(mismatches)}.')

  def copy_and_update(self, **new_column_by_column_name: np.ndarray) -> Self:
    """Returns a copy of this table with the specified changes applied.

    Args:
      **new_column_by_column_name: New values for the specified columns.

    Raises:
      ValueError: If a specified column name is not a column in this table.
    """
    self._validate_keys_are_column_names(new_column_by_column_name)
    return dataclasses.replace(self, **new_column_by_column_name)

  def copy_and_remap(
      self, **mapping_by_col: Mapping[TableEntry, TableEntry]
  ) -> Self:
    """Returns a copy of the table with the specified columns remapped.

    Args:
      **mapping_by_col: Each kwarg key should be the name of one of this table's
        columns, and each value should be a mapping. The values in the column
        will be looked up in the mapping and replaced with the result if one is
        found.

    Raises:
      ValueError: If a specified column name is not a column in this table.
    """
    self._validate_keys_are_column_names(mapping_by_col)
    if not self.size:
      return self
    remapped_cols = {}
    for column_name, mapping in mapping_by_col.items():
      col_arr = self.get_column(column_name)
      if col_arr.dtype == object:
        remapped = string_array.remap(col_arr, mapping)
      else:
        remapped = np.vectorize(lambda x: mapping.get(x, x))(col_arr)  # pylint: disable=cell-var-from-loop
      remapped_cols[column_name] = remapped
    return self.copy_and_update(**remapped_cols)


class RenamedTableView:
  """View of a table with renamed column names."""

  def __init__(self, table: Table, column_name_map: Mapping[str, str]):
    self._table = table
    self._column_name_map = column_name_map

  def get_row_by_key(
      self,
      key: int,
      column_name_map: Mapping[str, str] | None = None,
  ) -> Mapping[str, Any]:
    del column_name_map
    return self._table.get_row_by_key(
        key, column_name_map=self._column_name_map
    )


_DatabaseT = TypeVar('_DatabaseT', bound='Database')


class Database(abc.ABC):
  """Relational database base class."""

  @property
  @abc.abstractmethod
  def tables(self) -> Collection[str]:
    """The names of the tables in this database."""

  @abc.abstractmethod
  def get_table(self, table_name: str) -> Table:
    """Gets the table with the given name."""

  @property
  @abc.abstractmethod
  def foreign_keys(self) -> Mapping[str, Collection[tuple[str, str]]]:
    """Describes the relationship between keys in the database.

    Returns:
      A map from table names to pairs of `(column_name, foreign_table_name)`
      where `column_name` is a column containing foreign keys in the table named
      by the key, and the `foreign_table_name` is the name of the table that
      those foreign keys refer to.
    """

  @abc.abstractmethod
  def copy_and_update(
      self: _DatabaseT,
      **new_field_by_field_name,
  ) -> _DatabaseT:
    """Returns a copy of this database with the specified changes applied."""


def table_dependency_order(db: Database) -> Iterable[str]:
  """Yields the names of the tables in the database in dependency order.

  This order guarantees that a table appears after all other tables that
  it refers to using foreign keys. Specifically A < B implies that A contains
  no column that refers to B.key as a foreign key.

  Args:
    db: The database that defines the table names and foreign keys.
  """
  connections: dict[str, set[str]] = {}
  for table_name in db.tables:
    connection_set = set()
    for _, foreign_table in db.foreign_keys.get(table_name, ()):
      connection_set.add(foreign_table)
    connections[table_name] = connection_set
  yield from graphlib.TopologicalSorter(connections).static_order()


def concat_databases(dbs: Sequence[_DatabaseT]) -> _DatabaseT:
  """Concatenates the tables across a sequence of databases.

  Args:
    dbs: A non-empty sequence of database instances of the same type.

  Returns:
    A new database containing the concatenated tables from the input databases.

  Raises:
    ValueError: If `dbs` is empty or `dbs` contains different Database
      types.
  """
  if not dbs:
    raise ValueError('Need at least one value to concatenate.')
  distinct_db_types = {type(db) for db in dbs}
  if len(distinct_db_types) > 1:
    raise ValueError(
        f'All `dbs` must be of the same type, got: {distinct_db_types}'
    )

  first_db, *other_dbs = dbs
  concatted_tables: dict[str, Table] = {}
  key_offsets: dict[str, list[int]] = {}
  for table_name in table_dependency_order(first_db):
    first_table = first_db.get_table(table_name)
    columns: dict[str, list[np.ndarray]] = {
        column_name: [first_table.get_column(column_name)]
        for column_name in first_table.columns
    }
    key_offsets[table_name] = [
        first_table.key.max() + 1 if first_table.size else 0
    ]

    for prev_index, db in enumerate(other_dbs):
      table = db.get_table(table_name)
      for col_name in table.columns:
        columns[col_name].append(table.get_column(col_name))
      key_offset = key_offsets[table_name][prev_index]
      offset_key = table.key + key_offset
      columns['key'][-1] = offset_key
      if table.size:
        key_offsets[table_name].append(offset_key.max() + 1)
      else:
        key_offsets[table_name].append(key_offsets[table_name][prev_index])
      for fkey_col_name, foreign_table_name in first_db.foreign_keys.get(
          table_name, []
      ):
        fkey_columns = columns[fkey_col_name]
        fkey_columns[-1] = (
            fkey_columns[-1] + key_offsets[foreign_table_name][prev_index]
        )

    concatted_columns = {
        column_name: np.concatenate(values, axis=-1)
        for column_name, values in columns.items()
    }
    concatted_tables[table_name] = (type(first_table))(**concatted_columns)
  return first_db.copy_and_update(**concatted_tables)
