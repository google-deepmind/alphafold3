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

from typing import Any, ClassVar, Iterable, Iterator, TypeVar, overload

import numpy as np

_T = TypeVar('_T')

class CifDict:
  class ItemView:
    def __iter__(self) -> Iterator[tuple[str, list[str]]]: ...
    def __len__(self) -> int: ...

  class KeyView:
    @overload
    def __contains__(self, key: str) -> bool: ...
    @overload
    def __contains__(self, key: object) -> bool: ...
    def __iter__(self) -> Iterator[str]: ...
    def __len__(self) -> int: ...

  class ValueView:
    def __iter__(self) -> Iterator[list[str]]: ...
    def __len__(self) -> int: ...

  def __init__(self, d: dict[str, Iterable[str]]) -> None: ...
  def copy_and_update(self, d: dict[str, Iterable[str]]) -> CifDict: ...
  def extract_loop_as_dict(self, prefix: str, index: str) -> dict:
    """Extracts loop associated with a prefix from mmCIF data as a dict.

    For instance for an mmCIF with these fields:
    '_a.ix': ['1', '2', '3']
    '_a.1': ['a.1.1', 'a.1.2', 'a.1.3']
    '_a.2': ['a.2.1', 'a.2.2', 'a.2.3']

    this function called with prefix='_a.', index='_a.ix' extracts:
    {'1': {'a.ix': '1', 'a.1': 'a.1.1', 'a.2': 'a.2.1'}
     '2': {'a.ix': '2', 'a.1': 'a.1.2', 'a.2': 'a.2.2'}
     '3': {'a.ix': '3', 'a.1': 'a.1.3', 'a.2': 'a.2.3'}}

    Args:
      prefix: Prefix shared by each of the data items in the loop. The prefix
        should include the trailing period.
      index: Which item of loop data should serve as the key.

    Returns:
      Dict of dicts; each dict represents 1 entry from an mmCIF loop,
      indexed by the index column.
    """

  def extract_loop_as_list(self, prefix: str) -> list:
    """Extracts loop associated with a prefix from mmCIF data as a list.

    Reference for loop_ in mmCIF:
    http://mmcif.wwpdb.org/docs/tutorials/mechanics/pdbx-mmcif-syntax.html

    For instance for an mmCIF with these fields:
    '_a.1': ['a.1.1', 'a.1.2', 'a.1.3']
    '_a.2': ['a.2.1', 'a.2.2', 'a.2.3']

    this function called with prefix='_a.' extracts:
    [{'_a.1': 'a.1.1', '_a.2': 'a.2.1'}
     {'_a.1': 'a.1.2', '_a.2': 'a.2.2'}
     {'_a.1': 'a.1.3', '_a.2': 'a.2.3'}]

    Args:
      prefix: Prefix shared by each of the data items in the loop. The prefix
        should include the trailing period.

    Returns:
      A list of dicts; each dict represents 1 entry from an mmCIF loop.
    """

  def get(self, key: str, default_value: _T = ...) -> list[str] | _T: ...
  def get_array(
      self, key: str, dtype: object = ..., gather: object = ...
  ) -> np.ndarray:
    """Returns values looked up in dict converted to a NumPy array.

    Args:
      key: Key in dictionary.
      dtype: Optional (default `object`) Specifies output dtype of array. One of
        [object, np.{int,uint}{8,16,32,64} np.float{32,64}]. As with NumPy use
        `object` to return a NumPy array of strings.
      gather: Optional one of [slice, np.{int,uint}{32,64}] non-intermediate
        version of get_array(key, dtype)[gather].

    Returns:
      A NumPy array of given dtype. An optimised equivalent to
      np.array(cif[key]).astype(dtype).  With support of '.' being treated
      as np.nan if dtype is one of np.float{32,64}.
      Identical strings will all reference the same object to save space.

    Raises:
      KeyError - if key is not found.
      TypeError - if dtype is not valid or supported.
      ValueError - if string cannot convert to dtype.
    """

  def get_data_name(self) -> str: ...
  def items(self) -> CifDict.ItemView: ...
  def keys(self) -> CifDict.KeyView: ...
  def to_string(self) -> str: ...
  def to_dict(self) -> dict[str, list[str]]: ...
  def value_length(self, key: str) -> int: ...
  def values(self) -> CifDict.ValueView: ...
  def __bool__(self) -> bool: ...
  def __contains__(self, key: str) -> bool: ...
  def __getitem__(self, key: str) -> list[str]: ...
  def __getstate__(self) -> tuple: ...
  def __iter__(self) -> Iterator[str]: ...
  def __len__(self) -> int: ...
  def __setstate__(self, state: tuple) -> None: ...

def tokenize(cif_string: str) -> list[str]: ...
def split_line(line: str) -> list[str]: ...
def from_string(mmcif_string: str | bytes) -> CifDict: ...
def parse_multi_data_cif(cif_string: str | bytes) -> dict[str, CifDict]: ...
