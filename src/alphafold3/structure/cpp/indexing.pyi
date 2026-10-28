# Copyright 2026 DeepMind Technologies Limited
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

import numpy as np

def index_of(keys: np.ndarray, table_keys: np.ndarray) -> np.ndarray[np.int64]:
  """Maps 1D signed-integer keys to their 0-based positions in 1D table_keys.

  Args:
    keys: 1D signed-integer NumPy array of query keys.
    table_keys: 1D signed-integer NumPy array of table keys (if duplicate keys
      exist, the last matching index is returned).

  Returns:
    A 1D int64 NumPy array with the same length as `keys` containing the 0-based
    index in `table_keys` for each element of `keys` (may alias `keys` when it
    is a read-only int64 array and `table_keys` is `0..m - 1`).

  Raises:
    IndexError: If `keys` or `table_keys` is not a 1D signed-integer array, or
      any element of `keys` is absent from `table_keys`.
  """
def nonzero_indices(mask: np.ndarray) -> np.ndarray[np.int64] | None:
  """Returns int64 indices of True elements in a 1D bool mask, or None if all True."""

def take_by_key(
    values: np.ndarray, keys: np.ndarray, table_keys: np.ndarray
) -> np.ndarray:
  """Gathers 1D values at the positions of 1D keys in 1D table_keys.

  Equivalent to `values[index_of(keys, table_keys)]`, fused to avoid
  materializing an intermediate index array.

  Args:
    values: 1D aligned C-contiguous native-endian NumPy array with the same
      length as `table_keys` and dtype `bool`, 1/2/4/8-byte integer or
      floating-point, or `object`.
    keys: 1D signed-integer NumPy array of query keys.
    table_keys: 1D signed-integer NumPy array of table keys (if duplicate keys
      exist, the last matching index is used).

  Returns:
    A new 1D NumPy array with the same length as `keys` and dtype as `values`.

  Raises:
    IndexError: If `keys` or `table_keys` is not a 1D signed-integer array, or
      any element of `keys` is absent from `table_keys`.
    ValueError: If `values` is not a 1D aligned C-contiguous array with the same
      length as `table_keys`.
    TypeError: If `values` has an unsupported or byte-swapped dtype.
  """
