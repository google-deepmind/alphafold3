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

from collections.abc import Set
from typing import Any
import numpy as np

def isin(
    array: np.ndarray[np.int64 | object],
    test_elements: Set[int | Any] | np.ndarray,
    *,
    invert: bool = ...,
) -> np.ndarray[bool]:
  """Computes whether each element is in test_elements.

  Same use as np.isin, but much faster. If len(array) = n, len(test_elements) =
  m:
  * This function has complexity O(n) for a set or O(n + m) for an array.
  * np.isin with kind='sort' has complexity O(m * log(m) + n * log(m)).

  Args:
    array: Input NumPy array. The type of the member elements must be compatible
      with the type of the elements in test_elements.
    test_elements: A Python set/frozenset, or a signed-integer NumPy array when
      `array` is a signed-integer array.
    invert: If True, the values in the returned array are inverted, as if
      calculating `element not in test_elements`. Default is False. `isin(a, b,
      invert=True)` is equivalent to but faster than `~isin(a, b)`.

  Returns
    A boolean array of the same shape as the input array. Each value `val` is:
    * `val in test_elements` if `invert=False`,
    * `val not in test_elements` if `invert=True`.
  """
