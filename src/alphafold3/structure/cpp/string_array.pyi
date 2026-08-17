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

from collections.abc import Sequence
from typing import Any, overload

import numpy as np


def format_float_array(
    values: Sequence[float], num_decimal_places: int
) -> list[str]: ...


def isin(
    array: np.ndarray[object],
    test_elements: set[str | bytes],
    *,
    invert: bool = ...,
) -> np.ndarray[bool]: ...


@overload
def remap(
    array: np.ndarray[object],
    mapping: dict[str, str],
    default_value: str,
    inplace: bool = ...,
) -> np.ndarray[object]: ...


@overload
def remap(
    array: np.ndarray[object],
    mapping: dict[str, str],
    inplace: bool = ...,
) -> np.ndarray[object]: ...


def remap_multiple(
    arrays: Sequence[np.ndarray[object]],
    mapping: dict[tuple[Any], int],
) -> np.ndarray[int]: ...
