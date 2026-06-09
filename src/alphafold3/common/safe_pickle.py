# Copyright 2025 DeepMind Technologies Limited
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

"""Restricted-safe wrapper around pickle for loading trusted data.

This prevents arbitrary object instantiation during unpickling by only
allowing a small allowlist of built-in, innocuous types.

Intended for loading pickled constant data that ships with the repository.
If the pickle is tampered with, an UnpicklingError will be raised instead
of silently executing attacker-controlled bytecode.
"""

from collections.abc import Collection
import pickle
from typing import Any, BinaryIO, Final


# Builtin types expected from AlphaFold 3 generated data.
_ALLOWED_BUILTINS: Final[Collection[str]] = frozenset({
    "NoneType",
    "bool",
    "bytes",
    "dict",
    "float",
    "frozenset",
    "int",
    "list",
    "set",
    "str",
    "tuple",
})


class _RestrictedUnpickler(pickle.Unpickler):
  """A pickle `Unpickler` that forbids loading arbitrary global classes."""

  def find_class(self, module: str, name: str) -> Any:
    """Returns the class for `module` and `name` if allowed."""
    if module == "builtins" and name in _ALLOWED_BUILTINS:
      return super().find_class(module, name)
    raise pickle.UnpicklingError(f"Can't unpickle disallowed '{module}.{name}'")


def load(file_obj: BinaryIO) -> Any:
  """Safely loads pickle data from an already-opened binary file handle.

  Only built-in container/primitive types listed in `_ALLOWED_BUILTINS` are
  permitted. Any attempt to load other types raises `pickle.UnpicklingError`.

  Args:
    file_obj: A binary file-like object open for reading.

  Returns:
    The unpickled data.
  """

  return _RestrictedUnpickler(file_obj).load()
