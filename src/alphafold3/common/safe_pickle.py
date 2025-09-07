"""Restricted-safe wrapper around pickle for loading trusted data.

This prevents arbitrary object instantiation during unpickling by only
allowing a small whitelist of built-in, innocuous types.

Intended for loading pickled constant data that ships with the repository.
If the pickle is tampered with, an UnpicklingError will be raised instead
of silently executing attacker-controlled bytecode.
"""

from __future__ import annotations

import pickle
from typing import Any, Set


_ALLOWED_BUILTINS: Set[str] = {
    # Core container and primitive types that appear in our constant pickles.
    "tuple",
    "list",
    "dict",
    "set",
    "frozenset",
    "str",
    "bytes",
    "int",
    "float",
    "bool",
    "NoneType",
}


class _RestrictedUnpickler(pickle.Unpickler):
    """A Pickle `Unpickler` that forbids loading arbitrary global classes."""

    def find_class(self, module: str, name: str) -> Any:  # noqa: D401, N802
        if module == "builtins" and name in _ALLOWED_BUILTINS:
            return super().find_class(module, name)
        raise pickle.UnpicklingError(
            f"Attempted to load disallowed global '{module}.{name}' via pickle"
        )


def safe_load(file_obj) -> Any:  # type: ignore[valid-type]
    """Safely loads pickle data from an already-opened binary file handle.

    Only built-in container/primitive types listed in `_ALLOWED_BUILTINS`
    are permitted.  Any attempt to load other types raises `pickle.UnpicklingError`.
    """

    return _RestrictedUnpickler(file_obj).load() 