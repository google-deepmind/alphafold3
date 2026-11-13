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

"""Utilities for thread-safe code."""

from collections.abc import Callable
import functools
import threading
from typing import Any, TypeVar

_T = TypeVar('_T')
_EMPTY = object()


def cache_once(func: Callable[[], _T]) -> Callable[[], _T]:
  """Caches the return value of a zero-argument callable across threads.

  Args:
    func: Zero-argument callable initialized lazily on the first successful
      call. Concurrent cold callers wait for the in-flight call rather than
      invoking `func` redundantly; raised exceptions are not cached.

  Returns:
    A thread-safe zero-argument callable returning the cached result of `func`.
  """
  result: Any = _EMPTY
  lock = threading.Lock()

  @functools.cache
  @functools.wraps(func)
  def wrapper() -> _T:
    nonlocal result
    with lock:
      if result is _EMPTY:
        result = func()
      return result

  return wrapper
