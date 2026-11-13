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
import contextlib
import functools
import threading
from typing import Any, Final, ParamSpec, Self, TypeVar, cast, overload

_P = ParamSpec('_P')
_T = TypeVar('_T')
_EMPTY: Final[Any] = object()
_CREATE_LOCK_LOCK: Final[threading.Lock] = threading.Lock()


class PicklableLock:
  """Mutex that can be pickled."""

  __slots__ = ('_lock',)

  def __init__(self):
    self._lock = threading.Lock()

  def __getstate__(self) -> dict[str, Any]:
    return {}

  def __setstate__(self, state: dict[str, Any]) -> None:
    del state
    self._lock = threading.Lock()

  def __enter__(self) -> None:
    self._lock.acquire()

  def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
    self._lock.release()


class CachedProperty(functools.cached_property[_T]):
  """Thread-safe functools.cached_property with per-instance locking."""

  __slots__ = ()

  def __init__(self, func: Callable[[Any], _T]):
    super().__init__(func)
    if hasattr(self, 'lock'):
      self.lock = contextlib.nullcontext()

  @overload
  def __get__(self, instance: None, owner: type[Any] | None = None) -> Self:
    ...

  @overload
  def __get__(self, instance: Any, owner: type[Any] | None = None) -> _T:
    ...

  def _get_or_create_lock(self, instance: Any) -> PicklableLock:
    lock_attr = f'{self.attrname}_lock'
    if (lock := instance.__dict__.get(lock_attr)) is None:
      with _CREATE_LOCK_LOCK:
        if (lock := instance.__dict__.get(lock_attr)) is None:
          lock = PicklableLock()
          instance.__dict__[lock_attr] = lock
    return lock

  def __get__(
      self, instance: Any = None, owner: type[Any] | None = None
  ) -> _T | Self:
    if instance is None:
      return self
    if (result := instance.__dict__.get(self.attrname, _EMPTY)) is not _EMPTY:  # pyrefly: ignore[no-matching-overload]
      return cast(_T, result)
    with self._get_or_create_lock(instance):
      if (result := instance.__dict__.get(self.attrname, _EMPTY)) is not _EMPTY:  # pyrefly: ignore[no-matching-overload]
        return cast(_T, result)
      return super().__get__(instance, owner)


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


def thread_local_lru_cache(
    *, maxsize: int
) -> Callable[[Callable[_P, _T]], Callable[_P, _T]]:
  """Caches a callable per thread with a bounded LRU cache.

  Args:
    maxsize: Maximum number of cached entries per thread.

  Returns:
    A decorator giving each calling thread its own independent
    `functools.lru_cache(maxsize=maxsize)` instance so cache hits avoid
    cross-thread lock and refcount contention under `PYTHON_GIL=0`.

  Raises:
    ValueError: If `maxsize <= 0`.
  """
  if maxsize <= 0:
    raise ValueError(f'maxsize must be > 0, got {maxsize}')

  def decorator(func: Callable[_P, _T]) -> Callable[_P, _T]:
    state = threading.local()

    @functools.wraps(func)
    def wrapper(*args: _P.args, **kwargs: _P.kwargs) -> _T:
      try:
        cached = state.cached
      except AttributeError:
        cached = state.cached = functools.lru_cache(maxsize=maxsize)(func)
      return cached(*args, **kwargs) if kwargs else cached(*args)

    return wrapper

  return decorator
