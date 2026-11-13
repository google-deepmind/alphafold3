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

import concurrent.futures
import pickle
import threading
from typing import Self

from absl.testing import absltest
from alphafold3.common import thread_safe


@thread_safe.cache_once
def _cached_module_sentinel() -> list[int]:
  return []


@thread_safe.thread_local_lru_cache(maxsize=4)
def _lru_cached_module_fn(value: int) -> list[int]:
  return [value]


class _CachedNode:

  def __init__(
      self,
      *,
      children: tuple[Self, ...] = (),
      executor: concurrent.futures.Executor | None = None,
  ):
    self.call_count = 0
    self._children = children
    self._executor = executor

  @thread_safe.CachedProperty
  def value(self) -> int:
    self.call_count += 1
    if self._executor is not None:
      return 1 + sum(self._executor.map(lambda c: c.value, self._children))
    return 1 + sum(c.value for c in self._children)


class ThreadSafeTest(absltest.TestCase):

  def test_cache_once_single_flight_concurrent_calls(self):
    num_threads = 10
    entered = threading.Event()
    release = threading.Event()
    call_count = 0
    expected: list[int] = []

    @thread_safe.cache_once
    def compute() -> list[int]:
      nonlocal call_count
      call_count += 1
      entered.set()
      self.assertTrue(release.wait(timeout=5.0))
      return expected

    with concurrent.futures.ThreadPoolExecutor(num_threads) as executor:
      futures = [executor.submit(compute) for _ in range(num_threads)]
      self.assertTrue(entered.wait(timeout=5.0))
      release.set()
      results = [future.result() for future in futures]

    self.assertEqual(call_count, 1)
    for result in results:
      self.assertIs(result, expected)

  def test_cache_once_caches_none(self):
    call_count = 0

    @thread_safe.cache_once
    def returns_none() -> None:
      nonlocal call_count
      call_count += 1
      return None

    self.assertIsNone(returns_none())
    self.assertIsNone(returns_none())
    self.assertEqual(call_count, 1)

  def test_cache_once_does_not_cache_exceptions(self):
    call_count = 0

    @thread_safe.cache_once
    def fail_once() -> str:
      nonlocal call_count
      call_count += 1
      if call_count == 1:
        raise ValueError('transient failure')
      return 'ok'

    with self.assertRaisesRegex(ValueError, 'transient failure'):
      fail_once()
    self.assertEqual(fail_once(), 'ok')
    self.assertEqual(fail_once(), 'ok')
    self.assertEqual(call_count, 2)

  def test_cache_once_pickles_module_function(self):
    unpickled = pickle.loads(pickle.dumps(_cached_module_sentinel))  # pylint: disable=g-unsafe-pickle-load
    self.assertIs(unpickled, _cached_module_sentinel)

  def test_thread_local_lru_cache_memoizes_within_thread_and_isolates_threads(
      self,
  ):
    num_threads = 8
    barrier = threading.Barrier(num_threads)
    call_count = 0
    count_lock = threading.Lock()

    @thread_safe.thread_local_lru_cache(maxsize=8)
    def compute(key: str) -> list[str]:
      nonlocal call_count
      with count_lock:
        call_count += 1
      return [key]

    def worker() -> list[str]:
      barrier.wait(timeout=5.0)
      first = compute('ALA')
      second = compute('ALA')
      self.assertIs(first, second)
      return first

    with concurrent.futures.ThreadPoolExecutor(num_threads) as executor:
      futures = [executor.submit(worker) for _ in range(num_threads)]
      results = [future.result(timeout=5.0) for future in futures]

    self.assertEqual(call_count, num_threads)
    self.assertLen({id(result) for result in results}, num_threads)

  def test_thread_local_lru_cache_evicts_least_recently_used(self):
    seen_keys: list[str] = []

    @thread_safe.thread_local_lru_cache(maxsize=2)
    def compute(key: str) -> str:
      seen_keys.append(key)
      return key.lower()

    self.assertEqual(compute('A'), 'a')
    self.assertEqual(compute('B'), 'b')
    self.assertEqual(compute('A'), 'a')
    self.assertEqual(compute('C'), 'c')
    self.assertEqual(compute('A'), 'a')
    self.assertEqual(compute('B'), 'b')
    self.assertEqual(seen_keys, ['A', 'B', 'C', 'B'])

  def test_thread_local_lru_cache_caches_none_and_kwargs(self):
    call_count = 0

    @thread_safe.thread_local_lru_cache(maxsize=8)
    def lookup(key: str, **kwargs: int) -> int | None:
      nonlocal call_count
      call_count += 1
      return None if key == 'missing' and not kwargs else call_count

    self.assertIsNone(lookup('missing'))
    self.assertIsNone(lookup('missing'))
    first = lookup('x', a=1)
    self.assertEqual(lookup('x', a=1), first)
    self.assertNotEqual(lookup('x', a=2), first)
    self.assertEqual(call_count, 3)

  def test_thread_local_lru_cache_does_not_cache_exceptions(self):
    call_count = 0

    @thread_safe.thread_local_lru_cache(maxsize=4)
    def flaky(key: str) -> str:
      nonlocal call_count
      call_count += 1
      if call_count == 1:
        raise ValueError('transient error')
      return f'{key}:{call_count}'

    with self.assertRaisesRegex(ValueError, 'transient error'):
      flaky('k')
    self.assertEqual(flaky('k'), 'k:2')
    self.assertEqual(flaky('k'), 'k:2')
    self.assertEqual(call_count, 2)

  def test_thread_local_lru_cache_recursion(self):
    call_count = 0

    @thread_safe.thread_local_lru_cache(maxsize=16)
    def fib(n: int) -> int:
      nonlocal call_count
      call_count += 1
      return n if n < 2 else fib(n - 1) + fib(n - 2)

    self.assertEqual(fib(6), 8)
    self.assertEqual(call_count, 7)

  def test_thread_local_lru_cache_rejects_non_positive_maxsize(self):
    for bad_maxsize in (0, -1):
      with self.assertRaisesRegex(ValueError, 'maxsize must be > 0'):
        thread_safe.thread_local_lru_cache(maxsize=bad_maxsize)

  def test_thread_local_lru_cache_pickles_module_function(self):
    unpickled = pickle.loads(pickle.dumps(_lru_cached_module_fn))  # pylint: disable=g-unsafe-pickle-load
    self.assertIs(unpickled, _lru_cached_module_fn)

  def test_cached_property_single_flight_concurrent_reads(self):
    num_threads = 10
    entered = threading.Event()
    release = threading.Event()
    call_count = 0
    expected: list[int] = []

    class Holder:

      @thread_safe.CachedProperty
      def items(self) -> list[int]:
        nonlocal call_count
        call_count += 1
        entered.set()
        if not release.wait(timeout=5.0):
          raise TimeoutError('Timed out waiting for release')
        return expected

    holder = Holder()
    with concurrent.futures.ThreadPoolExecutor(num_threads) as executor:
      futures = [
          executor.submit(lambda: holder.items) for _ in range(num_threads)
      ]
      self.assertTrue(entered.wait(timeout=5.0))
      release.set()
      results = [future.result(timeout=5.0) for future in futures]

    self.assertEqual(call_count, 1)
    for result in results:
      self.assertIs(result, expected)

  def test_cached_property_nested_properties_do_not_deadlock(self):
    with concurrent.futures.ThreadPoolExecutor(8) as executor:
      root = _CachedNode(
          children=tuple(_CachedNode() for _ in range(8)),
          executor=executor,
      )
      self.assertEqual(root.value, 9)

  def test_cached_property_caches_none(self):
    call_count = 0

    class Holder:

      @thread_safe.CachedProperty
      def maybe_value(self) -> None:
        nonlocal call_count
        call_count += 1
        return None

    holder = Holder()
    self.assertIsNone(holder.maybe_value)
    self.assertIsNone(holder.maybe_value)
    self.assertEqual(call_count, 1)

  def test_cached_property_pickles_instance_and_lock(self):
    node = _CachedNode()
    self.assertEqual(node.value, 1)
    unpickled = pickle.loads(pickle.dumps(node))  # pylint: disable=g-unsafe-pickle-load
    self.assertEqual(unpickled.value, 1)
    self.assertEqual(unpickled.call_count, 1)
    del unpickled.__dict__['value']
    self.assertEqual(unpickled.value, 1)
    self.assertEqual(unpickled.call_count, 2)

  def test_cached_property_class_access_returns_descriptor(self):
    self.assertIsInstance(_CachedNode.value, thread_safe.CachedProperty)


if __name__ == '__main__':
  absltest.main()
