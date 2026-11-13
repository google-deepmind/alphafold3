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

from absl.testing import absltest
from alphafold3.common import thread_safe


@thread_safe.cache_once
def _cached_module_sentinel() -> list[int]:
  return []


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


if __name__ == '__main__':
  absltest.main()
