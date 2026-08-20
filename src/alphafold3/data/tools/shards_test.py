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

"""Tests for database shard path handling."""

from absl.testing import absltest
from alphafold3.data.tools import shards
import os
import tempfile


class ShardsTest(absltest.TestCase):

  def test_parse_shard_spec_explicit_count(self):
    spec = shards.parse_shard_spec('/tmp/db@20')
    self.assertIsNotNone(spec)
    self.assertEqual(spec.prefix, '/tmp/db')
    self.assertEqual(spec.num_shards, 20)

  def test_parse_shard_spec_star_from_filesystem(self):
    with tempfile.TemporaryDirectory() as d:
      prefix = os.path.join(d, 'db')
      for i in range(5):
        with open(f'{prefix}-{i:05d}-of-00005', 'w'):
          pass
      spec = shards.parse_shard_spec(f'{prefix}@*')
      self.assertIsNotNone(spec)
      self.assertEqual(spec.num_shards, 5)

  def test_parse_shard_spec_star_glob_chars_in_prefix(self):
    with tempfile.TemporaryDirectory() as d:
      prefix = os.path.join(d, 'db[abc')
      for i in range(3):
        with open(f'{prefix}-{i:05d}-of-00003', 'w'):
          pass
      spec = shards.parse_shard_spec(f'{prefix}@*')
      self.assertIsNotNone(spec)
      self.assertEqual(spec.num_shards, 3)

  def test_parse_shard_spec_star_no_files(self):
    with tempfile.TemporaryDirectory() as d:
      self.assertIsNone(shards.parse_shard_spec(os.path.join(d, 'nope@*')))

  def test_get_sharded_paths(self):
    paths = shards.get_sharded_paths('/tmp/db@3')
    self.assertEqual(len(paths), 3)
    self.assertTrue(paths[0].endswith('-00000-of-00003'))


if __name__ == '__main__':
  absltest.main()