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

from collections.abc import Sequence
import os
import pathlib

from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.data.tools import shards


class ShardsTest(parameterized.TestCase):

  @parameterized.named_parameters(
      ('no_shard_spec', '/path/to/file'),
      ('negative_number_of_shards', '/path/to/file@-12'),
      ('letter_in_num_shard', '/path/to/file@a12'),
      ('num_shards_too_high', '/path/t@o/file@100000'),
  )
  def test_invalid_shard(self, shard_spec: str):
    self.assertIsNone(shards.get_sharded_paths(shard_spec))

  @parameterized.named_parameters(
      ('zero_shards', '/path/to/file@0', []),
      ('one_shard', '/path/to/file@1', ['/path/to/file-00000-of-00001']),
      (
          'two_shards',
          '/path/to/file@2',
          ['/path/to/file-00000-of-00002', '/path/to/file-00001-of-00002'],
      ),
      (
          'ten_shards',
          '/path/to/file@10',
          [
              '/path/to/file-00000-of-00010',
              '/path/to/file-00001-of-00010',
              '/path/to/file-00002-of-00010',
              '/path/to/file-00003-of-00010',
              '/path/to/file-00004-of-00010',
              '/path/to/file-00005-of-00010',
              '/path/to/file-00006-of-00010',
              '/path/to/file-00007-of-00010',
              '/path/to/file-00008-of-00010',
              '/path/to/file-00009-of-00010',
          ],
      ),
      (
          'with_suffix',
          '/path/to/file@1_test.txt',
          ['/path/to/file-00000-of-00001_test.txt'],
      ),
  )
  def test_shard_spec(self, shard_spec: str, expected_shards: Sequence[str]):
    self.assertEqual(shards.get_sharded_paths(shard_spec), expected_shards)

  @parameterized.named_parameters(
      (
          'zero_shards',
          '/path/to/file@0',
          shards.ShardSpec('/path/to/file', 0, ''),
      ),
      (
          'valid_shard',
          '/path/to/file@12345',
          shards.ShardSpec('/path/to/file', 12345, ''),
      ),
      (
          'spec_with_extension_suffix',
          '/path/to/file@20.sstable',
          shards.ShardSpec('/path/to/file', 20, '.sstable'),
      ),
      (
          'spec_with_other_suffix',
          '/path/to/file@20_keys.sstable',
          shards.ShardSpec('/path/to/file', 20, '_keys.sstable'),
      ),
      ('no_shard', '/path/to/file', None),
      ('invalid_spec', '/path/to/file@a.txt', None),
  )
  def test_parse_shard_spec(
      self, shard_spec: str, expected_parsed_shard: shards.ShardSpec | None
  ):
    self.assertEqual(
        shards.parse_shard_spec(shard_spec), expected_parsed_shard
    )

  @parameterized.named_parameters(
      (
          'with_suffix',
          'file-00000-of-00039_test.txt',
          'file@*_test.txt',
          shards.ShardSpec('file', 39, '_test.txt'),
      ),
      (
          'without_suffix',
          'file-00000-of-00039',
          'file@*',
          shards.ShardSpec('file', 39, ''),
      ),
      (
          'min_shard_without_suffix',
          'file-00000-of-00001',
          'file@*',
          shards.ShardSpec('file', 1, ''),
      ),
      (
          'min_shard_with_suffix',
          'file-00000-of-00001_test.txt',
          'file@*_test.txt',
          shards.ShardSpec('file', 1, '_test.txt'),
      ),
      (
          'max_shard_without_suffix',
          'file-00000-of-99999',
          'file@*',
          shards.ShardSpec('file', 99_999, ''),
      ),
      (
          'max_shard_with_suffix',
          'file-00000-of-99999_test.txt',
          'file@*_test.txt',
          shards.ShardSpec('file', 99_999, '_test.txt'),
      ),
      ('no_star', 'file-00000-of-00039', 'file@', None),
      ('first_file_missing', 'file-00003-00039', 'file@*', None),
      ('num_shards_too_high', 'file-00000-of-200039', 'file@*', None),
      ('too_many_digits_in_file_num', 'file-0000-of-39', 'file@*', None),
      ('not_enough_digits_in_file_num', 'file-000000-of-039', 'file@*', None),
      ('not_enough_digits_in_count', 'file-00000-of-039', 'file@*', None),
      ('too_many_digits_in_count', 'file-00000-of-000039', 'file@*', None),
  )
  def test_path_with_star(
      self,
      filename: str,
      shard_spec: str,
      expected_parsed_shard: shards.ShardSpec,
  ):
    tmp_dir = self.create_tempdir()
    pathlib.Path(tmp_dir).joinpath(filename).write_text('')

    if expected_parsed_shard:
      # Prefix shard spec with tmp_dir.
      expected_parsed_shard = shards.ShardSpec(
          os.path.join(tmp_dir, expected_parsed_shard.prefix),
          num_shards=expected_parsed_shard.num_shards,
          suffix=expected_parsed_shard.suffix,
      )

    self.assertEqual(
        shards.parse_shard_spec(os.path.join(tmp_dir, shard_spec)),
        expected_parsed_shard,
    )


if __name__ == '__main__':
  absltest.main()
