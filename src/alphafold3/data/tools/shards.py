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

"""A library to handle shards of the format file_path@NUM_SHARDS.

For instance, /path/to/file@20 will generate the following shards:

- /path/to/file-00000-of-00020
- /path/to/file-00001-of-00020
- ...
- /path/to/file-00019-of-00020

This also supports @* pattern, which will determine the number of shards based
on the filesystem content.
"""

from collections.abc import Sequence
import dataclasses
import pathlib
import re


_MAX_NUM_SHARDS = 99_999
_SHARD_RE = re.compile(
    r"""
    ^(?P<prefix>[^\?\],\*]+)@
     (?P<shards>(\d{1,5})|\*)
     (?P<suffix>[\._][^\?\]@\*\/]*)?
    $""",
    re.X,
)


@dataclasses.dataclass(frozen=True)
class ShardSpec:
  prefix: str
  num_shards: int
  suffix: str


def parse_shard_spec(path: str) -> ShardSpec | None:
  """Returns the shard spec or None if the path is not a shard spec.

  For instance, if the shard spec is '/path/to/file@20', the output will be
  ('/path/to/file', 20).

  Args:
    path: the path to parse, e.g. /path/to/file@20 or /path/to/file@*.
  """
  parsed = re.fullmatch(_SHARD_RE, path)
  if not parsed:
    return None
  prefix = parsed.group('prefix')
  shards = parsed.group('shards')
  suffix = parsed.group('suffix') or ''

  if shards != '*':
    return ShardSpec(prefix=prefix, num_shards=int(shards), suffix=suffix)
  shard_slice = slice(len(prefix) + 10, len(prefix) + 15)
  shard_path = pathlib.Path(f'{prefix}-00000-of-?????{suffix}')
  for shard in sorted(shard_path.parent.glob(shard_path.name), reverse=True):
    try:
      num_shards = int(str(shard)[shard_slice])
      return ShardSpec(prefix=prefix, num_shards=num_shards, suffix=suffix)
    except ValueError:
      continue
  return None


def get_sharded_paths(shard_spec: str) -> Sequence[str] | None:
  """Returns a list of file path or None if the input is not a shard spec.

  Args:
    shard_spec: the specifications of the shard, e.g. /path/to/file@20.
  """
  parsed_spec = parse_shard_spec(shard_spec)
  if not parsed_spec:
    return None

  prefix = parsed_spec.prefix
  num_shards = parsed_spec.num_shards
  suffix = parsed_spec.suffix
  if num_shards > _MAX_NUM_SHARDS:
    raise ValueError(f'Shard count for {shard_spec} exceeds {_MAX_NUM_SHARDS}')
  return [
      f'{prefix}-{i:05d}-of-{num_shards:05d}{suffix}' for i in range(num_shards)
  ]
