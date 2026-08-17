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

"""Library for loading structure data from various sources."""

from collections.abc import Mapping, Sequence
import functools
import tarfile
from etils import epath


class NotFoundError(KeyError):
  """Raised when the structure store doesn't contain the requested target."""


class StructureStore:
  """Handles the retrieval of mmCIF files from a filesystem."""

  def __init__(
      self,
      structures: epath.PathLike | Mapping[str, str],
  ):
    """Initialises the instance.

    Args:
      structures: Path of the directory where the mmCIF files are or a Mapping
        from target name to mmCIF string.
    """
    if isinstance(structures, Mapping):
      self._structure_mapping = structures
      self._structure_path = None
      self._structure_tar = None

    else:
      self._structure_mapping = None
      structures = epath.Path(structures)
      if structures.suffix == '.tar':
        self._structure_tar = tarfile.open(
            fileobj=structures.open('rb'),
            mode='r',
        )
        self._structure_path = None

      else:
        self._structure_path = structures
        self._structure_tar = None

  @functools.cached_property
  def _tar_members(self) -> Mapping[str, tarfile.TarInfo]:
    assert self._structure_tar is not None
    return {
        path.stem: tarinfo
        for tarinfo in self._structure_tar.getmembers()
        if tarinfo.isfile()
        and (path := epath.Path(tarinfo.path.lower())).suffix == '.cif'
    }

  def get_mmcif_str(self, target_name: str) -> str:
    """Returns an mmCIF for a given `target_name`.

    Args:
      target_name: Name specifying the target mmCIF.

    Raises:
      NotFoundError: If the target is not found.
    """
    if self._structure_mapping is not None:
      try:
        return self._structure_mapping[target_name]
      except KeyError as e:
        raise NotFoundError(f'{target_name=} not found') from e

    if self._structure_tar is not None:
      try:
        member = self._tar_members[target_name]
        if struct_file := self._structure_tar.extractfile(member):
          return struct_file.read().decode()
        else:
          raise NotFoundError(f'{target_name=} not found')
      except KeyError:
        raise NotFoundError(f'{target_name=} not found') from None

    filepath = self._structure_path / f'{target_name}.cif'  # pyrefly: ignore[unsupported-operation]
    try:
      return filepath.read_text()
    except Exception as e:
      # Unfortunately, we can't predict which error type will be raised from the
      # underlying storage library (e.g. tensorflow or gcsfs), so this is an
      # attempt to stay backward compatible with file not found without
      # obscuring other possible error conditions
      exc_str = str(e).lower()
      if 'no such file' in exc_str or 'not found' in exc_str:
        raise NotFoundError(f'{target_name=} not found at {filepath=}') from e

      raise IOError(f'Error reading file {filepath}: {e}') from e

  def target_names(self) -> Sequence[str]:
    """Returns all targets in the store."""
    if self._structure_mapping is not None:
      return [*self._structure_mapping.keys()]
    elif self._structure_tar is not None:
      return sorted(self._tar_members.keys())
    elif self._structure_path is not None:
      return sorted([path.stem for path in self._structure_path.glob('*.cif')])
    return ()
