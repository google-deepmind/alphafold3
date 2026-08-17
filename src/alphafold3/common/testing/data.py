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

"""Module that provides an abstraction for accessing test data."""

import os
import pathlib
from typing import Literal, overload

from absl.testing import absltest


class Data:
  """Provides an abstraction for accessing test data."""

  def __init__(self, data_dir: os.PathLike[str] | str):
    """Initiailizes data wrapper, providing users with high level data access.

    Args:
      data_dir: Directory containing test data.
    """
    self._data_dir = pathlib.Path(data_dir)

  def path(self, data_name: str | os.PathLike[str] | None = None) -> str:
    """Returns the path to a given test data.

    Args:
      data_name: the name of the test data file relative to data_dir. If not
        set, this will return the absolute path to the data directory.
    """
    data_dir_path = (
        pathlib.Path(absltest.get_default_test_srcdir()) / self._data_dir
    )

    if data_name:
      return str(data_dir_path / data_name)

    return str(data_dir_path)

  @overload
  def load(
      self, data_name: str | os.PathLike[str], mode: Literal['rt'] = 'rt'
  ) -> str:
    ...

  @overload
  def load(
      self, data_name: str | os.PathLike[str], mode: Literal['rb'] = 'rb'
  ) -> bytes:
    ...

  def load(
      self, data_name: str | os.PathLike[str], mode: str = 'rt'
  ) -> str | bytes:
    """Returns the contents of a given test data.

    Args:
      data_name: the name of the test data file relative to data_dir.
      mode: the mode in which to read the data file. Defaults to text ('rt').
    """
    with open(self.path(data_name), mode=mode) as f:
      return f.read()
