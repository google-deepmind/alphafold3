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

from collections.abc import Iterator, Mapping
import contextlib
import io
import os
import tarfile
import tempfile
from absl.testing import absltest
from alphafold3.data import structure_stores


@contextlib.contextmanager
def _create_store(num_targets: int) -> Iterator[str]:
  tempdir = absltest.get_default_test_tmpdir()
  with tempfile.TemporaryDirectory(dir=tempdir) as root:
    for i in range(num_targets):
      with open(os.path.join(root, f'{i:04}.cif'), 'w') as f:
        f.write(f'data_{i:04}')
    yield root


@contextlib.contextmanager
def _create_tar_store(num_targets: int) -> Iterator[str]:
  tempdir = absltest.get_default_test_tmpdir()
  with tempfile.TemporaryDirectory(dir=tempdir) as root:
    tar_path = os.path.join(root, 'test.tar')
    with tarfile.open(tar_path, 'w') as tar:
      tarinfo = tarfile.TarInfo('mmcif_files')
      tarinfo.size = 0
      tarinfo.type = tarfile.DIRTYPE
      tar.addfile(tarinfo)
      for i in range(num_targets):
        tarinfo = tarfile.TarInfo(f'mmcif_files/{i:04}.cif')
        data = f'data_{i:04}'.encode()
        tarinfo.size = len(data)
        tar.addfile(tarinfo, io.BytesIO(data))

    yield tar_path


def _create_in_memory_store(num_targets: int) -> Mapping[str, str]:
  return {f'{i:04}': f'data_{i:04}' for i in range(num_targets)}


class StructureStoreFilesystemTest(absltest.TestCase):

  def test_get_mmcif_not_found(self):
    with _create_store(5) as store_path:
      store = structure_stores.StructureStore(store_path)
      with self.assertRaisesRegex(
          structure_stores.NotFoundError,
          "target_name='00010' not found at filepath=.*00010.cif",
      ):
        store.get_mmcif_str('00010')

  def test_get_mmcifs(self):
    with _create_store(5) as store_path:
      store = structure_stores.StructureStore(store_path)
      self.assertEqual(
          [store.get_mmcif_str(f'{i:04}') for i in range(5)],
          [f'data_{i:04}' for i in range(5)],
      )

  def test_target_names(self):
    with _create_store(5) as store_path:
      store = structure_stores.StructureStore(store_path)
      self.assertEqual(store.target_names(), [f'{i:04}' for i in range(5)])

  def test_tar_get_mmcif_not_found(self):
    with _create_tar_store(5) as store_path:
      store = structure_stores.StructureStore(store_path)
      with self.assertRaisesRegex(
          structure_stores.NotFoundError,
          "target_name='00010' not found",
      ):
        store.get_mmcif_str('00010')

  def test_tar_get_mmcifs(self):
    with _create_tar_store(5) as store_path:
      store = structure_stores.StructureStore(store_path)
      self.assertEqual(
          [store.get_mmcif_str(f'{i:04}') for i in range(5)],
          [f'data_{i:04}' for i in range(5)],
      )

  def test_tar_target_names(self):
    with _create_tar_store(5) as store_path:
      store = structure_stores.StructureStore(store_path)
      self.assertEqual(store.target_names(), [f'{i:04}' for i in range(5)])

  def test_in_memory_get_mmcif_not_found(self):
    store = structure_stores.StructureStore(_create_in_memory_store(5))
    with self.assertRaisesRegex(
        structure_stores.NotFoundError,
        "target_name='00010' not found",
    ):
      store.get_mmcif_str('00010')

  def test_in_memory_get_mmcifs(self):
    store = structure_stores.StructureStore(_create_in_memory_store(5))
    self.assertEqual(
        [store.get_mmcif_str(f'{i:04}') for i in range(5)],
        [f'data_{i:04}' for i in range(5)],
    )

  def test_in_memory_target_names(self):
    store = structure_stores.StructureStore(_create_in_memory_store(5))
    self.assertSequenceEqual(
        store.target_names(), [f'{i:04}' for i in range(5)]
    )


if __name__ == '__main__':
  absltest.main()
