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

"""Tests for the structure stores."""

from absl.testing import absltest
from alphafold3.data import structure_stores
import gzip
import io
import os
import tarfile
import tempfile


class StructureStoreTest(absltest.TestCase):

  def test_directory_store(self):
    with tempfile.TemporaryDirectory() as d:
      with open(os.path.join(d, '1abc.cif'), 'w') as f:
        f.write('data_1abc')
      store = structure_stores.StructureStore(d)
      self.assertEqual(store.get_mmcif_str('1abc'), 'data_1abc')
      self.assertIn('1abc', store.target_names())
      with self.assertRaises(structure_stores.NotFoundError):
        store.get_mmcif_str('missing')

  def test_tar_store_with_gzipped_mmcif(self):
    buf = io.BytesIO()
    gzdata = gzip.compress(b'data_gz')
    plain = b'data_plain'
    with tarfile.open(fileobj=buf, mode='w') as tar:
      for name, data in [('2xyz.cif.gz', gzdata), ('dir/3wxy.cif', plain)]:
        info = tarfile.TarInfo(name)
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    buf.seek(0)
    with tempfile.TemporaryDirectory() as d:
      tar_path = os.path.join(d, 'store.tar')
      with open(tar_path, 'wb') as f:
        f.write(buf.getvalue())
      store = structure_stores.StructureStore(tar_path)
      self.assertEqual(store.get_mmcif_str('2xyz'), 'data_gz')
      self.assertEqual(store.get_mmcif_str('3wxy'), 'data_plain')
      self.assertIn('2xyz', store.target_names())
      self.assertIn('3wxy', store.target_names())
      with self.assertRaises(structure_stores.NotFoundError):
        store.get_mmcif_str('missing')
      store.close()

  def test_directory_store_raises_not_found_not_io_error(self):
    with tempfile.TemporaryDirectory() as d:
      store = structure_stores.StructureStore(d)
      with self.assertRaises(structure_stores.NotFoundError):
        store.get_mmcif_str('does_not_exist')


if __name__ == '__main__':
  absltest.main()