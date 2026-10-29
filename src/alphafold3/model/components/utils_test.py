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

from absl.testing import absltest
from alphafold3.model.components import utils
import jax
import numpy as np


class UtilsTest(absltest.TestCase):

  def test_remove_invalidly_typed_feats(self):
    batch = {
        'feat_float': np.ones(4, dtype=np.float32),
        'feat_str': np.array('Remove me!'),
        'feat_str2': np.array('Remove me too!', dtype=object),
        'feat_int': np.ones(5, dtype=np.int64),
        'feat_bool': np.array(True),
    }
    new_batch_dtypes = jax.tree.map(
        lambda x: x.dtype, utils.remove_invalidly_typed_feats(batch)  # pyrefly: ignore[bad-argument-type]
    )
    expected_dtypes = {
        'feat_float': np.float32,
        'feat_int': np.int64,
        'feat_bool': bool,
    }
    self.assertDictEqual(expected_dtypes, new_batch_dtypes)


if __name__ == '__main__':
  absltest.main()
