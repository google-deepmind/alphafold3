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
from absl.testing import parameterized
from alphafold3.jax.geometry import utils
import jax
import jax.numpy as jnp
import numpy as np


class UtilsTest(parameterized.TestCase):

  @parameterized.parameters((None, (16, 32)), (-1, (16, 32)), (-1, (32,)))
  def test_weighted_mean(self, axis, weights_shape):
    key = jax.random.PRNGKey(42)
    weights = jax.random.uniform(
        key, shape=weights_shape, dtype=jnp.float32, minval=0.5, maxval=1.0
    )
    values = jax.random.normal(key, shape=(16, 32), dtype=jnp.float32) + 3.0

    weighted_mean = utils.weighted_mean(
        weights=weights, value=values, eps=0.0, axis=axis
    )
    jnp_average = jnp.average(values, weights=weights, axis=axis)
    # Fairly loose tolerance due to jax not providing the best numerical
    # guarantees.
    np.testing.assert_almost_equal(
        np.asarray(weighted_mean, dtype=np.float64),
        np.asarray(jnp_average, dtype=np.float64),
        decimal=5,
    )


if __name__ == '__main__':
  absltest.main()
