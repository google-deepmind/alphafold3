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

import functools
from absl.testing import absltest
from absl.testing import parameterized
from alphafold3.model.components import mapping
import haiku as hk
import jax
import jax.numpy as jnp
import numpy as np


class TestModule(hk.Module):

  def computation_fn(self, x, y, z):
    return jnp.einsum('qhc,khc->qkh', x, y) + z

  def __call__(self, x, y, z):
    return mapping.inference_subbatch(
        module=self.computation_fn,
        subbatch_size=2,
        batched_args=[x],
        nonbatched_args=[y, z],
    )


class TestInferenceSubbatch(parameterized.TestCase):

  def test_inference_subbatch_shape(self):
    """Tests shape of inference_subbatch output."""

    key1, key2, key3, key4 = jax.random.split(jax.random.PRNGKey(0), 4)
    x = jax.random.normal(key1, (8, 3, 2))
    y = jax.random.normal(key2, (7, 3, 2))
    z = jax.random.normal(key3, (1, 7, 3))

    @hk.transform
    def test_fn(x, y, z):
      return jnp.einsum('qhc,khc->qkh', x, y) + z

    apply_fn = functools.partial(test_fn.apply, {}, key4)
    output = apply_fn(x, y, z)

    self.assertEqual(output.shape, (8, 7, 3))

if __name__ == '__main__':
  absltest.main()
