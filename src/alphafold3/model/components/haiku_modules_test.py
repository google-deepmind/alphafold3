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
from alphafold3.model.components import haiku_modules
import chex
import haiku as hk
import jax
import jax.numpy as jnp

FLOAT_DTYPES = [
    jnp.float16,
    jnp.bfloat16,
    jnp.float32,
    jnp.float16,
    jnp.float64,
]


class ModulesTest(parameterized.TestCase):

  @parameterized.product(
      axis=[1, -1],
      upcast=[True, False],
      dtype=FLOAT_DTYPES,
  )
  def test_layernorm_upcast(self, axis, upcast, dtype):
    rng = jax.random.PRNGKey(0)
    x = jax.random.normal(rng, (2, 3, 4), dtype=dtype)

    def f(x):
      return haiku_modules.LayerNorm(
          axis=axis,
          create_scale=True,
          create_offset=True,
          upcast=upcast,
          use_fast_variance=False,
          name='layer_norm',
      )(x)

    def f_hk(x):
      dtype = x.dtype
      if upcast and dtype in [jnp.float16, jnp.bfloat16]:
        x = x.astype(jnp.float32)
      return hk.LayerNorm(
          axis=axis,
          create_scale=True,
          create_offset=True,
          use_fast_variance=False,
      )(x).astype(dtype)

    params = hk.transform(f).init(rng, x)

    out = jax.jit(hk.transform(f).apply)(params, rng, x)
    out_hk = jax.jit(hk.transform(f_hk).apply)(params, rng, x)

    chex.assert_trees_all_equal(out, out_hk)

  def test_layernorm_axis(self):
    rng = jax.random.PRNGKey(0)
    x = jax.random.normal(rng, (2, 3, 4), dtype=jnp.float32)

    def f(x):
      x = jnp.transpose(x, (2, 0, 1))
      x = haiku_modules.LayerNorm(
          axis=0,
          create_scale=True,
          create_offset=True,
          param_axis=0,
          name='layer_norm',
      )(x)
      return jnp.transpose(x, (1, 2, 0))

    def f_hk(x):
      return hk.LayerNorm(
          axis=-1, create_scale=True, create_offset=True, name='layer_norm'
      )(x)

    params = hk.transform(f).init(rng, x)
    params_hk = hk.transform(f_hk).init(rng, x)

    chex.assert_trees_all_close(params, params_hk)

    out = jax.jit(hk.transform(f).apply)(params, rng, x)
    out_hk = jax.jit(hk.transform(f_hk).apply)(params, rng, x)

    chex.assert_trees_all_close(out, out_hk)

  def test_linear_simple(self):
    bias = 1.4

    def f(x):
      return haiku_modules.Linear(
          bias_init=bias, use_bias=True, num_output=5, name='linear'
      )(x)

    f_trans = hk.transform(f)

    x = jnp.ones((3, 2))
    rng = jax.random.PRNGKey(0)
    params = f_trans.init(rng, x)
    out = f_trans.apply(params, rng, x)

    out_expected = jnp.dot(x, params['linear']['weights']) + bias
    chex.assert_trees_all_equal(out, out_expected)


if __name__ == '__main__':
  absltest.main()
