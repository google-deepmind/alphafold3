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
from alphafold3.model import model_config
from alphafold3.model.network import modules
import haiku as hk
import jax
import jax.numpy as jnp
import numpy as np


def get_padded_array(rng_key, shape, padding_mask, dtype, zeros, pad_fill=1e7):
  if zeros:
    array = jnp.zeros(shape, dtype)
  else:
    array = jax.random.uniform(rng_key, shape, dtype=dtype)

  pad_values = jax.lax.full_like(array, pad_fill)

  return jnp.where(padding_mask, pad_values, array)


def make_pad_mask(shape, padding):

  assert len(padding) == len(shape)
  shape_without_padding = tuple(
      (shape_dim - pad_dim) for shape_dim, pad_dim in zip(shape, padding)
  )
  mask = jnp.zeros(shape_without_padding, dtype=jnp.bool_)
  pad_config = tuple((0, pad_dim, 0) for pad_dim in padding)
  return jax.lax.pad(mask, True, pad_config)


class OuterProductMeanTest(parameterized.TestCase):

  def test_zero_init(self):
    num_output_channel = 5
    num_msa = 12
    num_res = 13
    num_channel = 4

    act_shape = (num_msa, num_res, num_channel)

    def apply_opm(act, mask):
      module = modules.OuterProductMean(
          modules.OuterProductMean.Config(chunk_size=4, num_outer_channel=4),
          global_config=model_config.GlobalConfig(final_init='zeros'),
          num_output_channel=num_output_channel,
          name='outer_product_mean',
      )
      return module(act, mask)

    act = jax.random.normal(
        jax.random.PRNGKey(42), act_shape, dtype=jnp.float32
    )

    mask = jnp.ones(act_shape[:-1], dtype=jnp.float32)

    init_fn, apply_fn = hk.transform(apply_opm)

    params = init_fn(jax.random.PRNGKey(428), act, mask)

    out = np.asarray(apply_fn(params, jax.random.PRNGKey(42), act, mask))

    np.testing.assert_equal(out, 0.0, 'Zero init should produce zero outputs')

  def test_padding(self):
    num_output_channel = 5
    num_msa = 12
    num_res = 13
    num_channel = 4

    pad_msa = 4
    pad_res = 5

    act_shape = (num_msa, num_res, num_channel)

    def apply_opm(act, mask):
      module = modules.OuterProductMean(
          modules.OuterProductMean.Config(chunk_size=4, num_outer_channel=4),
          global_config=model_config.GlobalConfig(final_init='linear'),
          num_output_channel=num_output_channel,
          name='outer_product_mean',
      )
      return module(act, mask)

    pad_bool_mask = make_pad_mask(act_shape[:-1], (pad_msa, pad_res))

    mask = pad_bool_mask.astype(jnp.float32)

    padded_array = get_padded_array(
        jax.random.PRNGKey(42),
        act_shape,
        pad_bool_mask[:, :, None],
        dtype=jnp.float32,
        zeros=False,
    )

    init_fn, apply_fn = hk.transform(apply_opm)

    params = init_fn(jax.random.PRNGKey(428), padded_array, mask)

    padded_out = np.asarray(
        apply_fn(params, jax.random.PRNGKey(42), padded_array, mask)
    )[:-pad_res, :-pad_res]

    unpadded_array = padded_array[:-pad_msa, :-pad_res]
    unpadded_mask = mask[:-pad_msa, :-pad_res]
    unpadded_out = np.asarray(
        apply_fn(params, jax.random.PRNGKey(42), unpadded_array, unpadded_mask)
    )

    np.testing.assert_array_almost_equal(padded_out, unpadded_out)


if __name__ == '__main__':
  absltest.main()
