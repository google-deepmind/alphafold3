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
from alphafold3.model import features
from alphafold3.model import model_config
from alphafold3.model.network import modules
from alphafold3.model.network import template_modules
import haiku as hk
import jax
import jax.numpy as jnp
import numpy as np

NUM_RES = 5
NUM_PAIR_CHANNEL = 3
NUM_TEMPLATES = 2
NUM_DENSE = 24


class TemplateModulesTest(absltest.TestCase):

  @hk.testing.transform_and_run
  def test_finite(self):
    config = template_modules.TemplateEmbedding.Config(
        template_stack=modules.PairFormerIteration.Config(
            num_layer=1,
            pair_attention=modules.GridSelfAttention.Config(
                num_head=2,
            ),
            pair_transition=modules.TransitionBlock.Config(
                num_intermediate_factor=2, use_glu_kernel=False
            ),
            triangle_multiplication_incoming=modules.TriangleMultiplication.Config(
                equation='kjc,kic->ijc', use_glu_kernel=False
            ),
            triangle_multiplication_outgoing=modules.TriangleMultiplication.Config(
                equation='ikc,jkc->ijc', use_glu_kernel=False
            ),
        ),
        dgram_features=template_modules.DistogramFeaturesConfig(
            max_bin=50.75,
            min_bin=3.25,
            num_bins=39,
        ),
        num_channels=8,
    )

    global_config = model_config.GlobalConfig(
        bfloat16='all',
        final_init='linear',
        pair_attention_chunk_size=[(None, 4)],
        pair_transition_shard_spec=[(None, None)],
        flash_attention_implementation='xla',
    )
    template_module = template_modules.TemplateEmbedding(config, global_config)

    keys = jax.random.split(jax.random.PRNGKey(42), 5)
    query_embedding = jax.random.normal(
        keys[0], [NUM_RES, NUM_RES, NUM_PAIR_CHANNEL]
    )

    template_batch = {
        'template_aatype': jax.random.randint(
            keys[2], [NUM_TEMPLATES, NUM_RES], 0, 21
        ),
        'template_atom_positions': jax.random.normal(
            keys[1], [NUM_TEMPLATES, NUM_RES, NUM_DENSE, 3]
        ),
        'template_atom_mask': jnp.ones([NUM_TEMPLATES, NUM_RES, NUM_DENSE]),
    }
    template_batch = features.Templates.from_data_dict(template_batch)

    padding_mask_2d = jnp.ones(shape=[NUM_RES, NUM_RES], dtype=jnp.float32)
    multichain_mask_2d = jnp.ones([NUM_RES, NUM_RES])

    result = template_module(
        query_embedding,
        template_batch,
        padding_mask_2d,
        multichain_mask_2d,
        key=keys[4],
    )
    self.assertTrue(np.isfinite(result).all())


if __name__ == '__main__':
  absltest.main()
