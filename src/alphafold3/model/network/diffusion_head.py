# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Diffusion Head."""

from collections.abc import Callable

from alphafold3.common import base_config
from alphafold3.model import feat_batch
from alphafold3.model import model_config
from alphafold3.model.components import haiku_modules as hm
from alphafold3.model.components import utils
from alphafold3.model.network import atom_cross_attention
from alphafold3.model.network import diffusion_transformer
from alphafold3.model.network import featurization
from alphafold3.model.network import noise_level_embeddings
import chex
import haiku as hk
import jax
import jax.numpy as jnp


# Carefully measured by averaging multimer training set.
SIGMA_DATA = 16.0


def random_rotation(key):
  # Create a random rotation (Gram-Schmidt orthogonalization of two
  # random normal vectors)
  v0, v1 = jax.random.normal(key, shape=(2, 3))
  e0 = v0 / jnp.maximum(1e-10, jnp.linalg.norm(v0))
  v1 = v1 - e0 * jnp.dot(v1, e0, precision=jax.lax.Precision.HIGHEST)
  e1 = v1 / jnp.maximum(1e-10, jnp.linalg.norm(v1))
  e2 = jnp.cross(e0, e1)
  return jnp.stack([e0, e1, e2])


def random_augmentation(
    rng_key: jnp.ndarray,
    positions: jnp.ndarray,
    mask: jnp.ndarray,
) -> jnp.ndarray:
  """Apply random rigid augmentation.

  Args:
    rng_key: random key
    positions: atom positions of shape (<common_axes>, 3)
    mask: per-atom mask of shape (<common_axes>,)

  Returns:
    Transformed positions with the same shape as input positions.
  """
  rotation_key, translation_key = jax.random.split(rng_key)

  center = utils.mask_mean(
      mask[..., None], positions, axis=(-2, -3), keepdims=True, eps=1e-6
  )
  rot = random_rotation(rotation_key)
  translation = jax.random.normal(translation_key, shape=(3,))

  augmented_positions = (
      jnp.einsum(
          '...i,ij->...j',
          positions - center,
          rot,
          precision=jax.lax.Precision.HIGHEST,
      )
      + translation
  )
  return augmented_positions * mask[..., None]


def noise_schedule(t, smin=0.0004, smax=160.0, p=7):
  return (
      SIGMA_DATA
      * (smax ** (1 / p) + t * (smin ** (1 / p) - smax ** (1 / p))) ** p
  )


class ConditioningConfig(base_config.BaseConfig):
  pair_channel: int
  seq_channel: int
  prob: float


class SampleConfig(base_config.BaseConfig):
  steps: int
  gamma_0: float = 0.8
  gamma_min: float = 1.0
  noise_scale: float = 1.003
  step_scale: float = 1.5
  num_samples: int = 1


class DiffusionHead(hk.Module):
  """Denoising Diffusion Head."""

  class Config(
      atom_cross_attention.AtomCrossAttEncoderConfig,
      atom_cross_attention.AtomCrossAttDecoderConfig,
  ):
    """Configuration for DiffusionHead."""

    eval_batch_size: int = 5
    eval_batch_dim_shard_size: int = 5
    conditioning: ConditioningConfig = base_config.autocreate(
        prob=0.8, pair_channel=128, seq_channel=384
    )
    eval: SampleConfig = base_config.autocreate(
        num_samples=5,
        steps=200,
    )
    transformer: diffusion_transformer.Transformer.Config = (
        base_config.autocreate()
    )

  def __init__(
      self,
      config: Config,
      global_config: model_config.GlobalConfig,
      name='diffusion_head',
  ):
    self.config = config
    self.global_config = global_config
    super().__init__(name=name)

  @hk.transparent
  def _conditioning(
      self,
      batch: feat_batch.Batch,
      embeddings: dict[str, jnp.ndarray],
      noise_level: jnp.ndarray,
      use_conditioning: bool,
  ) -> tuple[jnp.ndarray, jnp.ndarray]:
    single_embedding = use_conditioning * embeddings['single']
    pair_embedding = use_conditioning * embeddings['pair']

    rel_features = featurization.create_relative_encoding(
        batch.token_features, max_relative_idx=32, max_relative_chain=2
    ).astype(pair_embedding.dtype)
    features_2d = jnp.concatenate([pair_embedding, rel_features], axis=-1)
    pair_cond = hm.Linear(
        self.config.conditioning.pair_channel,
        precision='highest',
        name='pair_cond_initial_projection',
    )(
        hm.LayerNorm(
            use_fast_variance=False,
            create_offset=False,
            name='pair_cond_initial_norm',
        )(features_2d)
    )

    for idx in range(2):
      pair_cond += diffusion_transformer.transition_block(
          pair_cond, 2, self.global_config, name=f'pair_transition_{idx}'
      )

    target_feat = embeddings['target_feat']
    features_1d = jnp.concatenate([single_embedding, target_feat], axis=-1)
    single_cond = hm.LayerNorm(
        use_fast_variance=False,
        create_offset=False,
        name='single_cond_initial_norm',
    )(features_1d)
    single_cond = hm.Linear(
        self.config.conditioning.seq_channel,
        precision='highest',
        name='single_cond_initial_projection',
    )(single_cond)

    noise_embedding = noise_level_embeddings.noise_embeddings(
        sigma_scaled_noise_level=noise_level / SIGMA_DATA
    )
    single_cond += hm.Linear(
        self.config.conditioning.seq_channel,
        precision='highest',
        name='noise_embedding_initial_projection',
    )(
        hm.LayerNorm(
            use_fast_variance=False,
            create_offset=False,
            name='noise_embedding_initial_norm',
        )(noise_embedding)
    )

    for idx in range(2):
      single_cond += diffusion_transformer.transition_block(
          single_cond, 2, self.global_config, name=f'single_transition_{idx}'
      )

    return single_cond, pair_cond

  def __call__(
      self,
      # positions_noisy.shape: (num_token, max_atoms_per_token, 3)
      positions_noisy: jnp.ndarray,
      noise_level: jnp.ndarray,
      batch: feat_batch.Batch,
      embeddings: dict[str, jnp.ndarray],
      use_conditioning: bool,
  ) -> jnp.ndarray:

    with utils.bfloat16_context():
      # Get conditioning
      trunk_single_cond, trunk_pair_cond = self._conditioning(
          batch=batch,
          embeddings=embeddings,
          noise_level=noise_level,
          use_conditioning=use_conditioning,
      )

      # Extract features
      sequence_mask = batch.token_features.mask
      atom_mask = batch.predicted_structure_info.atom_mask

      # Position features
      act = positions_noisy * atom_mask[..., None]
      act = act / jnp.sqrt(noise_level**2 + SIGMA_DATA**2)

      enc = atom_cross_attention.atom_cross_att_encoder(
          token_atoms_act=act,
          trunk_single_cond=embeddings['single'],
          trunk_pair_cond=trunk_pair_cond,
          config=self.config,
          global_config=self.global_config,
          batch=batch,
          name='diffusion',
      )
      act = enc.token_act

      # Token-token attention
      chex.assert_shape(act, (None, self.config.per_token_channels))
      act = jnp.asarray(act, dtype=jnp.float32)

      act += hm.Linear(
          act.shape[-1],
          precision='highest',
          initializer=self.global_config.final_init,
          name='single_cond_embedding_projection',
      )(
          hm.LayerNorm(
              use_fast_variance=False,
              create_offset=False,
              name='single_cond_embedding_norm',
          )(trunk_single_cond)
      )

      act = jnp.asarray(act, dtype=jnp.float32)
      trunk_single_cond = jnp.asarray(trunk_single_cond, dtype=jnp.float32)
      trunk_pair_cond = jnp.asarray(trunk_pair_cond, dtype=jnp.float32)
      sequence_mask = jnp.asarray(sequence_mask, dtype=jnp.float32)

      transformer = diffusion_transformer.Transformer(
          self.config.transformer, self.global_config
      )
      act = transformer(
          act=act,
          single_cond=trunk_single_cond,
          mask=sequence_mask,
          pair_cond=trunk_pair_cond,
      )
      act = hm.LayerNorm(
          use_fast_variance=False, create_offset=False, name='output_norm'
      )(act)
      # (n_tokens, per_token_channels)

      # (Possibly) atom-granularity decoder
      assert isinstance(enc, atom_cross_attention.AtomCrossAttEncoderOutput)
      position_update = atom_cross_attention.atom_cross_att_decoder(
          token_act=act,
          enc=enc,
          config=self.config,
          global_config=self.global_config,
          batch=batch,
          name='diffusion',
      )

      skip_scaling = SIGMA_DATA**2 / (noise_level**2 + SIGMA_DATA**2)
      out_scaling = (
          noise_level * SIGMA_DATA / jnp.sqrt(noise_level**2 + SIGMA_DATA**2)
      )
    # End `with utils.bfloat16_context()`.

    return (
        skip_scaling * positions_noisy + out_scaling * position_update
    ) * atom_mask[..., None]


def sample(
    denoising_step: Callable[[jnp.ndarray, jnp.ndarray], jnp.ndarray],
    batch: feat_batch.Batch,
    key: jnp.ndarray,
    config: SampleConfig,
) -> dict[str, jnp.ndarray]:
  """Sample using denoiser on batch.

  Args:
    denoising_step: the denoising function.
    batch: the batch
    key: random key
    config: config for the sampling process (e.g. number of denoising steps,
      etc.)

  Returns:
    a dict
      {
         'atom_positions': jnp.array(...)       # shape (<common_axes>, 3)
         'mask': jnp.array(...)                 # shape (<common_axes>,)
      }
    where the <common_axes> are
    (num_samples, num_tokens, max_atoms_per_token)
  """

  mask = batch.predicted_structure_info.atom_mask

  def apply_denoising_step(carry, noise_level):
    key, positions, noise_level_prev = carry
    key, key_noise, key_aug = jax.random.split(key, 3)

    positions = random_augmentation(
        rng_key=key_aug, positions=positions, mask=mask
    )

    gamma = config.gamma_0 * (noise_level > config.gamma_min)
    t_hat = noise_level_prev * (1 + gamma)

    noise_scale = config.noise_scale * jnp.sqrt(t_hat**2 - noise_level_prev**2)
    noise = noise_scale * jax.random.normal(key_noise, positions.shape)
    positions_noisy = positions + noise

    positions_denoised = denoising_step(positions_noisy, t_hat)
    grad = (positions_noisy - positions_denoised) / t_hat

    d_t = noise_level - t_hat
    positions_out = positions_noisy + config.step_scale * d_t * grad

    return (key, positions_out, noise_level), positions_out

  num_samples = config.num_samples

  noise_levels = noise_schedule(jnp.linspace(0, 1, config.steps + 1))

  key, noise_key = jax.random.split(key)
  positions = jax.random.normal(noise_key, (num_samples,) + mask.shape + (3,))
  positions *= noise_levels[0]

  init = (
      jax.random.split(key, num_samples),
      positions,
      jnp.tile(noise_levels[None, 0], (num_samples,)),
  )

  apply_denoising_step = hk.vmap(
      apply_denoising_step, in_axes=(0, None), split_rng=(not hk.running_init())
  )
  result, _ = hk.scan(apply_denoising_step, init, noise_levels[1:], unroll=4)
  _, positions_out, _ = result

  final_dense_atom_mask = jnp.tile(mask[None], (num_samples, 1, 1))

  return {'atom_positions': positions_out, 'mask': final_dense_atom_mask}
