# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""AlphaFold3 model."""

from collections.abc import Iterable, Mapping
import concurrent
import dataclasses
import functools
from typing import Any, TypeAlias

from absl import logging
from alphafold3 import structure
from alphafold3.common import base_config
from alphafold3.model import confidences
from alphafold3.model import feat_batch
from alphafold3.model import features
from alphafold3.model import model_config
from alphafold3.model.atom_layout import atom_layout
from alphafold3.model.components import mapping
from alphafold3.model.components import utils
from alphafold3.model.network import atom_cross_attention
from alphafold3.model.network import confidence_head
from alphafold3.model.network import diffusion_head
from alphafold3.model.network import distogram_head
from alphafold3.model.network import evoformer as evoformer_network
from alphafold3.model.network import featurization
import haiku as hk
import jax
import jax.numpy as jnp
import numpy as np


ModelResult: TypeAlias = Mapping[str, Any]


@dataclasses.dataclass(frozen=True, kw_only=True)
class InferenceResult:
  """Postprocessed model result.

  Attributes:
    predicted_structure: Predicted protein structure.
    numerical_data: Useful numerical data (scalars or arrays) to be saved at
      inference time.
    metadata: Smaller numerical data (usually scalar) to be saved as inference
      metadata.
    debug_outputs: Additional dict for debugging, e.g. raw outputs of a model
      forward pass.
    model_id: Model identifier.
  """

  predicted_structure: structure.Structure = dataclasses.field()
  numerical_data: Mapping[str, float | int | np.ndarray] = dataclasses.field(
      default_factory=dict
  )
  metadata: Mapping[str, float | int | np.ndarray] = dataclasses.field(
      default_factory=dict
  )
  debug_outputs: Mapping[str, Any] = dataclasses.field(default_factory=dict)
  model_id: bytes = b''


def get_predicted_structure(
    result: ModelResult, batch: feat_batch.Batch
) -> structure.Structure:
  """Creates the predicted structure and ion preditions.

  Args:
    result: model output in a model specific layout
    batch: model input batch

  Returns:
    Predicted structure.
  """
  model_output_coords = result['diffusion_samples']['atom_positions']

  # Rearrange model output coordinates to the flat output layout.
  model_output_to_flat = atom_layout.compute_gather_idxs(
      source_layout=batch.convert_model_output.token_atoms_layout,
      target_layout=batch.convert_model_output.flat_output_layout,
  )
  pred_flat_atom_coords = atom_layout.convert(
      gather_info=model_output_to_flat,
      arr=model_output_coords,
      layout_axes=(-3, -2),
  )

  predicted_lddt = result.get('predicted_lddt')

  if predicted_lddt is not None:
    pred_flat_b_factors = atom_layout.convert(
        gather_info=model_output_to_flat,
        arr=predicted_lddt,
        layout_axes=(-2, -1),
    )
  else:
    # Handle models which don't have predicted_lddt outputs.
    pred_flat_b_factors = np.zeros(pred_flat_atom_coords.shape[:-1])

  (missing_atoms_indices,) = np.nonzero(model_output_to_flat.gather_mask == 0)
  if missing_atoms_indices.shape[0] > 0:
    missing_atoms_flat_layout = batch.convert_model_output.flat_output_layout[
        missing_atoms_indices
    ]
    missing_atoms_uids = list(
        zip(
            missing_atoms_flat_layout.chain_id,
            missing_atoms_flat_layout.res_id,
            missing_atoms_flat_layout.res_name,
            missing_atoms_flat_layout.atom_name,
        )
    )
    logging.warning(
        'Target %s: warning: %s atoms were not predicted by the '
        'model, setting their coordinates to (0, 0, 0). '
        'Missing atoms: %s',
        batch.convert_model_output.empty_output_struc.name,
        missing_atoms_indices.shape[0],
        missing_atoms_uids,
    )

  # Put them into a structure
  pred_struc = batch.convert_model_output.empty_output_struc
  pred_struc = pred_struc.copy_and_update_atoms(
      atom_x=pred_flat_atom_coords[..., 0],
      atom_y=pred_flat_atom_coords[..., 1],
      atom_z=pred_flat_atom_coords[..., 2],
      atom_b_factor=pred_flat_b_factors,
      atom_occupancy=np.ones(pred_flat_atom_coords.shape[:-1]),  # Always 1.0.
  )
  # Set manually/differently when adding metadata.
  pred_struc = pred_struc.copy_and_update_globals(release_date=None)
  return pred_struc


def create_target_feat_embedding(
    batch: feat_batch.Batch,
    config: evoformer_network.Evoformer.Config,
    global_config: model_config.GlobalConfig,
) -> jnp.ndarray:
  """Create target feature embedding."""

  dtype = jnp.bfloat16 if global_config.bfloat16 == 'all' else jnp.float32

  with utils.bfloat16_context():
    target_feat = featurization.create_target_feat(
        batch,
        append_per_atom_features=False,
    ).astype(dtype)

    enc = atom_cross_attention.atom_cross_att_encoder(
        token_atoms_act=None,
        trunk_single_cond=None,
        trunk_pair_cond=None,
        config=config.per_atom_conditioning,
        global_config=global_config,
        batch=batch,
        name='evoformer_conditioning',
    )
    target_feat = jnp.concatenate([target_feat, enc.token_act], axis=-1).astype(
        dtype
    )

  return target_feat


def _compute_ptm(
    result: ModelResult,
    num_tokens: int,
    asym_id: np.ndarray,
    pae_single_mask: np.ndarray,
    interface: bool,
) -> np.ndarray:
  """Computes the pTM metrics from PAE."""
  return np.stack(
      [
          confidences.predicted_tm_score(
              tm_adjusted_pae=tm_adjusted_pae[:num_tokens, :num_tokens],
              asym_id=asym_id,
              pair_mask=pae_single_mask[:num_tokens, :num_tokens],
              interface=interface,
          )
          for tm_adjusted_pae in result['tmscore_adjusted_pae_global']
      ],
      axis=0,
  )


def _compute_chain_pair_iptm(
    num_tokens: int,
    asym_ids: np.ndarray,
    mask: np.ndarray,
    tm_adjusted_pae: np.ndarray,
) -> np.ndarray:
  """Computes the chain pair ipTM metrics from PAE."""
  return np.stack(
      [
          confidences.chain_pairwise_predicted_tm_scores(
              tm_adjusted_pae=sample_tm_adjusted_pae[:num_tokens],
              asym_id=asym_ids[:num_tokens],
              pair_mask=mask[:num_tokens, :num_tokens],
          )
          for sample_tm_adjusted_pae in tm_adjusted_pae
      ],
      axis=0,
  )


class Model(hk.Module):
  """Full model. Takes in data batch and returns model outputs."""

  class HeadsConfig(base_config.BaseConfig):
    diffusion: diffusion_head.DiffusionHead.Config = base_config.autocreate()
    confidence: confidence_head.ConfidenceHead.Config = base_config.autocreate()
    distogram: distogram_head.DistogramHead.Config = base_config.autocreate()

  class Config(base_config.BaseConfig):
    evoformer: evoformer_network.Evoformer.Config = base_config.autocreate()
    global_config: model_config.GlobalConfig = base_config.autocreate()
    heads: 'Model.HeadsConfig' = base_config.autocreate()
    num_recycles: int = 10
    return_embeddings: bool = False
    return_distogram: bool = False

  def __init__(self, config: Config, name: str = 'diffuser'):
    super().__init__(name=name)
    self.config = config
    self.global_config = config.global_config
    self.diffusion_module = diffusion_head.DiffusionHead(
        self.config.heads.diffusion, self.global_config
    )

  @hk.transparent
  def _sample_diffusion(
      self,
      batch: feat_batch.Batch,
      embeddings: dict[str, jnp.ndarray],
      *,
      sample_config: diffusion_head.SampleConfig,
  ) -> dict[str, jnp.ndarray]:
    denoising_step = functools.partial(
        self.diffusion_module,
        batch=batch,
        embeddings=embeddings,
        use_conditioning=True,
    )

    sample = diffusion_head.sample(
        denoising_step=denoising_step,
        batch=batch,
        key=hk.next_rng_key(),
        config=sample_config,
    )
    return sample

  def __call__(
      self, batch: features.BatchDict, key: jax.Array | None = None
  ) -> ModelResult:
    if key is None:
      key = hk.next_rng_key()

    batch = feat_batch.Batch.from_data_dict(batch)

    embedding_module = evoformer_network.Evoformer(
        self.config.evoformer, self.global_config
    )
    target_feat = create_target_feat_embedding(
        batch=batch,
        config=embedding_module.config,
        global_config=self.global_config,
    )

    def recycle_body(_, args):
      prev, key = args
      key, subkey = jax.random.split(key)
      embeddings = embedding_module(
          batch=batch,
          prev=prev,
          target_feat=target_feat,
          key=subkey,
      )
      embeddings['pair'] = embeddings['pair'].astype(jnp.float32)
      embeddings['single'] = embeddings['single'].astype(jnp.float32)
      return embeddings, key

    num_res = batch.num_res

    embeddings = {
        'pair': jnp.zeros(
            [num_res, num_res, self.config.evoformer.pair_channel],
            dtype=jnp.float32,
        ),
        'single': jnp.zeros(
            [num_res, self.config.evoformer.seq_channel], dtype=jnp.float32
        ),
        'target_feat': target_feat,
    }
    if hk.running_init():
      embeddings, _ = recycle_body(None, (embeddings, key))
    else:
      # Number of recycles is number of additional forward trunk passes.
      num_iter = self.config.num_recycles + 1
      embeddings, _ = hk.fori_loop(0, num_iter, recycle_body, (embeddings, key))

    samples = self._sample_diffusion(
        batch,
        embeddings,
        sample_config=self.config.heads.diffusion.eval,
    )

    # Compute dist_error_fn over all samples for distance error logging.
    confidence_output = mapping.sharded_map(
        lambda dense_atom_positions: confidence_head.ConfidenceHead(
            self.config.heads.confidence, self.global_config
        )(
            dense_atom_positions=dense_atom_positions,
            embeddings=embeddings,
            seq_mask=batch.token_features.mask,
            token_atoms_to_pseudo_beta=batch.pseudo_beta_info.token_atoms_to_pseudo_beta,
            asym_id=batch.token_features.asym_id,
        ),
        in_axes=0,
    )(samples['atom_positions'])

    distogram = distogram_head.DistogramHead(
        self.config.heads.distogram, self.global_config
    )(batch, embeddings, return_distogram=self.config.return_distogram)

    output = {
        'diffusion_samples': samples,
        'distogram': distogram,
        **confidence_output,
    }
    if self.config.return_embeddings:
      output['single_embeddings'] = embeddings['single']
      output['pair_embeddings'] = embeddings['pair']
    return output

  @classmethod
  def get_inference_result(
      cls,
      batch: features.BatchDict,
      result: ModelResult,
      target_name: str = '',
  ) -> Iterable[InferenceResult]:
    """Get the predicted structure, scalars, and arrays for inference.

    This function also computes any inference-time quantities, which are not a
    part of the forward-pass, e.g. additional confidence scores. Note that this
    function is not serialized, so it should be slim if possible.

    Args:
      batch: data batch used for model inference, incl. TPU invalid types.
      result: output dict from the model's forward pass.
      target_name: target name to be saved within structure.

    Yields:
      inference_result: dataclass object that contains a predicted structure,
      important inference-time scalars and arrays, as well as a slightly trimmed
      dictionary of raw model result from the forward pass (for debugging).
    """
    del target_name
    batch = feat_batch.Batch.from_data_dict(batch)

    # Retrieve structure and construct a predicted structure.
    pred_structure = get_predicted_structure(result=result, batch=batch)

    num_tokens = batch.token_features.seq_length.item()

    pae_single_mask = np.tile(
        batch.frames.mask[:, None],
        [1, batch.frames.mask.shape[0]],
    )
    ptm = _compute_ptm(
        result=result,
        num_tokens=num_tokens,
        asym_id=batch.token_features.asym_id[:num_tokens],
        pae_single_mask=pae_single_mask,
        interface=False,
    )
    iptm = _compute_ptm(
        result=result,
        num_tokens=num_tokens,
        asym_id=batch.token_features.asym_id[:num_tokens],
        pae_single_mask=pae_single_mask,
        interface=True,
    )
    ptm_iptm_average = 0.8 * iptm + 0.2 * ptm

    asym_ids = batch.token_features.asym_id[:num_tokens]
    # Map asym IDs back to chain IDs. Asym IDs are constructed from chain IDs by
    # iterating over the chain IDs, and for each unique chain ID incrementing
    # the asym ID by 1 and mapping it to the particular chain ID. Asym IDs are
    # 1-indexed, so subtract 1 to get back to the chain ID.
    chain_ids = [pred_structure.chains[asym_id - 1] for asym_id in asym_ids]
    res_ids = batch.token_features.residue_index[:num_tokens]

    if len(np.unique(asym_ids[:num_tokens])) > 1:
      # There is more than one chain, hence interface pTM (i.e. ipTM) defined,
      # so use it.
      ranking_confidence = ptm_iptm_average
    else:
      # There is only one chain, hence ipTM=NaN, so use just pTM.
      ranking_confidence = ptm

    contact_probs = result['distogram']['contact_probs']
    # Compute PAE related summaries.
    _, chain_pair_pae_min, _ = confidences.chain_pair_pae(
        num_tokens=num_tokens,
        asym_ids=batch.token_features.asym_id,
        full_pae=result['full_pae'],
        mask=pae_single_mask,
    )
    chain_pair_pde_mean, chain_pair_pde_min = confidences.chain_pair_pde(
        num_tokens=num_tokens,
        asym_ids=batch.token_features.asym_id,
        full_pde=result['full_pde'],
    )
    intra_chain_single_pde, cross_chain_single_pde, _ = confidences.pde_single(
        num_tokens,
        batch.token_features.asym_id,
        result['full_pde'],
        contact_probs,
    )
    pae_metrics = confidences.pae_metrics(
        num_tokens=num_tokens,
        asym_ids=batch.token_features.asym_id,
        full_pae=result['full_pae'],
        mask=pae_single_mask,
        contact_probs=contact_probs,
        tm_adjusted_pae=result['tmscore_adjusted_pae_interface'],
    )
    ranking_confidence_pae = confidences.rank_metric(
        result['full_pae'],
        contact_probs * batch.frames.mask[:, None].astype(float),
    )
    chain_pair_iptm = _compute_chain_pair_iptm(
        num_tokens=num_tokens,
        asym_ids=batch.token_features.asym_id,
        mask=pae_single_mask,
        tm_adjusted_pae=result['tmscore_adjusted_pae_interface'],
    )
    # iptm_ichain is a vector of per-chain ptm values. iptm_ichain[0],
    # for example, is just the zeroth diagonal entry of the chain pair iptm
    # matrix:
    # [[x, , ],
    #  [ , , ],
    #  [ , , ]]]
    iptm_ichain = chain_pair_iptm.diagonal(axis1=-2, axis2=-1)
    # iptm_xchain is a vector of cross-chain interactions for each chain.
    # iptm_xchain[0], for example, is an average of chain 0's interactions with
    # other chains:
    # [[ ,x,x],
    #  [x, , ],
    #  [x, , ]]]
    iptm_xchain = confidences.get_iptm_xchain(chain_pair_iptm)

    predicted_distance_errors = result['average_pde']

    # Computing solvent accessible area with dssp can be slow for large
    # structures with lots of chains, so we parallelize the call.
    pred_structures = pred_structure.unstack()
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=min(len(pred_structures), 32)
    ) as executor:
      has_clash = list(executor.map(confidences.has_clash, pred_structures))
      fraction_disordered = list(
          executor.map(confidences.fraction_disordered, pred_structures)
      )

    for idx, pred_structure in enumerate(pred_structures):
      ranking_score = confidences.get_ranking_score(
          ptm=ptm[idx],
          iptm=iptm[idx],
          fraction_disordered_=fraction_disordered[idx],
          has_clash_=has_clash[idx],
      )
      yield InferenceResult(
          predicted_structure=pred_structure,
          numerical_data={
              'full_pde': result['full_pde'][idx, :num_tokens, :num_tokens],
              'full_pae': result['full_pae'][idx, :num_tokens, :num_tokens],
              'contact_probs': contact_probs[:num_tokens, :num_tokens],
          },
          metadata={
              'predicted_distance_error': predicted_distance_errors[idx],
              'ranking_score': ranking_score,
              'fraction_disordered': fraction_disordered[idx],
              'has_clash': has_clash[idx],
              'predicted_tm_score': ptm[idx],
              'interface_predicted_tm_score': iptm[idx],
              'chain_pair_pde_mean': chain_pair_pde_mean[idx],
              'chain_pair_pde_min': chain_pair_pde_min[idx],
              'chain_pair_pae_min': chain_pair_pae_min[idx],
              'ptm': ptm[idx],
              'iptm': iptm[idx],
              'ptm_iptm_average': ptm_iptm_average[idx],
              'intra_chain_single_pde': intra_chain_single_pde[idx],
              'cross_chain_single_pde': cross_chain_single_pde[idx],
              'pae_ichain': pae_metrics['pae_ichain'][idx],
              'pae_xchain': pae_metrics['pae_xchain'][idx],
              'ranking_confidence': ranking_confidence[idx],
              'ranking_confidence_pae': ranking_confidence_pae[idx],
              'chain_pair_iptm': chain_pair_iptm[idx],
              'iptm_ichain': iptm_ichain[idx],
              'iptm_xchain': iptm_xchain[idx],
              'token_chain_ids': chain_ids,
              'token_res_ids': res_ids,
          },
          model_id=result['__identifier__'],
          debug_outputs={},
      )
