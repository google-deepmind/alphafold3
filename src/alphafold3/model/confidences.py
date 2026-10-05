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

"""Functions for extracting and processing confidences from model outputs."""
import warnings

from absl import logging
from alphafold3 import structure
from alphafold3.constants import residue_names
from alphafold3.cpp import mkdssp
import jax.numpy as jnp
import numpy as np
from scipy import spatial


# From Sander & Rost 1994 https://doi.org/10.1002/prot.340200303
MAX_ACCESSIBLE_SURFACE_AREA = {
    'ALA': 106.0,
    'ARG': 248.0,
    'ASN': 157.0,
    'ASP': 163.0,
    'CYS': 135.0,
    'GLN': 198.0,
    'GLU': 194.0,
    'GLY': 84.0,
    'HIS': 184.0,
    'ILE': 169.0,
    'LEU': 164.0,
    'LYS': 205.0,
    'MET': 188.0,
    'PHE': 197.0,
    'PRO': 136.0,
    'SER': 130.0,
    'THR': 142.0,
    'TRP': 227.0,
    'TYR': 222.0,
    'VAL': 142.0,
}

# Weights for ranking confidence.
_IPTM_WEIGHT = 0.8
_FRACTION_DISORDERED_WEIGHT = 0.5
_CLASH_PENALIZATION_WEIGHT = 100.0


def windowed_solvent_accessible_area(cif: str, window: int = 25) -> np.ndarray:
  """Implementation of AlphaFold-RSA.

  AlphaFold-RSA defined in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9601767.

  Args:
    cif: Raw cif string.
    window: The window over which to average accessible surface area

  Returns:
    An array of size num_res that predicts disorder by using windowed solvent
    accessible surface area.
  """
  result = mkdssp.get_dssp(cif, calculate_surface_accessibility=True)
  parse_row = False
  rasa = []
  for row in result.splitlines():
    if parse_row:
      aa = row[13:14]
      if aa == '!':
        continue
      aa3 = residue_names.PROTEIN_COMMON_ONE_TO_THREE.get(aa, 'ALA')
      max_acc = MAX_ACCESSIBLE_SURFACE_AREA[aa3]
      acc = int(row[34:38])
      norm_acc = acc / max_acc
      if norm_acc > 1.0:
        norm_acc = 1.0
      rasa.append(norm_acc)
    if row.startswith('  #  RESIDUE'):
      parse_row = True

  half_w = (window - 1) // 2
  pad_rasa = np.pad(rasa, (half_w, half_w), 'reflect')
  rasa = np.convolve(pad_rasa, np.ones(window), 'valid') / window
  return rasa


def fraction_disordered(
    struc: structure.Structure, rasa_disorder_cutoff: float = 0.581
) -> float:
  """Compute fraction of protein residues that are disordered.

  Args:
    struc: A structure to compute rASA metrics on.
    rasa_disorder_cutoff: The threshold at which residues are considered
      disordered. Default value taken from
      https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9601767.

  Returns:
    The fraction of protein residues that are disordered
    (rasa > rasa_disorder_cutoff).
  """
  struc = struc.filter_to_entity_type(protein=True)
  rasa = []
  seq_rasa = {}
  for chain_id, chain_seq in struc.chain_single_letter_sequence().items():
    if chain_seq in seq_rasa:
      # We assume that identical sequences have approximately similar rasa
      # values to speed up the computation.
      rasa.extend(seq_rasa[chain_seq])
      continue
    chain_struc = struc.filter(chain_id=chain_id)
    # Rename the chain to 'A' as MKDSSP supports only single letter chain IDs.
    chain_struc = chain_struc.rename_chain_ids(new_id_by_old_id={chain_id: 'A'})
    try:
      rasa_per_residue = windowed_solvent_accessible_area(
          chain_struc.to_mmcif()
      )
      seq_rasa[chain_seq] = rasa_per_residue
      rasa.extend(rasa_per_residue)
    except (ValueError, RuntimeError) as e:
      logging.warning('%s: rasa calculation failed: %s', struc.name, e)

  if not rasa:
    return 0.0
  return np.mean(np.array(rasa) > rasa_disorder_cutoff)


def has_clash(
    struc: structure.Structure,
    cutoff_radius: float = 1.1,
    min_clashes_for_overlap: int = 100,
    min_fraction_for_overlap: float = 0.5,
) -> bool:
  """Determine whether the structure has at least one clashing chain.

  A clashing chain is defined as having greater than 100 polymer atoms within
  1.1A of another polymer atom, or having more than 50% of the chain with
  clashing atoms.

  Args:
    struc: A structure to get clash metrics for.
    cutoff_radius: atom distances under this threshold are considered a clash.
    min_clashes_for_overlap: The minimum number of atom-atom clashes for a chain
      to be considered overlapping.
    min_fraction_for_overlap: The minimum fraction of atoms within a chain that
      are clashing for the chain to be considered overlapping.

  Returns:
    True if the structure has at least one clashing chain.
  """
  struc = struc.filter_to_entity_type(protein=True, rna=True, dna=True)
  if not struc.chains:
    return False
  coords = struc.coords
  coord_kdtree = spatial.cKDTree(coords)
  clashes_per_atom = coord_kdtree.query_ball_point(
      coords, p=2.0, r=cutoff_radius
  )
  per_atom_has_clash = np.zeros(len(coords), dtype=np.int32)
  for atom_idx, clashing_indices in enumerate(clashes_per_atom):
    for clashing_idx in clashing_indices:
      if np.abs(struc.res_id[atom_idx] - struc.res_id[clashing_idx]) > 1 or (
          struc.chain_id[atom_idx] != struc.chain_id[clashing_idx]
      ):
        per_atom_has_clash[atom_idx] = True
        break
  for chain_id in struc.chains:
    mask = struc.chain_id == chain_id
    num_atoms = np.sum(mask)
    if num_atoms == 0:
      continue
    num_clashes = np.sum(per_atom_has_clash * mask)
    frac_clashes = num_clashes / num_atoms
    if (
        num_clashes > min_clashes_for_overlap
        or frac_clashes > min_fraction_for_overlap
    ):
      return True
  return False


def get_ranking_score(
    ptm: float, iptm: float, fraction_disordered_: float, has_clash_: bool
) -> float:
  # ipTM is NaN for single chain structures. Use pTM for such cases.
  if np.isnan(iptm):
    ptm_iptm_average = ptm
  else:
    ptm_iptm_average = _IPTM_WEIGHT * iptm + (1.0 - _IPTM_WEIGHT) * ptm
  return (
      ptm_iptm_average
      + _FRACTION_DISORDERED_WEIGHT * fraction_disordered_
      - _CLASH_PENALIZATION_WEIGHT * has_clash_
  )


def rank_metric(
    full_pde: jnp.ndarray | np.ndarray, contact_probs: jnp.ndarray | np.ndarray
) -> jnp.ndarray | np.ndarray:
  """Compute the metric that will be used to rank predictions, higher is better.

  Args:
    full_pde: A [num_samples, num_tokens,num_tokens] matrix of predicted
      distance errors between pairs of tokens.
    contact_probs: A [num_tokens, num_tokens] matrix consisting of the
      probability of contact (<8A) that is returned from the distogram head.

  Returns:
    A scalar that can be used to rank (higher is better).
  """
  if not isinstance(full_pde, type(contact_probs)):
    raise ValueError('full_pde and contact_probs must be of the same type.')

  if isinstance(full_pde, np.ndarray):
    sum_fn = np.sum
  elif isinstance(full_pde, jnp.ndarray):
    sum_fn = jnp.sum
  else:
    raise ValueError('full_pde must be a numpy array or a jax array.')
  # It was found that taking the contact_map weighted average was better than
  # just the predicted distance error on its own.
  return -sum_fn(full_pde * contact_probs[None, :, :], axis=(-2, -1)) / (
      sum_fn(contact_probs) + 1e-6
  )


def weighted_mean(mask, value, axis):
  return np.mean(mask * value, axis=axis) / (1e-8 + np.mean(mask, axis=axis))


def pde_single(
    num_tokens: int,
    asym_ids: np.ndarray,
    full_pde: np.ndarray,
    contact_probs: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
  """Compute 1D PDE summaries.

  Args:
    num_tokens: The number of tokens (not including padding).
    asym_ids: The asym_ids (array of shape num_tokens).
    full_pde: A [num_samples, num_tokens, num_tokens] matrix of predicted
      distance errors.
    contact_probs: A [num_tokens, num_tokens] matrix consisting of the
      probability of contact (<8A) that is returned from the distogram head.

  Returns:
    A tuple (ichain, xchain, full_chain) where:
      `ichain` is a [num_samples, num_chains] matrix where the
      value assigned to each chain is an average of the full PDE matrix over all
      its within-chain interactions, weighted by `contact_probs`.
      `xchain` is a [num_samples, num_chains] matrix where the
      value assigned to each chain is an average of the full PDE matrix over all
      its cross-chain interactions, weighted by `contact_probs`.
      `full_chain` is a [num_samples, num_tokens] matrix where the
      value assigned to each token is an average of it PDE against all tokens,
      weighted by `contact_probs`.
  """

  full_pde = full_pde[:, :num_tokens, :num_tokens]
  contact_probs = contact_probs[:num_tokens, :num_tokens]
  asym_ids = asym_ids[:num_tokens]
  unique_asym_ids = np.unique(asym_ids)
  num_chains = len(unique_asym_ids)
  num_samples = full_pde.shape[0]

  asym_ids = asym_ids[None]
  contact_probs = contact_probs[None]

  ichain = np.zeros((num_samples, num_chains))
  xchain = np.zeros((num_samples, num_chains))

  for idx, asym_id in enumerate(unique_asym_ids):
    my_asym_id = asym_ids == asym_id
    imask = my_asym_id[:, :, None] * my_asym_id[:, None, :]
    xmask = my_asym_id[:, :, None] * ~my_asym_id[:, None, :]
    imask = imask * contact_probs
    xmask = xmask * contact_probs
    ichain[:, idx] = weighted_mean(mask=imask, value=full_pde, axis=(-2, -1))
    xchain[:, idx] = weighted_mean(mask=xmask, value=full_pde, axis=(-2, -1))

  full_chain = weighted_mean(mask=contact_probs, value=full_pde, axis=(-1,))

  return ichain, xchain, full_chain


def chain_pair_pde(
    num_tokens: int, asym_ids: np.ndarray, full_pde: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
  """Compute predicted distance errors for all pairs of chains.

  Args:
    num_tokens: The number of tokens (not including padding).
    asym_ids: The asym_ids (array of shape num_tokens).
    full_pde: A [num_samples, num_tokens, num_tokens] matrix of predicted
      distance errors.

  Returns:
    chain_pair_pred_err_mean - a [num_chains, num_chains] matrix with average
      per chain-pair predicted distance error.
    chain_pair_pred_err_min - a [num_chains, num_chains] matrix with min
      per chain-pair predicted distance error.
  """
  full_pde = full_pde[:, :num_tokens, :num_tokens]
  asym_ids = asym_ids[:num_tokens]
  unique_asym_ids = np.unique(asym_ids)
  num_chains = len(unique_asym_ids)
  num_samples = full_pde.shape[0]
  chain_pair_pred_err_mean = np.zeros((num_samples, num_chains, num_chains))
  chain_pair_pred_err_min = np.zeros((num_samples, num_chains, num_chains))

  for idx1, asym_id_1 in enumerate(unique_asym_ids):
    subset = full_pde[:, asym_ids == asym_id_1, :]
    for idx2, asym_id_2 in enumerate(unique_asym_ids):
      subsubset = subset[:, :, asym_ids == asym_id_2]
      chain_pair_pred_err_mean[:, idx1, idx2] = np.mean(subsubset, axis=(1, 2))
      chain_pair_pred_err_min[:, idx1, idx2] = np.min(subsubset, axis=(1, 2))
  return chain_pair_pred_err_mean, chain_pair_pred_err_min


def weighted_nanmean(
    value: np.ndarray, mask: np.ndarray, axis: int
) -> np.ndarray:
  """Nan-mean with weighting -- empty slices return NaN."""
  assert mask.shape == value.shape
  assert not np.isnan(mask).all()

  nan_idxs = np.where(np.isnan(value))
  # Need to NaN the mask to get the correct denominator weighting.
  mask_with_nan = mask.copy()
  mask_with_nan[nan_idxs] = np.nan
  with warnings.catch_warnings():
    # Mean of empty slice is ok and should return a NaN.
    warnings.filterwarnings(action='ignore', message='Mean of empty slice')
    warnings.filterwarnings(
        action='ignore', message='invalid value encountered in (scalar )?divide'
    )
    return np.nanmean(value * mask_with_nan, axis=axis) / np.nanmean(
        mask_with_nan, axis=axis
    )


def chain_pair_pae(
    *,
    num_tokens: int,
    asym_ids: np.ndarray,
    full_pae: np.ndarray,
    mask: np.ndarray | None = None,
    contact_probs: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
  """Compute predicted errors for all pairs of chains.

  Args:
    num_tokens: The number of tokens (not including padding).
    asym_ids: The asym_ids (array of shape num_tokens).
    full_pae: A [num_samples, num_tokens, num_tokens] matrix of predicted
      errors.
    mask: A [num_tokens, num_tokens] mask matrix.
    contact_probs: A [num_tokens, num_tokens] matrix consisting of the
      probability of contact (<8A) that is returned from the distogram head.

  Returns:
    chain_pair_pred_err_mean - a [num_chains, num_chains] matrix with average
      per chain-pair predicted error.
  """
  if mask is None:
    mask = np.ones(shape=full_pae.shape[1:], dtype=bool)
  if contact_probs is None:
    contact_probs = np.ones(shape=full_pae.shape[1:], dtype=float)
  assert mask.shape == full_pae.shape[1:]

  full_pae = full_pae[:, :num_tokens, :num_tokens]
  mask = mask[:num_tokens, :num_tokens]
  asym_ids = asym_ids[:num_tokens]
  contact_probs = contact_probs[:num_tokens, :num_tokens]
  unique_asym_ids = np.unique(asym_ids)
  num_chains = len(unique_asym_ids)
  num_samples = full_pae.shape[0]
  chain_pair_pred_err_mean = np.zeros((num_samples, num_chains, num_chains))
  chain_pair_pred_err_min = np.zeros((num_samples, num_chains, num_chains))

  for idx1, asym_id_1 in enumerate(unique_asym_ids):
    subset = full_pae[:, asym_ids == asym_id_1, :]
    subset_mask = mask[asym_ids == asym_id_1, :]
    subset_contact_probs = contact_probs[asym_ids == asym_id_1, :]
    for idx2, asym_id_2 in enumerate(unique_asym_ids):
      subsubset = subset[:, :, asym_ids == asym_id_2]
      subsubset_mask = subset_mask[:, asym_ids == asym_id_2]
      subsubset_contact_probs = subset_contact_probs[:, asym_ids == asym_id_2]
      (flat_mask_idxs,) = np.where(subsubset_mask.flatten() > 0)
      flat_subsubset = subsubset.reshape([num_samples, -1])
      flat_contact_probs = subsubset_contact_probs.flatten()
      # A ligand chain will have no valid frames if it contains fewer than
      # three non-colinear atoms (e.g. a sodium ion).
      if not flat_mask_idxs.size:
        chain_pair_pred_err_mean[:, idx1, idx2] = np.nan
        chain_pair_pred_err_min[:, idx1, idx2] = np.nan
      else:
        chain_pair_pred_err_min[:, idx1, idx2] = np.min(
            flat_subsubset[:, flat_mask_idxs], axis=1
        )
        chain_pair_pred_err_mean[:, idx1, idx2] = weighted_mean(
            mask=flat_contact_probs[flat_mask_idxs],
            value=flat_subsubset[:, flat_mask_idxs],
            axis=-1,
        )
  return chain_pair_pred_err_mean, chain_pair_pred_err_min, unique_asym_ids


def reduce_chain_pair(
    *,
    chain_pair_met: np.ndarray,
    num_chain_tokens: np.ndarray,
    agg_over_col: bool,
    agg_type: str,
    weight_method: str,
) -> tuple[np.ndarray, np.ndarray]:
  """Compute 1D summaries from a chain-pair summary.

  Args:
    chain_pair_met: A [num_samples, num_chains, num_chains] aggregate matrix.
    num_chain_tokens: A [num_chains] array of number of tokens for each chain.
      Used for 'per_token' weighting.
    agg_over_col: Whether to aggregate the PAE over rows (i.e. average error
      when aligned to me) or columns (i.e. my average error when aligned to all
      others.)
    agg_type: The type of aggregation to use, 'mean' or 'min'.
    weight_method: The method to use for weighting the PAE, 'per_token' or
      'per_chain'.

  Returns:
    A tuple (ichain, xchain) where:
      `ichain` is a [num_samples, num_chains] matrix where the
      value assigned to each chain is an average of the full PAE matrix over all
      its within-chain interactions, weighted by `contact_probs`.
      `xchain` is a [num_samples, num_chains] matrix where the
      value assigned to each chain is an average of the full PAE matrix over all
      its cross-chain interactions, weighted by `contact_probs`.
  """
  num_samples, num_chains, _ = chain_pair_met.shape

  ichain = chain_pair_met.diagonal(axis1=-2, axis2=-1)

  if weight_method == 'per_chain':
    chain_weight = np.ones((num_chains,), dtype=float)
  elif weight_method == 'per_token':
    chain_weight = num_chain_tokens
  else:
    raise ValueError(f'Unknown weight method: {weight_method}')

  if agg_over_col:
    agg_axis = -1
  else:
    agg_axis = -2

  if agg_type == 'mean':
    weight = np.ones((num_samples, num_chains, num_chains), dtype=float)
    weight -= np.eye(num_chains, dtype=float)
    weight *= chain_weight[None] * chain_weight[:, None]
    xchain = weighted_nanmean(chain_pair_met, mask=weight, axis=agg_axis)
  elif agg_type == 'min':
    is_self = np.eye(num_chains)
    with warnings.catch_warnings():
      # Min over empty slice is ok and should return a NaN.
      warnings.filterwarnings('ignore', message='All-NaN slice encountered')
      xchain = np.nanmin(chain_pair_met + 1e8 * is_self, axis=agg_axis)
  else:
    raise ValueError(f'Unknown aggregation method: {agg_type}')

  return ichain, xchain


def pae_metrics(
    num_tokens: int,
    asym_ids: np.ndarray,
    full_pae: np.ndarray,
    mask: np.ndarray,
    contact_probs: np.ndarray,
    tm_adjusted_pae: np.ndarray,
):
  """PAE aggregate metrics."""
  assert mask.shape == full_pae.shape[1:]
  assert contact_probs.shape == full_pae.shape[1:]

  chain_pair_contact_weighted, _, unique_asym_ids = chain_pair_pae(
      num_tokens=num_tokens,
      asym_ids=asym_ids,
      full_pae=full_pae,
      mask=mask,
      contact_probs=contact_probs,
  )

  ret = {}
  ret['chain_pair_pae_mean'], ret['chain_pair_pae_min'], _ = chain_pair_pae(
      num_tokens=num_tokens,
      asym_ids=asym_ids,
      full_pae=full_pae,
      mask=mask,
  )
  chain_pair_iptm = np.stack(
      [
          chain_pairwise_predicted_tm_scores(
              tm_adjusted_pae=sample_tm_adjusted_pae[:num_tokens],
              asym_id=asym_ids[:num_tokens],
              pair_mask=mask[:num_tokens, :num_tokens],
          )
          for sample_tm_adjusted_pae in tm_adjusted_pae
      ],
      axis=0,
  )

  num_chain_tokens = np.array(
      [sum(asym_ids == asym_id) for asym_id in unique_asym_ids]
  )

  def reduce_chain_pair_fn(chain_pair: np.ndarray):
    def inner(agg_over_col):
      ichain_pae, xchain_pae = reduce_chain_pair(
          num_chain_tokens=num_chain_tokens,
          chain_pair_met=chain_pair,
          agg_over_col=agg_over_col,
          agg_type='mean',
          weight_method='per_chain',
      )
      return ichain_pae, xchain_pae

    ichain, xchain_row_agg = inner(False)
    _, xchain_col_agg = inner(True)
    with warnings.catch_warnings():
      # Mean of empty slice is ok and should return a NaN.
      warnings.filterwarnings(action='ignore', message='Mean of empty slice')
      xchain = np.nanmean(
          np.stack([xchain_row_agg, xchain_col_agg], axis=0), axis=0
      )
    return ichain, xchain

  pae_ichain, pae_xchain = reduce_chain_pair_fn(chain_pair_contact_weighted)
  iptm_ichain, iptm_xchain = reduce_chain_pair_fn(chain_pair_iptm)

  ret.update({
      'chain_pair_iptm': chain_pair_iptm,
      'iptm_ichain': iptm_ichain,
      'iptm_xchain': iptm_xchain,
      'pae_ichain': pae_ichain,
      'pae_xchain': pae_xchain,
  })

  return ret


def get_iptm_xchain(chain_pair_iptm: np.ndarray) -> np.ndarray:
  """Cross chain aggregate ipTM."""
  num_samples, num_chains, _ = chain_pair_iptm.shape
  weight = np.ones((num_samples, num_chains, num_chains), dtype=float)
  weight -= np.eye(num_chains, dtype=float)
  xchain_row_agg = weighted_nanmean(chain_pair_iptm, mask=weight, axis=-2)
  xchain_col_agg = weighted_nanmean(chain_pair_iptm, mask=weight, axis=-1)
  with warnings.catch_warnings():
    # Mean of empty slice is ok and should return a NaN.
    warnings.filterwarnings(action='ignore', message='Mean of empty slice')
    iptm_xchain = np.nanmean(
        np.stack([xchain_row_agg, xchain_col_agg], axis=0), axis=0
    )
  return iptm_xchain


def predicted_tm_score(
    tm_adjusted_pae: np.ndarray,
    pair_mask: np.ndarray,
    asym_id: np.ndarray,
    interface: bool = False,
) -> float:
  """Computes predicted TM alignment or predicted interface TM alignment score.

  Args:
    tm_adjusted_pae: [num_res, num_res] Relevant tensor for computing TMScore
      values.
    pair_mask: A [num_res, num_res] mask. The TM score will only aggregate over
      masked-on entries.
    asym_id: [num_res] asymmetric unit ID (the chain ID). Only needed for ipTM
      calculation, i.e. when interface=True.
    interface: If True, the interface predicted TM score is computed. If False,
      the predicted TM score without any residue pair restrictions is computed.

  Returns:
   score: pTM or ipTM score.
  """
  num_tokens, _ = tm_adjusted_pae.shape
  if tm_adjusted_pae.shape != (num_tokens, num_tokens):
    raise ValueError(
        f'Bad tm_adjusted_pae shape, expected ({num_tokens, num_tokens}), got '
        f'{tm_adjusted_pae.shape}.'
    )

  if pair_mask.shape != (num_tokens, num_tokens):
    raise ValueError(
        f'Bad pair_mask shape, expected ({num_tokens, num_tokens}), got '
        f'{pair_mask.shape}.'
    )
  if pair_mask.dtype != bool:
    raise TypeError(f'Bad pair mask type, expected bool, got {pair_mask.dtype}')
  if asym_id.shape[0] != num_tokens:
    raise ValueError(
        f'Bad asym_id shape, expected ({num_tokens},), got {asym_id.shape}.'
    )

  # Create pair mask.
  if interface:
    pair_mask = pair_mask * (asym_id[:, None] != asym_id[None, :])

  # Ions and other ligands with colinear atoms have ill-defined frames.
  if pair_mask.sum() == 0:
    return np.nan

  normed_residue_mask = pair_mask / (
      1e-8 + np.sum(pair_mask, axis=-1, keepdims=True)
  )
  per_alignment = np.sum(tm_adjusted_pae * normed_residue_mask, axis=-1)
  return per_alignment.max()


def chain_pairwise_predicted_tm_scores(
    tm_adjusted_pae: np.ndarray,
    pair_mask: np.ndarray,
    asym_id: np.ndarray,
) -> np.ndarray:
  """Compute predicted TM (pTM) between each pair of chains independently.

  Args:
    tm_adjusted_pae: [num_res, num_res] Relevant tensor for computing TMScore
      values.
    pair_mask: A [num_res, num_res] mask specifying which frames are valid.
      Invalid frames can be the result of chains with not enough atoms (e.g.
      ions).
    asym_id: [num_res] asymmetric unit ID (the chain ID).

  Returns:
    A [num_chains, num_chains] matrix, where row i, column j indicates the
    predicted TM-score for the interface between chain i and chain j.
  """
  unique_chains = list(np.unique(asym_id))
  num_chains = len(unique_chains)
  all_pairs_iptms = np.zeros((num_chains, num_chains))
  for i, chain_i in enumerate(unique_chains):
    chain_i_mask = asym_id == chain_i
    for j, chain_j in enumerate(unique_chains[i:]):
      chain_j_mask = asym_id == chain_j
      mask = chain_i_mask | chain_j_mask
      (indices,) = np.where(mask)
      is_interface = chain_i != chain_j
      indices = np.ix_(indices, indices)
      iptm = predicted_tm_score(
          tm_adjusted_pae=tm_adjusted_pae[indices],
          pair_mask=pair_mask[indices],
          asym_id=asym_id[mask],
          interface=is_interface,
      )
      all_pairs_iptms[i, i + j] = iptm
      all_pairs_iptms[i + j, i] = iptm
  return all_pairs_iptms
