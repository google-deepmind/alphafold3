# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Data-side of the input features processing."""

import dataclasses
import datetime
import itertools
from typing import Any, Self, TypeAlias

from absl import logging
from alphafold3 import structure
from alphafold3.common import folding_input
from alphafold3.constants import chemical_components
from alphafold3.constants import mmcif_names
from alphafold3.constants import periodic_table
from alphafold3.constants import residue_names
from alphafold3.cpp import cif_dict
from alphafold3.data import msa as msa_module
from alphafold3.data import templates
from alphafold3.data.tools import rdkit_utils
from alphafold3.model import data3
from alphafold3.model import data_constants
from alphafold3.model import merging_features
from alphafold3.model import msa_pairing
from alphafold3.model.atom_layout import atom_layout
from alphafold3.structure import chemical_components as struc_chem_comps
import chex
import jax.numpy as jnp
import numpy as np
from rdkit import Chem


xnp_ndarray: TypeAlias = np.ndarray | jnp.ndarray  # pylint: disable=invalid-name
BatchDict: TypeAlias = dict[str, xnp_ndarray]

_STANDARD_RESIDUES = frozenset({
    *residue_names.PROTEIN_TYPES_WITH_UNKNOWN,
    *residue_names.NUCLEIC_TYPES_WITH_2_UNKS,
})


@dataclasses.dataclass(frozen=True)
class PaddingShapes:
  num_tokens: int
  msa_size: int
  num_chains: int
  num_templates: int
  num_atoms: int


def _pad_to(
    arr: np.ndarray, shape: tuple[int | None, ...], **kwargs
) -> np.ndarray:
  """Pads an array to a given shape. Wrapper around np.pad().

  Args:
    arr: numpy array to pad
    shape: target shape, use None for axes that should stay the same
    **kwargs: additional args for np.pad, e.g. constant_values=-1

  Returns:
    the padded array

  Raises:
    ValueError if arr and shape have a different number of axes.
  """
  if arr.ndim != len(shape):
    raise ValueError(
        f'arr and shape have different number of axes. {arr.shape=}, {shape=}'
    )

  num_pad = []
  for axis, width in enumerate(shape):
    if width is None:
      num_pad.append((0, 0))
    else:
      if width >= arr.shape[axis]:
        num_pad.append((0, width - arr.shape[axis]))
      else:
        raise ValueError(
            f'Can not pad to a smaller shape. {arr.shape=}, {shape=}'
        )
  padded_arr = np.pad(arr, pad_width=num_pad, **kwargs)
  return padded_arr


def _unwrap(obj):
  """Unwrap an object from a zero-dim np.ndarray."""
  if isinstance(obj, np.ndarray) and obj.ndim == 0:
    return obj.item()
  else:
    return obj


@chex.dataclass(mappable_dataclass=False, frozen=True)
class Chains:
  chain_id: np.ndarray
  asym_id: np.ndarray
  entity_id: np.ndarray
  sym_id: np.ndarray


def _compute_asym_entity_and_sym_id(
    all_tokens: atom_layout.AtomLayout,
) -> Chains:
  """Compute asym_id, entity_id and sym_id.

  Args:
    all_tokens: atom layout containing a representative atom for each token.

  Returns:
    A Chains object
  """

  # Find identical sequences and assign entity_id and sym_id to every chain.
  seq_to_entity_id_sym_id = {}
  seen_chain_ids = set()
  chain_ids = []
  asym_ids = []
  entity_ids = []
  sym_ids = []
  for chain_id in all_tokens.chain_id:
    if chain_id not in seen_chain_ids:
      asym_id = len(seen_chain_ids) + 1
      seen_chain_ids.add(chain_id)
      seq = ','.join(all_tokens.res_name[all_tokens.chain_id == chain_id])
      if seq not in seq_to_entity_id_sym_id:
        entity_id = len(seq_to_entity_id_sym_id) + 1
        sym_id = 1
      else:
        entity_id, sym_id = seq_to_entity_id_sym_id[seq]
        sym_id += 1
      seq_to_entity_id_sym_id[seq] = (entity_id, sym_id)

      chain_ids.append(chain_id)
      asym_ids.append(asym_id)
      entity_ids.append(entity_id)
      sym_ids.append(sym_id)

  return Chains(
      chain_id=np.array(chain_ids),
      asym_id=np.array(asym_ids),
      entity_id=np.array(entity_ids),
      sym_id=np.array(sym_ids),
  )


def tokenizer(
    flat_output_layout: atom_layout.AtomLayout,
    ccd: chemical_components.Ccd,
    max_atoms_per_token: int,
    flatten_non_standard_residues: bool,
    logging_name: str,
) -> tuple[atom_layout.AtomLayout, atom_layout.AtomLayout, np.ndarray]:
  """Maps a flat atom layout to tokens for evoformer.

  Creates the evoformer tokens as one token per polymer residue and one token
  per ligand atom. The tokens are represented as AtomLayouts all_tokens
  (1 representative atom per token) atoms per residue, and
  all_token_atoms_layout (num_tokens, max_atoms_per_token). The atoms in a
  residue token use the layout of the corresponding CCD entry

  Args:
    flat_output_layout: flat AtomLayout containing all atoms that the model
      wants to predict.
    ccd: The chemical components dictionary.
    max_atoms_per_token: number of slots per token.
    flatten_non_standard_residues: whether to flatten non-standard residues,
      i.e. whether to use one token per atom for non-standard residues.
    logging_name: logging name for debugging (usually the mmcif_id).

  Returns:
    A tuple (all_tokens, all_tokens_atoms_layout) with
      all_tokens: AtomLayout shape (num_tokens,) containing one representative
        atom per token.
      all_token_atoms_layout: AtomLayout with shape
        (num_tokens, max_atoms_per_token) containing all atoms per token.
      standard_token_idxs: The token index that each token would have if not
        flattening non standard resiudes.
  """
  # Select  the representative atom for each token.
  token_idxs = []
  single_atom_token = []
  standard_token_idxs = []
  current_standard_token_id = 0
  # Iterate over residues, and provide a group_iter over the atoms of each
  # residue.
  for key, group_iter in itertools.groupby(
      zip(
          flat_output_layout.chain_type,
          flat_output_layout.chain_id,
          flat_output_layout.res_id,
          flat_output_layout.res_name,
          flat_output_layout.atom_name,
          np.arange(flat_output_layout.shape[0]),
      ),
      key=lambda x: x[:3],
  ):

    # Get chain type and chain id of this residue
    chain_type, chain_id, _ = key

    # Get names and global idxs for all atoms of this residue
    _, _, _, res_names, atom_names, idxs = zip(*group_iter)

    # As of March 2023, all OTHER CHAINs in pdb are artificial nucleics.
    is_nucleic_backbone = (
        chain_type in mmcif_names.NUCLEIC_ACID_CHAIN_TYPES
        or chain_type == mmcif_names.OTHER_CHAIN
    )
    if chain_type in mmcif_names.PEPTIDE_CHAIN_TYPES:
      res_name = res_names[0]
      if (
          flatten_non_standard_residues
          and res_name not in residue_names.PROTEIN_TYPES_WITH_UNKNOWN
          and res_name != residue_names.MSE
      ):
        # For non-standard protein residues take all atoms.
        # NOTE: This may get very large if we include hydrogens.
        token_idxs.extend(idxs)
        single_atom_token += [True] * len(idxs)
        standard_token_idxs.extend([current_standard_token_id] * len(idxs))
      else:
        # For standard protein residues take 'CA' if it exists, else first atom.
        if 'CA' in atom_names:
          token_idxs.append(idxs[atom_names.index('CA')])
        else:
          token_idxs.append(idxs[0])
        single_atom_token += [False]
        standard_token_idxs.append(current_standard_token_id)
      current_standard_token_id += 1
    elif is_nucleic_backbone:
      res_name = res_names[0]
      if (
          flatten_non_standard_residues
          and res_name not in residue_names.NUCLEIC_TYPES_WITH_2_UNKS
      ):
        # For non-standard nucleic residues take all atoms.
        token_idxs.extend(idxs)
        single_atom_token += [True] * len(idxs)
        standard_token_idxs.extend([current_standard_token_id] * len(idxs))
      else:
        # For standard nucleic residues take C1' if it exists, else first atom.
        if "C1'" in atom_names:
          token_idxs.append(idxs[atom_names.index("C1'")])
        else:
          token_idxs.append(idxs[0])
        single_atom_token += [False]
        standard_token_idxs.append(current_standard_token_id)
      current_standard_token_id += 1
    elif chain_type in mmcif_names.NON_POLYMER_CHAIN_TYPES:
      # For non-polymers take all atoms
      token_idxs.extend(idxs)
      single_atom_token += [True] * len(idxs)
      standard_token_idxs.extend([current_standard_token_id] * len(idxs))
      current_standard_token_id += len(idxs)
    else:
      # Chain type that we don't handle yet.
      logging.warning(
          '%s: ignoring chain %s with chain type %s.',
          logging_name,
          chain_id,
          chain_type,
      )

  assert len(token_idxs) == len(single_atom_token)
  assert len(token_idxs) == len(standard_token_idxs)
  standard_token_idxs = np.array(standard_token_idxs, dtype=np.int32)

  # Create the list of all tokens, represented as a flat AtomLayout with 1
  # representative atom per token.
  all_tokens = flat_output_layout[token_idxs]

  # Create the 2D atoms_per_token layout
  num_tokens = all_tokens.shape[0]

  # Target lists.
  target_atom_names = []
  target_atom_elements = []
  target_res_ids = []
  target_res_names = []
  target_chain_ids = []
  target_chain_types = []

  # uids of all atoms in the flat layout, to check whether the dense atoms
  # exist -- This is necessary for terminal atoms (e.g. 'OP3' or 'OXT')
  all_atoms_uids = set(
      zip(
          flat_output_layout.chain_id,
          flat_output_layout.res_id,
          flat_output_layout.atom_name,
      )
  )

  for idx, single_atom in enumerate(single_atom_token):
    if not single_atom:
      # Standard protein and nucleic residues have many atoms per token
      chain_id = all_tokens.chain_id[idx]
      res_id = all_tokens.res_id[idx]
      res_name = all_tokens.res_name[idx]
      atom_names = []
      atom_elements = []

      res_atoms = struc_chem_comps.get_all_atoms_in_entry(
          ccd=ccd, res_name=res_name
      )
      atom_names_elements = list(
          zip(
              res_atoms['_chem_comp_atom.atom_id'],
              res_atoms['_chem_comp_atom.type_symbol'],
              strict=True,
          )
      )

      for atom_name, atom_element in atom_names_elements:
        # Remove hydrogens if they are not in flat layout.
        if atom_element in ['H', 'D'] and (
            (chain_id, res_id, atom_name) not in all_atoms_uids
        ):
          continue
        elif (chain_id, res_id, atom_name) in all_atoms_uids:
          atom_names.append(atom_name)
          atom_elements.append(atom_element)
        # Leave spaces for OXT etc.
        else:
          atom_names.append('')
          atom_elements.append('')

      if len(atom_names) > max_atoms_per_token:
        logging.warning(
            'Atom list for chain %s '
            'residue %s %s is too long and will be truncated: '
            '%s to the max atoms limit %s. Dropped atoms: %s',
            chain_id,
            res_id,
            res_name,
            len(atom_names),
            max_atoms_per_token,
            list(
                zip(
                    atom_names[max_atoms_per_token:],
                    atom_elements[max_atoms_per_token:],
                    strict=True,
                )
            ),
        )
        atom_names = atom_names[:max_atoms_per_token]
        atom_elements = atom_elements[:max_atoms_per_token]

      num_pad = max_atoms_per_token - len(atom_names)
      atom_names.extend([''] * num_pad)
      atom_elements.extend([''] * num_pad)

    else:
      # ligands have only 1 atom per token
      padding = [''] * (max_atoms_per_token - 1)
      atom_names = [all_tokens.atom_name[idx]] + padding
      atom_elements = [all_tokens.atom_element[idx]] + padding

    # Append the atoms to the target lists.
    target_atom_names.append(atom_names)
    target_atom_elements.append(atom_elements)
    target_res_names.append([all_tokens.res_name[idx]] * max_atoms_per_token)
    target_res_ids.append([all_tokens.res_id[idx]] * max_atoms_per_token)
    target_chain_ids.append([all_tokens.chain_id[idx]] * max_atoms_per_token)
    target_chain_types.append(
        [all_tokens.chain_type[idx]] * max_atoms_per_token
    )

  # Make sure to get the right shape also for 0 tokens
  trg_shape = (num_tokens, max_atoms_per_token)
  all_token_atoms_layout = atom_layout.AtomLayout(
      atom_name=np.array(target_atom_names, dtype=object).reshape(trg_shape),
      atom_element=np.array(target_atom_elements, dtype=object).reshape(
          trg_shape
      ),
      res_name=np.array(target_res_names, dtype=object).reshape(trg_shape),
      res_id=np.array(target_res_ids, dtype=int).reshape(trg_shape),
      chain_id=np.array(target_chain_ids, dtype=object).reshape(trg_shape),
      chain_type=np.array(target_chain_types, dtype=object).reshape(trg_shape),
  )

  return all_tokens, all_token_atoms_layout, standard_token_idxs


@chex.dataclass(mappable_dataclass=False, frozen=True)
class MSA:
  """Dataclass containing MSA."""

  rows: xnp_ndarray
  mask: xnp_ndarray
  deletion_matrix: xnp_ndarray
  # Occurrence of each residue type along the sequence, averaged over MSA rows.
  profile: xnp_ndarray
  # Occurrence of deletions along the sequence, averaged over MSA rows.
  deletion_mean: xnp_ndarray
  # Number of MSA alignments.
  num_alignments: xnp_ndarray

  @classmethod
  def compute_features(
      cls,
      *,
      all_tokens: atom_layout.AtomLayout,
      standard_token_idxs: np.ndarray,
      padding_shapes: PaddingShapes,
      fold_input: folding_input.Input,
      logging_name: str,
      max_paired_sequence_per_species: int,
  ) -> Self:
    """Compute the msa features."""
    seen_entities = {}

    substruct = atom_layout.make_structure(
        flat_layout=all_tokens,
        atom_coords=np.zeros(all_tokens.shape + (3,)),
        name=logging_name,
    )
    prot = substruct.filter_to_entity_type(protein=True)
    num_unique_chains = len(set(prot.chain_single_letter_sequence().values()))
    need_msa_pairing = num_unique_chains > 1

    np_chains_list = []
    input_chains_by_id = {chain.id: chain for chain in fold_input.chains}
    nonempty_chain_ids = set(all_tokens.chain_id)
    for asym_id, chain_info in enumerate(substruct.iter_chains(), start=1):
      b_chain_id = chain_info['chain_id']
      chain_type = chain_info['chain_type']
      chain = input_chains_by_id[b_chain_id]

      # Generalised "sequence" for ligands (can't trust residue name)
      chain_tokens = all_tokens[all_tokens.chain_id == b_chain_id]
      assert chain_tokens.res_name is not None
      three_letter_sequence = ','.join(chain_tokens.res_name.tolist())
      chain_num_tokens = len(chain_tokens.atom_name)
      if chain_type in mmcif_names.POLYMER_CHAIN_TYPES:
        sequence = substruct.chain_single_letter_sequence()[b_chain_id]
        if chain_type in mmcif_names.NUCLEIC_ACID_CHAIN_TYPES:
          # Only allow nucleic residue types for nucleic chains (can have some
          # protein residues in e.g. tRNA, but that causes MSA search failures).
          # Replace non nucleic residue types by UNK_NUCLEIC.
          nucleic_types_one_letter = (
              residue_names.DNA_TYPES_ONE_LETTER
              + residue_names.RNA_TYPES_ONE_LETTER_WITH_UNKNOWN
          )
          sequence = ''.join([
              base
              if base in nucleic_types_one_letter
              else residue_names.UNK_NUCLEIC_ONE_LETTER
              for base in sequence
          ])
      else:
        sequence = 'X' * chain_num_tokens

      skip_chain = (
          chain_type not in mmcif_names.STANDARD_POLYMER_CHAIN_TYPES
          or len(sequence) <= 4
          or b_chain_id not in nonempty_chain_ids
      )
      if three_letter_sequence in seen_entities:
        entity_id = seen_entities[three_letter_sequence]
      else:
        entity_id = len(seen_entities) + 1

      if chain_type in mmcif_names.STANDARD_POLYMER_CHAIN_TYPES:
        unpaired_a3m = ''
        paired_a3m = ''
        if not skip_chain:
          if need_msa_pairing and isinstance(chain, folding_input.ProteinChain):
            paired_a3m = chain.paired_msa
          if isinstance(
              chain, folding_input.RnaChain | folding_input.ProteinChain
          ):
            unpaired_a3m = chain.unpaired_msa
        # If we generated the MSA ourselves, it is already deduplicated. If it
        # is user-provided, keep it as is to prevent destroying desired pairing.
        unpaired_msa = msa_module.Msa.from_a3m(
            query_sequence=sequence,
            chain_poly_type=chain_type,
            a3m=unpaired_a3m,
            deduplicate=False,
        )

        paired_msa = msa_module.Msa.from_a3m(
            query_sequence=sequence,
            chain_poly_type=chain_type,
            a3m=paired_a3m,
            deduplicate=False,
        )
      else:
        unpaired_msa = msa_module.Msa.from_empty(
            query_sequence='-' * len(sequence),
            chain_poly_type=mmcif_names.PROTEIN_CHAIN,
        )
        paired_msa = msa_module.Msa.from_empty(
            query_sequence='-' * len(sequence),
            chain_poly_type=mmcif_names.PROTEIN_CHAIN,
        )

      msa_features = unpaired_msa.featurize()
      all_seqs_msa_features = paired_msa.featurize()

      msa_features = data3.fix_features(msa_features)
      all_seqs_msa_features = data3.fix_features(all_seqs_msa_features)

      msa_features = msa_features | {
          f'{k}_all_seq': v for k, v in all_seqs_msa_features.items()
      }
      feats = msa_features
      feats['chain_id'] = b_chain_id
      feats['asym_id'] = np.full(chain_num_tokens, asym_id)
      feats['entity_id'] = entity_id
      np_chains_list.append(feats)

    # Add profile features to each chain.
    for chain in np_chains_list:
      chain.update(
          data3.get_profile_features(chain['msa'], chain['deletion_matrix'])
      )

    # Allow 50% of the MSA to come from MSA pairing.
    max_paired_sequences = padding_shapes.msa_size // 2
    if need_msa_pairing:
      np_chains_list = list(map(dict, np_chains_list))
      np_chains_list = msa_pairing.create_paired_features(
          np_chains_list,
          max_paired_sequences=max_paired_sequences,
          nonempty_chain_ids=nonempty_chain_ids,
          max_hits_per_species=max_paired_sequence_per_species,
      )
      np_chains_list = msa_pairing.deduplicate_unpaired_sequences(
          np_chains_list
      )

    # Remove all gapped rows from all seqs.
    nonempty_asym_ids = []
    for chain in np_chains_list:
      if chain['chain_id'] in nonempty_chain_ids:
        nonempty_asym_ids.append(chain['asym_id'][0])
    if 'msa_all_seq' in np_chains_list[0]:
      np_chains_list = msa_pairing.remove_all_gapped_rows_from_all_seqs(
          np_chains_list, asym_ids=nonempty_asym_ids
      )

    # Crop MSA rows.
    cropped_chains_list = []
    for chain in np_chains_list:
      unpaired_msa_size, paired_msa_size = (
          msa_pairing.choose_paired_unpaired_msa_crop_sizes(
              unpaired_msa=chain['msa'],
              paired_msa=chain.get('msa_all_seq'),
              total_msa_crop_size=padding_shapes.msa_size,
              max_paired_sequences=max_paired_sequences,
          )
      )
      cropped_chain = {
          'asym_id': chain['asym_id'],
          'chain_id': chain['chain_id'],
          'profile': chain['profile'],
          'deletion_mean': chain['deletion_mean'],
      }
      for feat in data_constants.NUM_SEQ_NUM_RES_MSA_FEATURES:
        if feat in chain:
          cropped_chain[feat] = chain[feat][:unpaired_msa_size]
        if feat + '_all_seq' in chain:
          cropped_chain[feat + '_all_seq'] = chain[feat + '_all_seq'][
              :paired_msa_size
          ]
      cropped_chains_list.append(cropped_chain)

    # Merge Chains.
    # Make sure the chain order is unaltered before slicing with tokens.
    curr_chain_order = [chain['chain_id'] for chain in cropped_chains_list]
    orig_chain_order = [chain['chain_id'] for chain in substruct.iter_chains()]
    assert curr_chain_order == orig_chain_order
    np_example = {
        'asym_id': np.concatenate(
            [c['asym_id'] for c in cropped_chains_list], axis=0
        ),
    }
    for feature in data_constants.NUM_SEQ_NUM_RES_MSA_FEATURES:
      for feat in [feature, feature + '_all_seq']:
        if feat in cropped_chains_list[0]:
          np_example[feat] = merging_features.merge_msa_features(
              feat, cropped_chains_list
          )
    for feature in ['profile', 'deletion_mean']:
      feature_list = [c[feature] for c in cropped_chains_list]
      np_example[feature] = np.concatenate(feature_list, axis=0)

    # Crop MSA rows to maximum size given by chains participating in the crop.
    max_allowed_unpaired = max([
        len(chain['msa'])
        for chain in cropped_chains_list
        if chain['asym_id'][0] in nonempty_asym_ids
    ])
    np_example['msa'] = np_example['msa'][:max_allowed_unpaired]
    if 'msa_all_seq' in np_example:
      max_allowed_paired = max([
          len(chain['msa_all_seq'])
          for chain in cropped_chains_list
          if chain['asym_id'][0] in nonempty_asym_ids
      ])
      np_example['msa_all_seq'] = np_example['msa_all_seq'][:max_allowed_paired]

    np_example = merging_features.merge_paired_and_unpaired_msa(np_example)

    # Crop MSA residues. Need to use the standard token indices, since msa does
    # not expand non-standard residues. This means that for expanded residues,
    # we get repeated msa columns.
    new_cropping_idxs = standard_token_idxs
    for feature in data_constants.NUM_SEQ_NUM_RES_MSA_FEATURES:
      if feature in np_example:
        np_example[feature] = np_example[feature][:, new_cropping_idxs].copy()
    for feature in ['profile', 'deletion_mean']:
      np_example[feature] = np_example[feature][new_cropping_idxs]

    # Make MSA mask.
    np_example['msa_mask'] = np.ones_like(np_example['msa'], dtype=np.float32)

    # Count MSA size before padding.
    num_alignments = np_example['msa'].shape[0]

    # Pad:
    msa_size, num_tokens = padding_shapes.msa_size, padding_shapes.num_tokens

    def safe_cast_int8(x):
      return np.clip(x, np.iinfo(np.int8).min, np.iinfo(np.int8).max).astype(
          np.int8
      )

    return MSA(
        rows=_pad_to(safe_cast_int8(np_example['msa']), (msa_size, num_tokens)),
        mask=_pad_to(
            np_example['msa_mask'].astype(bool), (msa_size, num_tokens)
        ),
        # deletion_matrix may be out of int8 range, but we mostly care about
        # small values since we arctan it in the model.
        deletion_matrix=_pad_to(
            safe_cast_int8(np_example['deletion_matrix']),
            (msa_size, num_tokens),
        ),
        profile=_pad_to(np_example['profile'], (num_tokens, None)),
        deletion_mean=_pad_to(np_example['deletion_mean'], (num_tokens,)),
        num_alignments=np.array(num_alignments, dtype=np.int32),
    )

  def index_msa_rows(self, indices: xnp_ndarray) -> Self:
    assert indices.ndim == 1

    return MSA(
        rows=self.rows[indices, :],
        mask=self.mask[indices, :],
        deletion_matrix=self.deletion_matrix[indices, :],
        profile=self.profile,
        deletion_mean=self.deletion_mean,
        num_alignments=self.num_alignments,
    )

  @classmethod
  def from_data_dict(cls, batch: BatchDict) -> Self:
    output = cls(
        rows=batch['msa'],
        mask=batch['msa_mask'],
        deletion_matrix=batch['deletion_matrix'],
        profile=batch['profile'],
        deletion_mean=batch['deletion_mean'],
        num_alignments=batch['num_alignments'],
    )
    return output

  def as_data_dict(self) -> BatchDict:
    return {
        'msa': self.rows,
        'msa_mask': self.mask,
        'deletion_matrix': self.deletion_matrix,
        'profile': self.profile,
        'deletion_mean': self.deletion_mean,
        'num_alignments': self.num_alignments,
    }


@chex.dataclass(mappable_dataclass=False, frozen=True)
class Templates:
  """Dataclass containing templates."""

  # aatype of templates, int32 w shape [num_templates, num_res]
  aatype: xnp_ndarray
  # atom positions of templates, float32 w shape [num_templates, num_res, 24, 3]
  atom_positions: xnp_ndarray
  # atom mask of templates, bool w shape [num_templates, num_res, 24]
  atom_mask: xnp_ndarray

  @classmethod
  def compute_features(
      cls,
      all_tokens: atom_layout.AtomLayout,
      standard_token_idxs: np.ndarray,
      padding_shapes: PaddingShapes,
      fold_input: folding_input.Input,
      max_templates: int,
      logging_name: str,
  ) -> Self:
    """Compute the template features."""

    seen_entities = {}
    polymer_entity_features = {True: {}, False: {}}

    substruct = atom_layout.make_structure(
        flat_layout=all_tokens,
        atom_coords=np.zeros(all_tokens.shape + (3,)),
        name=logging_name,
    )
    np_chains_list = []

    input_chains_by_id = {chain.id: chain for chain in fold_input.chains}

    nonempty_chain_ids = set(all_tokens.chain_id)
    for chain_info in substruct.iter_chains():
      chain_id = chain_info['chain_id']
      chain_type = chain_info['chain_type']
      chain = input_chains_by_id[chain_id]

      # Generalised "sequence" for ligands (can't trust residue name)
      chain_tokens = all_tokens[all_tokens.chain_id == chain_id]
      assert chain_tokens.res_name is not None
      three_letter_sequence = ','.join(chain_tokens.res_name.tolist())
      chain_num_tokens = len(chain_tokens.atom_name)

      # Don't compute features for chains not included in the crop, or ligands.
      skip_chain = (
          chain_type != mmcif_names.PROTEIN_CHAIN
          or chain_num_tokens <= 4  # not cache filled
          or chain_id not in nonempty_chain_ids
      )

      if three_letter_sequence in seen_entities:
        entity_id = seen_entities[three_letter_sequence]
      else:
        entity_id = len(seen_entities) + 1

      if entity_id not in polymer_entity_features[skip_chain]:
        if skip_chain:
          template_features = data3.empty_template_features(chain_num_tokens)
        else:
          assert isinstance(chain, folding_input.ProteinChain)

          sorted_features = []
          for template in chain.templates:
            struc = structure.from_mmcif(
                template.mmcif,
                fix_mse_residues=True,
                fix_arginines=True,
                include_bonds=False,
                include_water=False,
                include_other=True,  # For non-standard polymer chains.
            )
            hit_features = templates.get_polymer_features(
                chain=struc,
                chain_poly_type=mmcif_names.PROTEIN_CHAIN,
                query_sequence_length=len(chain.sequence),
                query_to_hit_mapping=dict(template.query_to_template_map),
            )
            sorted_features.append(hit_features)

          template_features = templates.package_template_features(
              hit_features=sorted_features,
              include_ligand_features=False,
          )

          template_features = data3.fix_template_features(
              template_features=template_features, num_res=len(chain.sequence)
          )

        template_features = _reduce_template_features(
            template_features, max_templates
        )
        polymer_entity_features[skip_chain][entity_id] = template_features

      seen_entities[three_letter_sequence] = entity_id
      feats = polymer_entity_features[skip_chain][entity_id].copy()
      feats['chain_id'] = chain_id
      np_chains_list.append(feats)

    # We pad the num_templates dimension before merging, so that different
    # chains can be concatenated on the num_res dimension.  Masking will be
    # applied so that each chains templates can't see each other.
    for chain in np_chains_list:
      chain['template_aatype'] = _pad_to(
          chain['template_aatype'], (max_templates, None)
      )
      chain['template_atom_positions'] = _pad_to(
          chain['template_atom_positions'], (max_templates, None, None, None)
      )
      chain['template_atom_mask'] = _pad_to(
          chain['template_atom_mask'], (max_templates, None, None)
      )

    # Merge on token dimension.
    np_example = {
        ft: np.concatenate([c[ft] for c in np_chains_list], axis=1)
        for ft in np_chains_list[0]
        if ft in data_constants.TEMPLATE_FEATURES
    }

    # Crop template data. Need to use the standard token indices, since msa does
    # not expand non-standard residues. This means that for expanded residues,
    # we get repeated template information.
    for feature_name, v in np_example.items():
      np_example[feature_name] = v[:max_templates, standard_token_idxs, ...]

    # Pad along the token dimension.
    templates_features = Templates(
        aatype=_pad_to(
            np_example['template_aatype'], (None, padding_shapes.num_tokens)
        ),
        atom_positions=_pad_to(
            np_example['template_atom_positions'],
            (None, padding_shapes.num_tokens, None, None),
        ),
        atom_mask=_pad_to(
            np_example['template_atom_mask'].astype(bool),
            (None, padding_shapes.num_tokens, None),
        ),
    )
    return templates_features

  @classmethod
  def from_data_dict(cls, batch: BatchDict) -> Self:
    """Make Template from batch dictionary."""
    return cls(
        aatype=batch['template_aatype'],
        atom_positions=batch['template_atom_positions'],
        atom_mask=batch['template_atom_mask'],
    )

  def as_data_dict(self) -> BatchDict:
    return {
        'template_aatype': self.aatype,
        'template_atom_positions': self.atom_positions,
        'template_atom_mask': self.atom_mask,
    }


def _reduce_template_features(
    template_features: data3.FeatureDict,
    max_templates: int,
) -> data3.FeatureDict:
  """Reduces template features to max num templates and defined feature set."""
  num_templates = template_features['template_aatype'].shape[0]
  template_keep_mask = np.arange(num_templates) < max_templates
  template_fields = data_constants.TEMPLATE_FEATURES + (
      'template_release_timestamp',
  )
  template_features = {
      k: v[template_keep_mask]
      for k, v in template_features.items()
      if k in template_fields
  }
  return template_features


@chex.dataclass(mappable_dataclass=False, frozen=True)
class TokenFeatures:
  """Dataclass containing features for tokens."""

  residue_index: xnp_ndarray
  token_index: xnp_ndarray
  aatype: xnp_ndarray
  mask: xnp_ndarray
  seq_length: xnp_ndarray

  # Chain symmetry identifiers
  # for an A3B2 stoichiometry the meaning of these features is as follows:
  # asym_id:    1 2 3 4 5
  # entity_id:  1 1 1 2 2
  # sym_id:     1 2 3 1 2
  asym_id: xnp_ndarray
  entity_id: xnp_ndarray
  sym_id: xnp_ndarray

  # token type features
  is_protein: xnp_ndarray
  is_rna: xnp_ndarray
  is_dna: xnp_ndarray
  is_ligand: xnp_ndarray
  is_nonstandard_polymer_chain: xnp_ndarray
  is_water: xnp_ndarray

  @classmethod
  def compute_features(
      cls,
      all_tokens: atom_layout.AtomLayout,
      padding_shapes: PaddingShapes,
  ) -> Self:
    """Compute the per-token features."""

    residue_index = all_tokens.res_id.astype(np.int32)

    token_index = np.arange(1, len(all_tokens.atom_name) + 1).astype(np.int32)

    aatype = []
    for res_name, chain_type in zip(all_tokens.res_name, all_tokens.chain_type):
      if chain_type in mmcif_names.POLYMER_CHAIN_TYPES:
        res_name = mmcif_names.fix_non_standard_polymer_res(
            res_name=res_name, chain_type=chain_type
        )
        if (
            chain_type == mmcif_names.DNA_CHAIN
            and res_name == residue_names.UNK_DNA
        ):
          res_name = residue_names.UNK_NUCLEIC_ONE_LETTER
      elif chain_type in mmcif_names.NON_POLYMER_CHAIN_TYPES:
        res_name = residue_names.UNK
      else:
        raise ValueError(f'Chain type {chain_type} not polymer or ligand.')
      aa = residue_names.POLYMER_TYPES_ORDER_WITH_UNKNOWN_AND_GAP[res_name]
      aatype.append(aa)
    aatype = np.array(aatype, dtype=np.int32)

    mask = np.ones(all_tokens.shape[0], dtype=bool)
    chains = _compute_asym_entity_and_sym_id(all_tokens)
    m = dict(zip(chains.chain_id, chains.asym_id))
    asym_id = np.array([m[c] for c in all_tokens.chain_id], dtype=np.int32)

    m = dict(zip(chains.chain_id, chains.entity_id))
    entity_id = np.array([m[c] for c in all_tokens.chain_id], dtype=np.int32)

    m = dict(zip(chains.chain_id, chains.sym_id))
    sym_id = np.array([m[c] for c in all_tokens.chain_id], dtype=np.int32)

    seq_length = np.array(all_tokens.shape[0], dtype=np.int32)

    is_protein = all_tokens.chain_type == mmcif_names.PROTEIN_CHAIN
    is_rna = all_tokens.chain_type == mmcif_names.RNA_CHAIN
    is_dna = all_tokens.chain_type == mmcif_names.DNA_CHAIN
    is_ligand = np.isin(
        all_tokens.chain_type, list(mmcif_names.LIGAND_CHAIN_TYPES)
    )
    standard_polymer_chain = list(mmcif_names.NON_POLYMER_CHAIN_TYPES) + list(
        mmcif_names.STANDARD_POLYMER_CHAIN_TYPES
    )
    is_nonstandard_polymer_chain = np.isin(
        all_tokens.chain_type, standard_polymer_chain, invert=True
    )
    is_water = all_tokens.chain_type == mmcif_names.WATER

    return TokenFeatures(
        residue_index=_pad_to(residue_index, (padding_shapes.num_tokens,)),
        token_index=_pad_to(token_index, (padding_shapes.num_tokens,)),
        aatype=_pad_to(aatype, (padding_shapes.num_tokens,)),
        mask=_pad_to(mask, (padding_shapes.num_tokens,)),
        asym_id=_pad_to(asym_id, (padding_shapes.num_tokens,)),
        entity_id=_pad_to(entity_id, (padding_shapes.num_tokens,)),
        sym_id=_pad_to(sym_id, (padding_shapes.num_tokens,)),
        seq_length=seq_length,
        is_protein=_pad_to(is_protein, (padding_shapes.num_tokens,)),
        is_rna=_pad_to(is_rna, (padding_shapes.num_tokens,)),
        is_dna=_pad_to(is_dna, (padding_shapes.num_tokens,)),
        is_ligand=_pad_to(is_ligand, (padding_shapes.num_tokens,)),
        is_nonstandard_polymer_chain=_pad_to(
            is_nonstandard_polymer_chain, (padding_shapes.num_tokens,)
        ),
        is_water=_pad_to(is_water, (padding_shapes.num_tokens,)),
    )

  @classmethod
  def from_data_dict(cls, batch: BatchDict) -> Self:
    return cls(
        residue_index=batch['residue_index'],
        token_index=batch['token_index'],
        aatype=batch['aatype'],
        mask=batch['seq_mask'],
        entity_id=batch['entity_id'],
        asym_id=batch['asym_id'],
        sym_id=batch['sym_id'],
        seq_length=batch['seq_length'],
        is_protein=batch['is_protein'],
        is_rna=batch['is_rna'],
        is_dna=batch['is_dna'],
        is_ligand=batch['is_ligand'],
        is_nonstandard_polymer_chain=batch['is_nonstandard_polymer_chain'],
        is_water=batch['is_water'],
    )

  def as_data_dict(self) -> BatchDict:
    return {
        'residue_index': self.residue_index,
        'token_index': self.token_index,
        'aatype': self.aatype,
        'seq_mask': self.mask,
        'entity_id': self.entity_id,
        'asym_id': self.asym_id,
        'sym_id': self.sym_id,
        'seq_length': self.seq_length,
        'is_protein': self.is_protein,
        'is_rna': self.is_rna,
        'is_dna': self.is_dna,
        'is_ligand': self.is_ligand,
        'is_nonstandard_polymer_chain': self.is_nonstandard_polymer_chain,
        'is_water': self.is_water,
    }


@chex.dataclass(mappable_dataclass=False, frozen=True)
class PredictedStructureInfo:
  """Contains information necessary to work with predicted structure."""

  atom_mask: xnp_ndarray
  residue_center_index: xnp_ndarray

  @classmethod
  def compute_features(
      cls,
      all_tokens: atom_layout.AtomLayout,
      all_token_atoms_layout: atom_layout.AtomLayout,
      padding_shapes: PaddingShapes,
  ) -> Self:
    """Compute the PredictedStructureInfo features.

    Args:
      all_tokens: flat AtomLayout with 1 representative atom per token, shape
        (num_tokens,)
      all_token_atoms_layout: AtomLayout for all atoms per token, shape
        (num_tokens, max_atoms_per_token)
      padding_shapes: padding shapes.

    Returns:
      A PredictedStructureInfo object.
    """
    atom_mask = _pad_to(
        all_token_atoms_layout.atom_name.astype(bool),
        (padding_shapes.num_tokens, None),
    )
    residue_center_index = np.zeros(padding_shapes.num_tokens, dtype=np.int32)
    for idx in range(all_tokens.shape[0]):
      repr_atom = all_tokens.atom_name[idx]
      atoms = list(all_token_atoms_layout.atom_name[idx, :])
      if repr_atom in atoms:
        residue_center_index[idx] = atoms.index(repr_atom)
      else:
        # Representative atoms can be missing if cropping the number of atoms
        # per residue.
        logging.warning(
            'The representative atom in all_tokens (%s) is not in '
            'all_token_atoms_layout (%s)',
            all_tokens[idx : idx + 1],
            all_token_atoms_layout[idx, :],
        )
        residue_center_index[idx] = 0
    return cls(atom_mask=atom_mask, residue_center_index=residue_center_index)

  @classmethod
  def from_data_dict(cls, batch: BatchDict) -> Self:
    return cls(
        atom_mask=batch['pred_dense_atom_mask'],
        residue_center_index=batch['residue_center_index'],
    )

  def as_data_dict(self) -> BatchDict:
    return {
        'pred_dense_atom_mask': self.atom_mask,
        'residue_center_index': self.residue_center_index,
    }


@chex.dataclass(mappable_dataclass=False, frozen=True)
class PolymerLigandBondInfo:
  """Contains information about polymer-ligand bonds."""

  tokens_to_polymer_ligand_bonds: atom_layout.GatherInfo
  # Gather indices to convert from cropped dense atom layout to bonds layout
  # (num_tokens, 2)
  token_atoms_to_bonds: atom_layout.GatherInfo

  @classmethod
  def compute_features(
      cls,
      all_tokens: atom_layout.AtomLayout,
      all_token_atoms_layout: atom_layout.AtomLayout,
      bond_layout: atom_layout.AtomLayout | None,
      padding_shapes: PaddingShapes,
  ) -> Self:
    """Computes the InterChainBondInfo features.

    Args:
      all_tokens: AtomLayout for tokens; shape (num_tokens,).
      all_token_atoms_layout: Atom Layout for all atoms (num_tokens,
        max_atoms_per_token)
      bond_layout: Bond layout for polymer-ligand bonds.
      padding_shapes: Padding shapes.

    Returns:
      A PolymerLigandBondInfo object.
    """

    if bond_layout is not None:
      # Must convert to list before calling np.isin, will not work raw.
      peptide_types = list(mmcif_names.PEPTIDE_CHAIN_TYPES)
      nucleic_types = list(mmcif_names.NUCLEIC_ACID_CHAIN_TYPES) + [
          mmcif_names.OTHER_CHAIN
      ]
      # These atom renames are so that we can use the atom layout code with
      # all_tokens, which only has a single atom per token.
      atom_names = bond_layout.atom_name.copy()
      atom_names[np.isin(bond_layout.chain_type, peptide_types)] = 'CA'
      atom_names[np.isin(bond_layout.chain_type, nucleic_types)] = "C1'"
      adjusted_bond_layout = atom_layout.AtomLayout(
          atom_name=atom_names,
          res_id=bond_layout.res_id,
          chain_id=bond_layout.chain_id,
          chain_type=bond_layout.chain_type,
      )
      # Remove bonds that are not in the crop.
      cropped_tokens_to_bonds = atom_layout.compute_gather_idxs(
          source_layout=all_tokens, target_layout=adjusted_bond_layout
      )
      bond_is_in_crop = np.all(
          cropped_tokens_to_bonds.gather_mask, axis=1
      ).astype(bool)
      adjusted_bond_layout = adjusted_bond_layout[bond_is_in_crop, :]
    else:
      # Create layout with correct shape when bond_layout is None.
      s = (0, 2)
      adjusted_bond_layout = atom_layout.AtomLayout(
          atom_name=np.array([], dtype=object).reshape(s),
          res_id=np.array([], dtype=int).reshape(s),
          chain_id=np.array([], dtype=object).reshape(s),
      )
    adjusted_bond_layout = adjusted_bond_layout.copy_and_pad_to(
        (padding_shapes.num_tokens, 2)
    )
    tokens_to_polymer_ligand_bonds = atom_layout.compute_gather_idxs(
        source_layout=all_tokens, target_layout=adjusted_bond_layout
    )

    # Stuff for computing the bond loss.
    if bond_layout is not None:
      # Pad to num_tokens (hoping that there are never more bonds than tokens).
      padded_bond_layout = bond_layout.copy_and_pad_to(
          (padding_shapes.num_tokens, 2)
      )
      token_atoms_to_bonds = atom_layout.compute_gather_idxs(
          source_layout=all_token_atoms_layout, target_layout=padded_bond_layout
      )
    else:
      token_atoms_to_bonds = atom_layout.GatherInfo(
          gather_idxs=np.zeros((padding_shapes.num_tokens, 2), dtype=int),
          gather_mask=np.zeros((padding_shapes.num_tokens, 2), dtype=bool),
          input_shape=np.array((
              padding_shapes.num_tokens,
              all_token_atoms_layout.shape[1],
          )),
      )

    return cls(
        tokens_to_polymer_ligand_bonds=tokens_to_polymer_ligand_bonds,
        token_atoms_to_bonds=token_atoms_to_bonds,
    )

  @classmethod
  def from_data_dict(cls, batch: BatchDict) -> Self:
    return cls(
        tokens_to_polymer_ligand_bonds=atom_layout.GatherInfo.from_dict(
            batch, key_prefix='tokens_to_polymer_ligand_bonds'
        ),
        token_atoms_to_bonds=atom_layout.GatherInfo.from_dict(
            batch, key_prefix='token_atoms_to_polymer_ligand_bonds'
        ),
    )

  def as_data_dict(self) -> BatchDict:
    return {
        **self.tokens_to_polymer_ligand_bonds.as_dict(
            key_prefix='tokens_to_polymer_ligand_bonds'
        ),
        **self.token_atoms_to_bonds.as_dict(
            key_prefix='token_atoms_to_polymer_ligand_bonds'
        ),
    }


@chex.dataclass(mappable_dataclass=False, frozen=True)
class LigandLigandBondInfo:
  """Contains information about the location of ligand-ligand bonds."""

  tokens_to_ligand_ligand_bonds: atom_layout.GatherInfo

  @classmethod
  def compute_features(
      cls,
      all_tokens: atom_layout.AtomLayout,
      bond_layout: atom_layout.AtomLayout | None,
      padding_shapes: PaddingShapes,
  ) -> Self:
    """Computes the InterChainBondInfo features.

    Args:
      all_tokens: AtomLayout for tokens; shape (num_tokens,).
      bond_layout: Bond layout for ligand-ligand bonds.
      padding_shapes: Padding shapes.

    Returns:
      A LigandLigandBondInfo object.
    """

    if bond_layout is not None:
      # Discard any bonds that do not join to an existing atom.
      keep_mask = []
      all_atom_ids = {
          uid
          for uid in zip(
              all_tokens.chain_id,
              all_tokens.res_id,
              all_tokens.atom_name,
              strict=True,
          )
      }
      for chain_id, res_id, atom_name in zip(
          bond_layout.chain_id,
          bond_layout.res_id,
          bond_layout.atom_name,
          strict=True,
      ):
        atom_a = (chain_id[0], res_id[0], atom_name[0])
        atom_b = (chain_id[1], res_id[1], atom_name[1])
        if atom_a in all_atom_ids and atom_b in all_atom_ids:
          keep_mask.append(True)
        else:
          keep_mask.append(False)
      keep_mask = np.array(keep_mask).astype(bool)
      bond_layout = bond_layout[keep_mask]
      # Remove any bonds to Hydrogen atoms.
      bond_layout = bond_layout[
          ~np.char.startswith(bond_layout.atom_name.astype(str), 'H').any(
              axis=1
          )
      ]
      atom_names = bond_layout.atom_name
      adjusted_bond_layout = atom_layout.AtomLayout(
          atom_name=atom_names,
          res_id=bond_layout.res_id,
          chain_id=bond_layout.chain_id,
          chain_type=bond_layout.chain_type,
      )
    else:
      # Create layout with correct shape when bond_layout is None.
      s = (0, 2)
      adjusted_bond_layout = atom_layout.AtomLayout(
          atom_name=np.array([], dtype=object).reshape(s),
          res_id=np.array([], dtype=int).reshape(s),
          chain_id=np.array([], dtype=object).reshape(s),
      )
    # 10 x num_tokens as max_inter_bonds_ratio + max_intra_bonds_ration = 2.061.
    adjusted_bond_layout = adjusted_bond_layout.copy_and_pad_to(
        (padding_shapes.num_tokens * 10, 2)
    )
    gather_idx = atom_layout.compute_gather_idxs(
        source_layout=all_tokens, target_layout=adjusted_bond_layout
    )
    return cls(tokens_to_ligand_ligand_bonds=gather_idx)

  @classmethod
  def from_data_dict(cls, batch: BatchDict) -> Self:
    return cls(
        tokens_to_ligand_ligand_bonds=atom_layout.GatherInfo.from_dict(
            batch, key_prefix='tokens_to_ligand_ligand_bonds'
        )
    )

  def as_data_dict(self) -> BatchDict:
    return {
        **self.tokens_to_ligand_ligand_bonds.as_dict(
            key_prefix='tokens_to_ligand_ligand_bonds'
        )
    }


@chex.dataclass(mappable_dataclass=False, frozen=True)
class PseudoBetaInfo:
  """Contains information for extracting pseudo-beta and equivalent atoms."""

  token_atoms_to_pseudo_beta: atom_layout.GatherInfo

  @classmethod
  def compute_features(
      cls,
      all_token_atoms_layout: atom_layout.AtomLayout,
      ccd: chemical_components.Ccd,
      padding_shapes: PaddingShapes,
      logging_name: str,
  ) -> Self:
    """Compute the PseudoBetaInfo features.

    Args:
      all_token_atoms_layout:  AtomLayout for all atoms per token, shape
        (num_tokens, max_atoms_per_token)
      ccd: The chemical components dictionary.
      padding_shapes: padding shapes.
      logging_name: logging name for debugging (usually the mmcif_id)

    Returns:
      A PseudoBetaInfo object.
    """
    token_idxs = []
    atom_idxs = []
    for token_idx in range(all_token_atoms_layout.shape[0]):
      chain_type = all_token_atoms_layout.chain_type[token_idx, 0]
      atom_names = list(all_token_atoms_layout.atom_name[token_idx, :])
      atom_idx = None
      is_nucleic_backbone = (
          chain_type in mmcif_names.NUCLEIC_ACID_CHAIN_TYPES
          or chain_type == mmcif_names.OTHER_CHAIN
      )
      if chain_type == mmcif_names.PROTEIN_CHAIN:
        # Protein chains
        if 'CB' in atom_names:
          atom_idx = atom_names.index('CB')
        elif 'CA' in atom_names:
          atom_idx = atom_names.index('CA')
      elif is_nucleic_backbone:
        # RNA / DNA chains
        res_name = all_token_atoms_layout.res_name[token_idx, 0]
        cifdict = ccd.get(res_name)
        if cifdict:
          parent = cifdict['_chem_comp.mon_nstd_parent_comp_id'][0]
          if parent != '?':
            res_name = parent
        if res_name in {'A', 'G', 'DA', 'DG'}:
          if 'C4' in atom_names:
            atom_idx = atom_names.index('C4')
        else:
          if 'C2' in atom_names:
            atom_idx = atom_names.index('C2')
      elif chain_type in mmcif_names.NON_POLYMER_CHAIN_TYPES:
        # Ligands: there is only one atom per token
        atom_idx = 0
      else:
        logging.warning(
            '%s: Unknown chain type for token %i. (%s)',
            logging_name,
            token_idx,
            all_token_atoms_layout[token_idx : token_idx + 1],
        )
        atom_idx = 0
      if atom_idx is None:
        (valid_atom_idxs,) = np.nonzero(
            all_token_atoms_layout.atom_name[token_idx, :]
        )
        if valid_atom_idxs.shape[0] > 0:
          atom_idx = valid_atom_idxs[0]
        else:
          atom_idx = 0
        logging.warning(
            '%s token %i (%s), does not contain a pseudo-beta atom.'
            'Using first valid atom (%s) instead.',
            logging_name,
            token_idx,
            all_token_atoms_layout[token_idx : token_idx + 1],
            all_token_atoms_layout.atom_name[token_idx, atom_idx],
        )

      token_idxs.append(token_idx)
      atom_idxs.append(atom_idx)

    pseudo_beta_layout = all_token_atoms_layout[token_idxs, atom_idxs]
    pseudo_beta_layout = pseudo_beta_layout.copy_and_pad_to((
        padding_shapes.num_tokens,
    ))
    token_atoms_to_pseudo_beta = atom_layout.compute_gather_idxs(
        source_layout=all_token_atoms_layout, target_layout=pseudo_beta_layout
    )

    return cls(
        token_atoms_to_pseudo_beta=token_atoms_to_pseudo_beta,
    )

  @classmethod
  def from_data_dict(cls, batch: BatchDict) -> Self:
    return cls(
        token_atoms_to_pseudo_beta=atom_layout.GatherInfo.from_dict(
            batch, key_prefix='token_atoms_to_pseudo_beta'
        ),
    )

  def as_data_dict(self) -> BatchDict:
    return {
        **self.token_atoms_to_pseudo_beta.as_dict(
            key_prefix='token_atoms_to_pseudo_beta'
        ),
    }


_DEFAULT_BLANK_REF = {
    'positions': np.zeros(3),
    'mask': 0,
    'element': 0,
    'charge': 0,
    'atom_name_chars': np.zeros(4),
}


def random_rotation(random_state: np.random.RandomState) -> np.ndarray:
  # Create a random rotation (Gram-Schmidt orthogonalization of two
  # random normal vectors)
  v0, v1 = random_state.normal(size=(2, 3))
  e0 = v0 / np.maximum(1e-10, np.linalg.norm(v0))
  v1 = v1 - e0 * np.dot(v1, e0)
  e1 = v1 / np.maximum(1e-10, np.linalg.norm(v1))
  e2 = np.cross(e0, e1)
  return np.stack([e0, e1, e2])


def random_augmentation(
    positions: np.ndarray,
    random_state: np.random.RandomState,
) -> np.ndarray:
  """Center then apply random translation and rotation."""

  center = np.mean(positions, axis=0)
  rot = random_rotation(random_state)
  positions_target = np.einsum('ij,kj->ki', rot, positions - center)

  translation = random_state.normal(size=(3,))
  positions_target = positions_target + translation
  return positions_target


def _get_reference_positions_from_ccd_cif(
    ccd_cif: cif_dict.CifDict,
    ref_max_modified_date: datetime.date,
    logging_name: str,
) -> np.ndarray:
  """Creates reference positions from a CCD mmcif data block."""
  num_atoms = len(ccd_cif['_chem_comp_atom.atom_id'])
  if '_chem_comp_atom.pdbx_model_Cartn_x_ideal' in ccd_cif:
    atom_x = ccd_cif['_chem_comp_atom.pdbx_model_Cartn_x_ideal']
    atom_y = ccd_cif['_chem_comp_atom.pdbx_model_Cartn_y_ideal']
    atom_z = ccd_cif['_chem_comp_atom.pdbx_model_Cartn_z_ideal']
  else:
    atom_x = np.array(['?'] * num_atoms)
    atom_y = np.array(['?'] * num_atoms)
    atom_z = np.array(['?'] * num_atoms)
  pos = np.array([[x, y, z] for x, y, z in zip(atom_x, atom_y, atom_z)])
  # Unknown reference coordinates are specified by '?' in chem comp dict.
  # Replace unknown reference coords with 0.
  if '?' in pos and '_chem_comp.pdbx_modified_date' in ccd_cif:
    # Use reference coordinates if modifed date is before cutoff.
    modified_dates = [
        datetime.date.fromisoformat(date)
        for date in ccd_cif['_chem_comp.pdbx_modified_date']
    ]
    max_modified_date = max(modified_dates)
    if max_modified_date < ref_max_modified_date:
      atom_x = ccd_cif['_chem_comp_atom.model_Cartn_x']
      atom_y = ccd_cif['_chem_comp_atom.model_Cartn_y']
      atom_z = ccd_cif['_chem_comp_atom.model_Cartn_z']
      pos = np.array([[x, y, z] for x, y, z in zip(atom_x, atom_y, atom_z)])
  if '?' in pos:
    if np.all(pos == '?'):
      logging.warning('All ref positions unknown for: %s', logging_name)
    else:
      logging.warning('Some ref positions unknown for: %s', logging_name)
    pos[pos == '?'] = 0
  return np.array(pos, dtype=np.float32)


def get_reference(
    res_name: str,
    chemical_components_data: struc_chem_comps.ChemicalComponentsData,
    ccd: chemical_components.Ccd,
    random_state: np.random.RandomState,
    ref_max_modified_date: datetime.date,
    conformer_max_iterations: int | None,
) -> tuple[dict[str, Any], Any, Any]:
  """Reference structure for residue from CCD or SMILES.

  Uses CCD entry if available, otherwise uses SMILES from chemical components
  data. Conformer generation is done using RDKit, with a fallback to CCD ideal
  or reference coordinates if RDKit fails and those coordinates are supplied.

  Args:
    res_name: ccd code of the residue.
    chemical_components_data: ChemicalComponentsData for making ref structure.
    ccd: The chemical components dictionary.
    random_state: Numpy RandomState
    ref_max_modified_date: date beyond which reference structures must not be
      modified to be allowed to use reference coordinates.
    conformer_max_iterations: Optional override for maximum number of iterations
      to run for RDKit conformer search.

  Returns:
    Mapping from atom names to features, from_atoms, dest_atoms.
  """

  ccd_cif = ccd.get(res_name)

  mol = None
  if ccd_cif:
    try:
      mol = rdkit_utils.mol_from_ccd_cif(ccd_cif, remove_hydrogens=False)
    except rdkit_utils.MolFromMmcifError:
      logging.warning('Failed to construct mol from ccd_cif for: %s', res_name)
  else:  # No CCD entry, use SMILES from chemical components data.
    if not (
        chemical_components_data.chem_comp
        and res_name in chemical_components_data.chem_comp
        and chemical_components_data.chem_comp[res_name].pdbx_smiles
    ):
      raise ValueError(f'No CCD entry or SMILES for {res_name}.')
    smiles_string = chemical_components_data.chem_comp[res_name].pdbx_smiles
    logging.info('Using SMILES for: %s - %s', res_name, smiles_string)

    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
      # In this case the model will not have any information about this molecule
      # and will not be able to predict anything about it.
      raise ValueError(
          f'Failed to construct RDKit Mol for {res_name} from SMILES string: '
          f'{smiles_string} . This is likely due to an issue with the SMILES '
          'string. Note that the userCCD input format provides an alternative '
          'way to define custom molecules directly without RDKit or SMILES.'
      )
    mol = Chem.AddHs(mol)
    # No existing names, we assign them from the graph.
    mol = rdkit_utils.assign_atom_names_from_graph(mol)
    # Temporary CCD cif with just atom and bond information, no coordinates.
    ccd_cif = rdkit_utils.mol_to_ccd_cif(mol, component_id='fake_cif')

  conformer = None
  atom_names = []
  elements = []
  charges = []
  pos = []

  # If mol is not None (must be True for SMILES case), then we try and generate
  # an RDKit conformer.
  if mol is not None:
    conformer_random_seed = int(random_state.randint(1, 1 << 31))
    conformer = rdkit_utils.get_random_conformer(
        mol=mol,
        random_seed=conformer_random_seed,
        max_iterations=conformer_max_iterations,
        logging_name=res_name,
    )
    if conformer:
      for idx, atom in enumerate(mol.GetAtoms()):
        atom_names.append(atom.GetProp('atom_name'))
        elements.append(atom.GetAtomicNum())
        charges.append(atom.GetFormalCharge())
        coords = conformer.GetAtomPosition(idx)
        pos.append([coords.x, coords.y, coords.z])
      pos = np.array(pos, dtype=np.float32)

  # If no mol could be generated (can only happen when using CCD), or no
  # conformer could be generated from the mol (can happen in either case), then
  # use CCD cif instead (which will have zero coordinates for SMILES case).
  if conformer is None:
    atom_names = ccd_cif['_chem_comp_atom.atom_id']
    charges = ccd_cif['_chem_comp_atom.charge']
    type_symbols = ccd_cif['_chem_comp_atom.type_symbol']
    elements = [
        periodic_table.ATOMIC_NUMBER.get(elem_type.capitalize(), 0)
        for elem_type in type_symbols
    ]
    pos = _get_reference_positions_from_ccd_cif(
        ccd_cif=ccd_cif,
        ref_max_modified_date=ref_max_modified_date,
        logging_name=res_name,
    )

  # Augment reference positions.
  pos = random_augmentation(pos, random_state)

  # Extract atom and bond information from CCD cif.
  from_atom = ccd_cif.get('_chem_comp_bond.atom_id_1', None)
  dest_atom = ccd_cif.get('_chem_comp_bond.atom_id_2', None)

  features = {}
  for atom_name in atom_names:
    features[atom_name] = {}
    idx = atom_names.index(atom_name)
    charge = 0 if charges[idx] == '?' else int(charges[idx])
    atom_name_chars = np.array([ord(c) - 32 for c in atom_name], dtype=int)
    atom_name_chars = _pad_to(atom_name_chars, (4,))
    features[atom_name]['positions'] = pos[idx]
    features[atom_name]['mask'] = 1
    features[atom_name]['element'] = elements[idx]
    features[atom_name]['charge'] = charge
    features[atom_name]['atom_name_chars'] = atom_name_chars
  return features, from_atom, dest_atom


@chex.dataclass(mappable_dataclass=False, frozen=True)
class RefStructure:
  """Contains ref structure information."""

  # Array with positions, float32, shape [num_res, max_atoms_per_token, 3]
  positions: xnp_ndarray
  # Array with masks, bool, shape [num_res, max_atoms_per_token]
  mask: xnp_ndarray
  # Array with elements, int32, shape [num_res, max_atoms_per_token]
  element: xnp_ndarray
  # Array with charges, float32, shape [num_res, max_atoms_per_token]
  charge: xnp_ndarray
  # Array with atom name characters, int32, [num_res, max_atoms_per_token, 4]
  atom_name_chars: xnp_ndarray
  # Array with reference space uids, int32, [num_res, max_atoms_per_token]
  ref_space_uid: xnp_ndarray

  @classmethod
  def compute_features(
      cls,
      all_token_atoms_layout: atom_layout.AtomLayout,
      ccd: chemical_components.Ccd,
      padding_shapes: PaddingShapes,
      chemical_components_data: struc_chem_comps.ChemicalComponentsData,
      random_state: np.random.RandomState,
      ref_max_modified_date: datetime.date,
      conformer_max_iterations: int | None,
      ligand_ligand_bonds: atom_layout.AtomLayout | None = None,
  ) -> tuple[Self, Any]:
    """Reference structure information for each residue."""

    # Get features per atom
    padded_shape = (padding_shapes.num_tokens, all_token_atoms_layout.shape[1])
    result = {
        'positions': np.zeros((*padded_shape, 3), 'float32'),
        'mask': np.zeros(padded_shape, 'bool'),
        'element': np.zeros(padded_shape, 'int32'),
        'charge': np.zeros(padded_shape, 'float32'),
        'atom_name_chars': np.zeros((*padded_shape, 4), 'int32'),
        'ref_space_uid': np.zeros((*padded_shape,), 'int32'),
    }

    atom_names_all = []
    chain_ids_all = []
    res_ids_all = []

    # Cache reference conformations for each residue.
    conformations = {}
    ref_space_uids = {}
    for idx in np.ndindex(all_token_atoms_layout.shape):
      chain_id = all_token_atoms_layout.chain_id[idx]
      res_id = all_token_atoms_layout.res_id[idx]
      res_name = all_token_atoms_layout.res_name[idx]
      is_non_standard = res_name not in _STANDARD_RESIDUES
      atom_name = all_token_atoms_layout.atom_name[idx]
      if not atom_name:
        ref = _DEFAULT_BLANK_REF
      else:
        if (chain_id, res_id) not in conformations:
          conf, from_atom, dest_atom = get_reference(
              res_name=res_name,
              chemical_components_data=chemical_components_data,
              ccd=ccd,
              random_state=random_state,
              ref_max_modified_date=ref_max_modified_date,
              conformer_max_iterations=conformer_max_iterations,
          )
          conformations[(chain_id, res_id)] = conf

          if (
              is_non_standard
              and (from_atom is not None)
              and (dest_atom is not None)
          ):
            # Add intra-ligand bond graph
            atom_names_ligand = np.stack(
                [from_atom, dest_atom], axis=1, dtype=object
            )
            atom_names_all.append(atom_names_ligand)
            res_ids_all.append(
                np.full_like(atom_names_ligand, res_id, dtype=int)
            )
            chain_ids_all.append(
                np.full_like(atom_names_ligand, chain_id, dtype=object)
            )

        conformation = conformations.get(
            (chain_id, res_id), {atom_name: _DEFAULT_BLANK_REF}
        )
        if atom_name not in conformation:
          logging.warning(
              'Missing atom "%s" for CCD "%s"',
              atom_name,
              all_token_atoms_layout.res_name[idx],
          )
        ref = conformation.get(atom_name, _DEFAULT_BLANK_REF)
      for k in ref:
        result[k][idx] = ref[k]

      # Assign a unique reference space id to each component, to determine which
      # reference positions live in the same reference space.
      space_str_id = (
          all_token_atoms_layout.chain_id[idx],
          all_token_atoms_layout.res_id[idx],
      )
      if space_str_id not in ref_space_uids:
        ref_space_uids[space_str_id] = len(ref_space_uids)
      result['ref_space_uid'][idx] = ref_space_uids[space_str_id]

    if atom_names_all:
      atom_names_all = np.concatenate(atom_names_all, axis=0)
      res_ids_all = np.concatenate(res_ids_all, axis=0)
      chain_ids_all = np.concatenate(chain_ids_all, axis=0)
      if ligand_ligand_bonds is not None:
        adjusted_ligand_ligand_bonds = atom_layout.AtomLayout(
            atom_name=np.concatenate(
                [ligand_ligand_bonds.atom_name, atom_names_all], axis=0
            ),
            chain_id=np.concatenate(
                [ligand_ligand_bonds.chain_id, chain_ids_all], axis=0
            ),
            res_id=np.concatenate(
                [ligand_ligand_bonds.res_id, res_ids_all], axis=0
            ),
        )
      else:
        adjusted_ligand_ligand_bonds = atom_layout.AtomLayout(
            atom_name=atom_names_all,
            chain_id=chain_ids_all,
            res_id=res_ids_all,
        )
    else:
      adjusted_ligand_ligand_bonds = ligand_ligand_bonds

    return cls(**result), adjusted_ligand_ligand_bonds

  @classmethod
  def from_data_dict(cls, batch: BatchDict) -> Self:
    return cls(
        positions=batch['ref_pos'],
        mask=batch['ref_mask'],
        element=batch['ref_element'],
        charge=batch['ref_charge'],
        atom_name_chars=batch['ref_atom_name_chars'],
        ref_space_uid=batch['ref_space_uid'],
    )

  def as_data_dict(self) -> BatchDict:
    return {
        'ref_pos': self.positions,
        'ref_mask': self.mask,
        'ref_element': self.element,
        'ref_charge': self.charge,
        'ref_atom_name_chars': self.atom_name_chars,
        'ref_space_uid': self.ref_space_uid,
    }


@chex.dataclass(mappable_dataclass=False, frozen=True)
class ConvertModelOutput:
  """Contains atom layout info."""

  cleaned_struc: structure.Structure
  token_atoms_layout: atom_layout.AtomLayout
  flat_output_layout: atom_layout.AtomLayout
  empty_output_struc: structure.Structure
  polymer_ligand_bonds: atom_layout.AtomLayout
  ligand_ligand_bonds: atom_layout.AtomLayout

  @classmethod
  def compute_features(
      cls,
      all_token_atoms_layout: atom_layout.AtomLayout,
      padding_shapes: PaddingShapes,
      cleaned_struc: structure.Structure,
      flat_output_layout: atom_layout.AtomLayout,
      empty_output_struc: structure.Structure,
      polymer_ligand_bonds: atom_layout.AtomLayout,
      ligand_ligand_bonds: atom_layout.AtomLayout,
  ) -> Self:
    """Pads the all_token_atoms_layout and stores other data."""
    # Crop and pad the all_token_atoms_layout.
    token_atoms_layout = all_token_atoms_layout.copy_and_pad_to(
        (padding_shapes.num_tokens, all_token_atoms_layout.shape[1])
    )

    return cls(
        cleaned_struc=cleaned_struc,
        token_atoms_layout=token_atoms_layout,
        flat_output_layout=flat_output_layout,
        empty_output_struc=empty_output_struc,
        polymer_ligand_bonds=polymer_ligand_bonds,
        ligand_ligand_bonds=ligand_ligand_bonds,
    )

  @classmethod
  def from_data_dict(cls, batch: BatchDict) -> Self:
    """Construct atom layout object from dictionary."""

    return cls(
        cleaned_struc=_unwrap(batch.get('cleaned_struc', None)),
        token_atoms_layout=_unwrap(batch.get('token_atoms_layout', None)),
        flat_output_layout=_unwrap(batch.get('flat_output_layout', None)),
        empty_output_struc=_unwrap(batch.get('empty_output_struc', None)),
        polymer_ligand_bonds=_unwrap(batch.get('polymer_ligand_bonds', None)),
        ligand_ligand_bonds=_unwrap(batch.get('ligand_ligand_bonds', None)),
    )

  def as_data_dict(self) -> BatchDict:
    return {
        'cleaned_struc': np.array(self.cleaned_struc, object),
        'token_atoms_layout': np.array(self.token_atoms_layout, object),
        'flat_output_layout': np.array(self.flat_output_layout, object),
        'empty_output_struc': np.array(self.empty_output_struc, object),
        'polymer_ligand_bonds': np.array(self.polymer_ligand_bonds, object),
        'ligand_ligand_bonds': np.array(self.ligand_ligand_bonds, object),
    }


@chex.dataclass(mappable_dataclass=False, frozen=True)
class AtomCrossAtt:
  """Operate on flat atoms."""

  token_atoms_to_queries: atom_layout.GatherInfo
  tokens_to_queries: atom_layout.GatherInfo
  tokens_to_keys: atom_layout.GatherInfo
  queries_to_keys: atom_layout.GatherInfo
  queries_to_token_atoms: atom_layout.GatherInfo

  @classmethod
  def compute_features(
      cls,
      all_token_atoms_layout: atom_layout.AtomLayout,  # (num_tokens, num_dense)
      queries_subset_size: int,
      keys_subset_size: int,
      padding_shapes: PaddingShapes,
  ) -> Self:
    """Computes gather indices and meta data to work with a flat atom list."""

    token_atoms_layout = all_token_atoms_layout.copy_and_pad_to(
        (padding_shapes.num_tokens, all_token_atoms_layout.shape[1])
    )
    token_atoms_mask = token_atoms_layout.atom_name.astype(bool)
    flat_layout = token_atoms_layout[token_atoms_mask]
    num_atoms = flat_layout.shape[0]

    padded_flat_layout = flat_layout.copy_and_pad_to((
        padding_shapes.num_atoms,
    ))

    # Create the layout for queries
    num_subsets = padding_shapes.num_atoms // queries_subset_size
    lay_arr = padded_flat_layout.to_array()
    queries_layout = atom_layout.AtomLayout.from_array(
        lay_arr.reshape((6, num_subsets, queries_subset_size))
    )

    # Create the layout for the keys (the key subsets are centered around the
    # query subsets)
    # Create initial gather indices (contain out-of-bound indices)
    subset_centers = np.arange(
        queries_subset_size / 2, padding_shapes.num_atoms, queries_subset_size
    )
    flat_to_key_gathers = (
        subset_centers[:, None]
        + np.arange(-keys_subset_size / 2, keys_subset_size / 2)[None, :]
    )
    flat_to_key_gathers = flat_to_key_gathers.astype(int)
    # Shift subsets with out-of-bound indices, such that they are fully within
    # the bounds.
    for row in range(flat_to_key_gathers.shape[0]):
      if flat_to_key_gathers[row, 0] < 0:
        flat_to_key_gathers[row, :] -= flat_to_key_gathers[row, 0]
      elif flat_to_key_gathers[row, -1] > num_atoms - 1:
        overflow = flat_to_key_gathers[row, -1] - (num_atoms - 1)
        flat_to_key_gathers[row, :] -= overflow
    # Create the keys layout.
    keys_layout = padded_flat_layout[flat_to_key_gathers]

    # Create gather indices for conversion between token atoms layout,
    # queries layout and keys layout.
    token_atoms_to_queries = atom_layout.compute_gather_idxs(
        source_layout=token_atoms_layout, target_layout=queries_layout
    )

    token_atoms_to_keys = atom_layout.compute_gather_idxs(
        source_layout=token_atoms_layout, target_layout=keys_layout
    )

    queries_to_keys = atom_layout.compute_gather_idxs(
        source_layout=queries_layout, target_layout=keys_layout
    )

    queries_to_token_atoms = atom_layout.compute_gather_idxs(
        source_layout=queries_layout, target_layout=token_atoms_layout
    )

    # Create gather indices for conversion of tokens layout to
    # queries and keys layout
    token_idxs = np.arange(padding_shapes.num_tokens).astype(np.int64)
    token_idxs = np.broadcast_to(token_idxs[:, None], token_atoms_layout.shape)
    tokens_to_queries = atom_layout.GatherInfo(
        gather_idxs=atom_layout.convert(
            token_atoms_to_queries, token_idxs, layout_axes=(0, 1)
        ),
        gather_mask=atom_layout.convert(
            token_atoms_to_queries, token_atoms_mask, layout_axes=(0, 1)
        ),
        input_shape=np.array((padding_shapes.num_tokens,)),
    )

    tokens_to_keys = atom_layout.GatherInfo(
        gather_idxs=atom_layout.convert(
            token_atoms_to_keys, token_idxs, layout_axes=(0, 1)
        ),
        gather_mask=atom_layout.convert(
            token_atoms_to_keys, token_atoms_mask, layout_axes=(0, 1)
        ),
        input_shape=np.array((padding_shapes.num_tokens,)),
    )

    return cls(
        token_atoms_to_queries=token_atoms_to_queries,
        tokens_to_queries=tokens_to_queries,
        tokens_to_keys=tokens_to_keys,
        queries_to_keys=queries_to_keys,
        queries_to_token_atoms=queries_to_token_atoms,
    )

  @classmethod
  def from_data_dict(cls, batch: BatchDict) -> Self:
    return cls(
        token_atoms_to_queries=atom_layout.GatherInfo.from_dict(
            batch, key_prefix='token_atoms_to_queries'
        ),
        tokens_to_queries=atom_layout.GatherInfo.from_dict(
            batch, key_prefix='tokens_to_queries'
        ),
        tokens_to_keys=atom_layout.GatherInfo.from_dict(
            batch, key_prefix='tokens_to_keys'
        ),
        queries_to_keys=atom_layout.GatherInfo.from_dict(
            batch, key_prefix='queries_to_keys'
        ),
        queries_to_token_atoms=atom_layout.GatherInfo.from_dict(
            batch, key_prefix='queries_to_token_atoms'
        ),
    )

  def as_data_dict(self) -> BatchDict:
    return {
        **self.token_atoms_to_queries.as_dict(
            key_prefix='token_atoms_to_queries'
        ),
        **self.tokens_to_queries.as_dict(key_prefix='tokens_to_queries'),
        **self.tokens_to_keys.as_dict(key_prefix='tokens_to_keys'),
        **self.queries_to_keys.as_dict(key_prefix='queries_to_keys'),
        **self.queries_to_token_atoms.as_dict(
            key_prefix='queries_to_token_atoms'
        ),
    }


@chex.dataclass(mappable_dataclass=False, frozen=True)
class Frames:
  """Features for backbone frames."""

  mask: xnp_ndarray

  @classmethod
  def compute_features(
      cls,
      all_tokens: atom_layout.AtomLayout,
      all_token_atoms_layout: atom_layout.AtomLayout,
      ref_structure: RefStructure,
      padding_shapes: PaddingShapes,
  ) -> Self:
    """Computes features for backbone frames."""
    num_tokens = padding_shapes.num_tokens
    all_token_atoms_layout = all_token_atoms_layout.copy_and_pad_to(
        (num_tokens, all_token_atoms_layout.shape[1])
    )

    all_token_atoms_to_all_tokens = atom_layout.compute_gather_idxs(
        source_layout=all_token_atoms_layout, target_layout=all_tokens
    )
    ref_coordinates = atom_layout.convert(
        all_token_atoms_to_all_tokens,
        ref_structure.positions.astype(np.float32),
        layout_axes=(0, 1),
    )
    ref_mask = atom_layout.convert(
        all_token_atoms_to_all_tokens,
        ref_structure.mask.astype(bool),
        layout_axes=(0, 1),
    )
    ref_mask = ref_mask & all_token_atoms_to_all_tokens.gather_mask.astype(bool)

    all_frame_mask = []

    # Iterate over tokens
    for idx, args in enumerate(
        zip(all_tokens.chain_type, all_tokens.chain_id, all_tokens.res_id)
    ):

      chain_type, chain_id, res_id = args

      if chain_type in list(mmcif_names.PEPTIDE_CHAIN_TYPES):
        frame_mask = True
      elif chain_type in list(mmcif_names.NUCLEIC_ACID_CHAIN_TYPES):
        frame_mask = True
      elif chain_type in list(mmcif_names.NON_POLYMER_CHAIN_TYPES):
        # For ligands, build frames from closest atoms from the same molecule.
        (local_token_idxs,) = np.where(
            (all_tokens.chain_type == chain_type)
            & (all_tokens.chain_id == chain_id)
            & (all_tokens.res_id == res_id)
        )

        if len(local_token_idxs) < 3:
          frame_mask = False

        else:
          # [local_tokens]
          local_dist = np.linalg.norm(
              ref_coordinates[idx] - ref_coordinates[local_token_idxs], axis=-1
          )
          local_mask = ref_mask[local_token_idxs]
          cost = local_dist + 1e8 * ~local_mask
          cost = cost + 1e8 * (idx == local_token_idxs)
          # [local_tokens]
          closest_idxs = np.argsort(cost, axis=0)

          # The closest indices index an array of local tokens. Convert this
          # to indices of the full (num_tokens,) array.
          global_closest_idxs = local_token_idxs[closest_idxs]

          # Construct frame by placing the current token at the origin and two
          # nearest atoms on either side.
          global_frame_idxs = np.array(
              (global_closest_idxs[0], idx, global_closest_idxs[1])
          )

          # Check that the frame atoms are not colinear.
          a, b, c = ref_coordinates[global_frame_idxs]
          vec1 = a - b
          vec2 = c - b
          # Reference coordinates can be all zeros, in which case we have
          # to explicitly set colinearity.
          if np.isclose(np.linalg.norm(vec1, axis=-1), 0) or np.isclose(
              np.linalg.norm(vec2, axis=-1), 0
          ):
            is_colinear = True
            logging.info('Found identical coordinates: Assigning as colinear.')
          else:
            vec1 = vec1 / np.linalg.norm(vec1, axis=-1)
            vec2 = vec2 / np.linalg.norm(vec2, axis=-1)
            cos_angle = np.einsum('...k,...k->...', vec1, vec2)
            # <25 degree deviation is considered colinear.
            is_colinear = 1 - np.abs(cos_angle) < 0.0937

          frame_mask = not is_colinear
      else:
        # No frame for other chain types.
        frame_mask = False

      all_frame_mask.append(frame_mask)

    all_frame_mask = np.array(all_frame_mask, dtype=bool)

    mask = _pad_to(all_frame_mask, (padding_shapes.num_tokens,))

    return cls(mask=mask)

  @classmethod
  def from_data_dict(cls, batch: BatchDict) -> Self:
    return cls(mask=batch['frames_mask'])

  def as_data_dict(self) -> BatchDict:
    return {'frames_mask': self.mask}
