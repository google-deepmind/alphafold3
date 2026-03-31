# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Functions for producing "paired" and "unpaired" MSA features for each chain.

The paired MSA:
- Is made from the result of the all_seqs MSA query.
- Is ordered such that you can concatenate features across chains and related
  sequences will end up on the same row. Related here means "from the same
  species". Gaps are added to facilitate this whenever a sequence has no
  suitable pair.

The unpaired MSA:
- Is made from the results of the remaining MSA queries.
- Has no special ordering properties.
- Is deduplicated such that it doesn't contain any sequences in the paired MSA.
"""

from collections.abc import Mapping, MutableMapping, Sequence
from alphafold3.model import data_constants
import numpy as np


def _align_species(
    all_species: Sequence[bytes],
    chains_species_to_rows: Sequence[Mapping[bytes, np.ndarray]],
    min_hits_per_species: Mapping[bytes, int],
) -> np.ndarray:
  """Aligns MSA row indices based on species.

  Within a species, MSAs are aligned based on their original order (the first
  sequence for a species in the first chain's MSA is aligned to the first
  sequence for the same species in the second chain's MSA).

  Args:
    all_species: A list of all unique species identifiers.
    chains_species_to_rows: A dictionary for each chain, that maps species to
      the set of MSA row indices from that species in that chain.
    min_hits_per_species: A mapping from species id, to the minimum MSA size
      across chains for that species (ignoring chains with zero hits).

  Returns:
    A matrix of size [num_msa_rows, num_chains], where the i,j element is an
    index into the jth chains MSA.  Each row consists of sequences from each
    chain for the same species (or -1 if that chain has no sequences for that
    species).
  """
  # Each species block is of size [num_seqs x num_chains] and consists of
  # indices into the respective MSAs that have been aligned and are all for the
  # same species.
  species_blocks = []
  for species in all_species:
    chain_row_indices = []
    for species_to_rows in chains_species_to_rows:
      min_msa_size = min_hits_per_species[species]
      if species not in species_to_rows:
        # If a given chain has no hits for a species then we pad it with -1's,
        # later on these values are used to make sure each feature is padded
        # with its appropriate pad value.
        row_indices = np.full(min_msa_size, fill_value=-1, dtype=np.int32)
      else:
        # We crop down to the smallest MSA for a given species across chains.
        row_indices = species_to_rows[species][:min_msa_size]
      chain_row_indices.append(row_indices)
    species_block = np.stack(chain_row_indices, axis=1)
    species_blocks.append(species_block)
  aligned_matrix = np.concatenate(species_blocks, axis=0)
  return aligned_matrix


def create_paired_features(
    chains: Sequence[MutableMapping[str, np.ndarray]],
    max_paired_sequences: int,
    nonempty_chain_ids: set[str],
    max_hits_per_species: int,
) -> Sequence[MutableMapping[str, np.ndarray]]:
  """Creates per-chain MSA features where the MSAs have been aligned.

  Args:
    chains: A list of feature dicts, one for each chain.
    max_paired_sequences: No more than this many paired sequences will be
      returned from this function.
    nonempty_chain_ids: A set of chain ids (str) that are included in the crop
      there is no reason to process chains not in this list.
    max_hits_per_species: No more than this number of sequences will be returned
      for a given species.

  Returns:
    An updated feature dictionary for each chain, where the {}_all_seq features
    have been aligned so that the nth row in chain 1 is aligned to the nth row
    in chain 2's features.
  """
  # The number of chains that the given species appears in - we rank hits
  # across more chains higher.
  species_num_chains = {}

  # For each chain we keep a mapping from species to the row indices in the
  # original MSA for that chain.
  chains_species_to_rows = []

  # Keep track of the minimum number of hits across chains for a given species.
  min_hits_per_species = {}

  for chain in chains:
    species_ids = chain['msa_species_identifiers_all_seq']

    # The query gets an empty species_id, so no pairing happens for this row.
    if (
        species_ids.size == 0
        or (species_ids.size == 1 and not species_ids[0])
        or chain['chain_id'] not in nonempty_chain_ids
    ):
      chains_species_to_rows.append({})
      continue

    # For each species keep track of which row indices in the original MSA are
    # from this species.
    row_indices = np.arange(len(species_ids))
    # The grouping np.split code requires that the input is already clustered
    # by species id.
    sort_idxs = species_ids.argsort()
    species_ids = species_ids[sort_idxs]
    row_indices = row_indices[sort_idxs]

    species, unique_row_indices = np.unique(species_ids, return_index=True)
    grouped_row_indices = np.split(row_indices, unique_row_indices[1:])
    species_to_rows = dict(zip(species, grouped_row_indices, strict=True))
    chains_species_to_rows.append(species_to_rows)

    for s in species:
      species_num_chains[s] = species_num_chains.get(s, 0) + 1

    for species, row_indices in species_to_rows.items():
      min_hits_per_species[species] = min(
          min_hits_per_species.get(species, max_hits_per_species),
          len(row_indices),
      )

  # Construct a mapping from the number of chains a species appears in to
  # the list of species with that count.
  num_chains_to_species = {}
  for species, num_chains in species_num_chains.items():
    if not species or num_chains <= 1:
      continue
    if num_chains not in num_chains_to_species:
      num_chains_to_species[num_chains] = []
    num_chains_to_species[num_chains].append(species)

  num_rows_seen = 0
  # We always keep the first row as it is the query sequence.
  all_rows = [np.array([[0] * len(chains)], dtype=np.int32)]

  # We prioritize species that have hits across more chains.
  for num_chains in sorted(num_chains_to_species, reverse=True):
    all_species = num_chains_to_species[num_chains]

    # Align all the per-chain row indices by species, so every paired row is
    # for a single species.
    rows = _align_species(
        all_species, chains_species_to_rows, min_hits_per_species
    )
    # Sort rows by the product of the original indices in the respective chain
    # MSAS, so as to rank hits that appear higher in the original MSAs higher.
    rank_metric = np.abs(np.prod(rows.astype(np.float32), axis=1))
    sorted_rows = rows[np.argsort(rank_metric), :]
    all_rows.append(sorted_rows)
    num_rows_seen += rows.shape[0]
    if num_rows_seen >= max_paired_sequences:
      break

  all_rows = np.concatenate(all_rows, axis=0)
  all_rows = all_rows[:max_paired_sequences, :]

  # Now we just have to select the relevant rows from the original msa and
  # deletion matrix features
  paired_chains = []
  for chain_idx, chain in enumerate(chains):
    out_chain = {k: v for k, v in chain.items() if 'all_seq' not in k}
    selected_row_indices = all_rows[:, chain_idx]
    for feat_name in {'msa', 'deletion_matrix'}:
      all_seq_name = f'{feat_name}_all_seq'
      feat_value = chain[all_seq_name]

      # The selected row indices are padded to be the same shape for each chain,
      # they are padded with -1's, so we add a single row onto the feature with
      # the appropriate pad value.  This has the effect that we correctly pad
      # each feature since all padded indices will select this padding row.
      pad_value = data_constants.MSA_PAD_VALUES[feat_name]
      feat_value = np.concatenate([
          feat_value,
          np.full((1, feat_value.shape[1]), pad_value, feat_value.dtype),
      ])

      feat_value = feat_value[selected_row_indices, :]
      out_chain[all_seq_name] = feat_value
    out_chain['num_alignments_all_seq'] = np.array(
        out_chain['msa_all_seq'].shape[0]
    )
    paired_chains.append(out_chain)
  return paired_chains


def deduplicate_unpaired_sequences(
    np_chains: Sequence[MutableMapping[str, np.ndarray]],
) -> Sequence[MutableMapping[str, np.ndarray]]:
  """Deduplicates unpaired sequences based on paired sequences."""

  feature_names = np_chains[0].keys()
  msa_features = (
      data_constants.NUM_SEQ_MSA_FEATURES
      + data_constants.NUM_SEQ_NUM_RES_MSA_FEATURES
  )

  for chain in np_chains:
    sequence_set = set(
        hash(s.data.tobytes()) for s in chain['msa_all_seq'].astype(np.int8)
    )
    keep_rows = []
    # Go through unpaired MSA seqs and remove any rows that correspond to the
    # sequences that are already present in the paired MSA.
    for row_num, seq in enumerate(chain['msa'].astype(np.int8)):
      if hash(seq.data.tobytes()) not in sequence_set:
        keep_rows.append(row_num)
    for feature_name in feature_names:
      if feature_name in msa_features:
        chain[feature_name] = chain[feature_name][keep_rows]
    chain['num_alignments'] = np.array(chain['msa'].shape[0], dtype=np.int32)
  return np_chains


def choose_paired_unpaired_msa_crop_sizes(
    unpaired_msa: np.ndarray,
    paired_msa: np.ndarray | None,
    total_msa_crop_size: int,
    max_paired_sequences: int,
) -> tuple[int, int | None]:
  """Returns the sizes of the MSA crop and MSA_all_seq crop.

  NOTE: Unpaired + paired MSA sizes can exceed total_msa_size when
    there are lots of gapped rows. Through the pairing logic another chain(s)
    will have fewer than total_msa_size.

  Args:
     unpaired_msa: The unpaired MSA array (not all_seq).
     paired_msa: The paired MSA array (all_seq).
     total_msa_crop_size: The maximum total number of sequences to crop to.
     max_paired_sequences: The maximum number of sequences that can come from
       MSA pairing.

  Returns:
    A tuple of:
      The size of the reduced MSA crop (not all_seq features).
      The size of the unreduced MSA crop (for all_seq features) or None, if
        paired_msa is None.
  """
  if paired_msa is not None:
    paired_crop_size = np.minimum(paired_msa.shape[0], max_paired_sequences)

    # We reduce the number of un-paired sequences, by the number of times a
    # sequence from this chains MSA is included in the paired MSA.  This keeps
    # the MSA size for each chain roughly constant.
    cropped_all_seq_msa = paired_msa[:max_paired_sequences]
    num_non_gapped_pairs = cropped_all_seq_msa.shape[0]

    assert num_non_gapped_pairs <= max_paired_sequences
    unpaired_crop_size = np.minimum(
        unpaired_msa.shape[0], total_msa_crop_size - num_non_gapped_pairs
    )
    assert unpaired_crop_size >= 0
  else:
    unpaired_crop_size = np.minimum(unpaired_msa.shape[0], total_msa_crop_size)
    paired_crop_size = None
  return unpaired_crop_size, paired_crop_size


def remove_all_gapped_rows_from_all_seqs(
    chains_list: Sequence[dict[str, np.ndarray]], asym_ids: Sequence[float]
) -> Sequence[dict[str, np.ndarray]]:
  """Removes all gapped rows from all_seq feat based on selected asym_ids."""

  merged_msa_all_seq = np.concatenate(
      [
          chain['msa_all_seq']
          for chain in chains_list
          if chain['asym_id'][0] in asym_ids
      ],
      axis=1,
  )

  non_gapped_keep_rows = np.any(
      merged_msa_all_seq != data_constants.MSA_GAP_IDX, axis=1
  )
  for chain in chains_list:
    for feat_name in list(chains_list)[0]:
      if '_all_seq' in feat_name:
        feat_name_split = feat_name.split('_all_seq')[0]
        if feat_name_split in (
            data_constants.NUM_SEQ_NUM_RES_MSA_FEATURES
            + data_constants.NUM_SEQ_MSA_FEATURES
        ):
          # For consistency we do this for all chains even though the
          # gapped rows are based on a selected set asym_ids.
          chain[feat_name] = chain[feat_name][non_gapped_keep_rows]
    chain['num_alignments_all_seq'] = np.sum(non_gapped_keep_rows)
  return chains_list
