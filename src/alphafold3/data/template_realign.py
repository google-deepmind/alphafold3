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

"""Realign sequences found in PDB seqres to the actual CIF sequences."""

from collections.abc import Mapping


class AlignmentError(Exception):
  """Failed alignment between the hit sequence and the actual mmCIF sequence."""


def realign_hit_to_structure(
    *,
    hit_sequence: str,
    hit_start_index: int,
    hit_end_index: int,
    full_length: int,
    structure_sequence: str,
    query_to_hit_mapping: Mapping[int, int],
) -> Mapping[int, int]:
  """Realigns the hit sequence to the Structure sequence.

  For example, for the given input:
    query_sequence : ABCDEFGHIJKL
    hit_sequence   : ---DEFGHIJK-
    struc_sequence : XDEFGHKL
  the mapping is {3: 0, 4: 1, 5: 2, 6: 3, 7: 4, 8: 5, 9: 6, 10: 7}. However, the
  actual Structure sequence has an extra X at the start as well as no IJ. So the
  alignment from the query to the Structure sequence will be:
    hit_sequence   : ---DEFGHIJK-
    struc_aligned  : --XDEFGH--KL
  and the new mapping will therefore be: {3: 1, 4: 2, 5: 3, 6: 4, 7: 5, 10: 6}.

  Args:
    hit_sequence: The PDB seqres hit sequence obtained from Hmmsearch, but
      without any gaps. This is not the full PDB seqres template sequence but
      rather just its subsequence from hit_start_index to hit_end_index.
    hit_start_index: The start index of the hit sequence in the full PDB seqres
      template sequence (inclusive).
    hit_end_index: The end index of the hit sequence in the full PDB seqres
      template sequence (exclusive).
    full_length: The length of the full PDB seqres template sequence.
    structure_sequence: The actual sequence extracted from the Structure
      corresponding to this template. In vast majority of cases this is the same
      as the PDB seqres sequence, but this function handles the cases when not.
    query_to_hit_mapping: The mapping from the query sequence to the
      hit_sequence.

  Raises:
    AlignmentError: if the alignment between the sequence returned by Hmmsearch
      differs from the actual sequence found in the mmCIF and can't be aligned
      using the simple alignment algorithm.

  Returns:
    A mapping from the query sequence to the actual Structure sequence.
  """
  max_num_gaps = full_length - len(structure_sequence)
  if max_num_gaps < 0:
    raise AlignmentError(
        f'The Structure sequence ({len(structure_sequence)}) '
        f'must be shorter than the PDB seqres sequence ({full_length}):\n'
        f'Structure sequence : {structure_sequence}\n'
        f'PDB seqres sequence: {hit_sequence}'
    )

  if len(hit_sequence) != hit_end_index - hit_start_index:
    raise AlignmentError(
        f'The difference of {hit_end_index=} and {hit_start_index=} does not '
        f'equal to the length of the {hit_sequence}: {len(hit_sequence)}'
    )

  best_score = -1
  best_start = 0
  best_query_to_hit_mapping = query_to_hit_mapping
  max_num_gaps_before_subseq = min(hit_start_index, max_num_gaps)
  # It is possible the gaps needed to align the PDB seqres subsequence and
  # the Structure subsequence need to be inserted before the match region.
  # Try and pick the alignment with the best number of aligned residues.
  for num_gaps_before_subseq in range(0, max_num_gaps_before_subseq + 1):
    start = hit_start_index - num_gaps_before_subseq
    end = hit_end_index - num_gaps_before_subseq
    structure_subseq = structure_sequence[start:end]

    new_query_to_hit_mapping, score = _remap_to_struc_seq(
        hit_seq=hit_sequence,
        struc_seq=structure_subseq,
        max_num_gaps=max_num_gaps - num_gaps_before_subseq,
        mapping=query_to_hit_mapping,
    )
    if score >= best_score:
      # Use >= to prefer matches with larger number of gaps before.
      best_score = score
      best_start = start
      best_query_to_hit_mapping = new_query_to_hit_mapping

  return {q: h + best_start for q, h in best_query_to_hit_mapping.items()}


def _remap_to_struc_seq(
    *,
    hit_seq: str,
    struc_seq: str,
    max_num_gaps: int,
    mapping: Mapping[int, int],
) -> tuple[Mapping[int, int], int]:
  """Remaps the query -> hit mapping to match the actual Structure sequence.

  Args:
    hit_seq: The hit sequence - a subsequence of the PDB seqres sequence without
      any Hmmsearch modifications like inserted gaps or lowercased residues.
    struc_seq: The actual sequence obtained from the corresponding Structure.
    max_num_gaps: The maximum number of gaps that can be inserted in the
      Structure sequence. In practice, this is the length difference between the
      PDB seqres sequence and the actual Structure sequence.
    mapping: The mapping from the query residues to the hit residues. This will
      be remapped to point to the actual Structure sequence using a simple
      realignment algorithm.

  Returns:
    A tuple of (mapping, score):
      * Mapping from the query to the actual Structure sequence.
      * Score which is the number of matching aligned residues.

  Raises:
    ValueError if the structure sequence isn't shorter than the seqres sequence.
    ValueError if the alignment fails.
  """
  hit_seq_idx = 0
  struc_seq_idx = 0
  hit_to_struc_seq_mapping = {}
  score = 0

  # This while loop is guaranteed to terminate since we increase both
  # struc_seq_idx and hit_seq_idx by at least 1 in each iteration.
  remaining_num_gaps = max_num_gaps
  while hit_seq_idx < len(hit_seq) and struc_seq_idx < len(struc_seq):
    if hit_seq[hit_seq_idx] != struc_seq[struc_seq_idx]:
      # Explore which alignment aligns the next residue (if present).
      best_shift = 0
      for shift in range(0, remaining_num_gaps + 1):
        next_hit_res = hit_seq[hit_seq_idx + shift : hit_seq_idx + shift + 1]
        next_struc_res = struc_seq[struc_seq_idx : struc_seq_idx + 1]
        if next_hit_res == next_struc_res:
          best_shift = shift
          break
      hit_seq_idx += best_shift
      remaining_num_gaps -= best_shift

    hit_to_struc_seq_mapping[hit_seq_idx] = struc_seq_idx
    score += hit_seq[hit_seq_idx] == struc_seq[struc_seq_idx]
    hit_seq_idx += 1
    struc_seq_idx += 1

  fixed_mapping = {}
  for query_idx, original_hit_idx in mapping.items():
    fixed_hit_idx = hit_to_struc_seq_mapping.get(original_hit_idx)
    if fixed_hit_idx is not None:
      fixed_mapping[query_idx] = fixed_hit_idx

  return fixed_mapping, score
