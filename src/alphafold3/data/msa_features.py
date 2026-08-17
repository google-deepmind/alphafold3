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

"""Utilities for computing MSA features."""

from collections.abc import Sequence
import re
from alphafold3.constants import mmcif_names
import numpy as np

_PROTEIN_TO_ID = {
    'A': 0,
    'B': 3,  # Same as D.
    'C': 4,
    'D': 3,
    'E': 6,
    'F': 13,
    'G': 7,
    'H': 8,
    'I': 9,
    'J': 20,  # Same as unknown (X).
    'K': 11,
    'L': 10,
    'M': 12,
    'N': 2,
    'O': 20,  # Same as unknown (X).
    'P': 14,
    'Q': 5,
    'R': 1,
    'S': 15,
    'T': 16,
    'U': 4,  # Same as C.
    'V': 19,
    'W': 17,
    'X': 20,
    'Y': 18,
    'Z': 6,  # Same as E.
    '-': 21,
}

_RNA_TO_ID = {
    # Map non-standard residues to UNK_NUCLEIC (N) -> 30
    **{chr(i): 30 for i in range(ord('A'), ord('Z') + 1)},
    # Continue the RNA indices from where Protein indices left off.
    '-': 21,
    'A': 22,
    'G': 23,
    'C': 24,
    'U': 25,
}

_DNA_TO_ID = {
    # Map non-standard residues to UNK_NUCLEIC (N) -> 30
    **{chr(i): 30 for i in range(ord('A'), ord('Z') + 1)},
    # Continue the DNA indices from where DNA indices left off.
    '-': 21,
    'A': 26,
    'G': 27,
    'C': 28,
    'T': 29,
}


def extract_msa_features(
    msa_sequences: Sequence[str], chain_poly_type: str
) -> tuple[np.ndarray, np.ndarray]:
  """Extracts MSA features.

  Example:
  The input raw MSA is: `[["AAAAAA"], ["Ai-CiDiiiEFa"]]`
  The output MSA will be: `[["AAAAAA"], ["A-CDEF"]]`
  The deletions will be: `[[0, 0, 0, 0, 0, 0], [0, 1, 0, 1, 3, 0]]`

  Args:
    msa_sequences: A list of strings, each string with one MSA sequence. Each
      string must have the same, constant number of non-lowercase (matching)
      residues.
    chain_poly_type: Either 'polypeptide(L)' (protein), 'polyribonucleotide'
      (RNA), or 'polydeoxyribonucleotide' (DNA). Use the appropriate string
      constant from mmcif_names.py.

  Returns:
    A tuple with:
    * MSA array of shape (num_seq, num_res) that contains only the uppercase
      characters or gaps (-) from the original MSA.
    * Deletions array of shape (num_seq, num_res) that contains the number
      of deletions (lowercase letters in the MSA) to the left from each
      non-deleted residue (uppercase letters in the MSA).

  Raises:
    ValueError if any of the preconditions are not met.
  """

  # Select the appropriate character map based on the chain type.
  if chain_poly_type == mmcif_names.RNA_CHAIN:
    char_map = _RNA_TO_ID
  elif chain_poly_type == mmcif_names.DNA_CHAIN:
    char_map = _DNA_TO_ID
  elif chain_poly_type == mmcif_names.PROTEIN_CHAIN:
    char_map = _PROTEIN_TO_ID
  else:
    raise ValueError(f'{chain_poly_type=} invalid.')

  # Handle empty MSA.
  if not msa_sequences:
    empty_msa = np.array([], dtype=np.int32).reshape((0, 0))
    empty_deletions = np.array([], dtype=np.int32).reshape((0, 0))
    return empty_msa, empty_deletions

  # Get the number of rows and columns in the MSA.
  num_rows = len(msa_sequences)
  num_cols = sum(1 for c in msa_sequences[0] if c in char_map)

  # Initialize the output arrays.
  msa_arr = np.zeros((num_rows, num_cols), dtype=np.int32)
  deletions_arr = np.zeros((num_rows, num_cols), dtype=np.int32)

  # Populate the output arrays.
  for problem_row, msa_sequence in enumerate(msa_sequences):
    deletion_count = 0
    upper_count = 0
    problem_col = 0
    problems = []
    for current in msa_sequence:
      msa_id = char_map.get(current, -1)
      if msa_id == -1:
        if not current.islower():
          problems.append(f'({problem_row}, {problem_col}):{current}')
        deletion_count += 1
      else:
        # Check the access is safe before writing to the array.
        # We don't need to check problem_row since it's guaranteed to be within
        # the array bounds, while upper_count is incremented in the loop.
        if upper_count < deletions_arr.shape[1]:
          deletions_arr[problem_row, upper_count] = deletion_count
          msa_arr[problem_row, upper_count] = msa_id
        deletion_count = 0
        upper_count += 1
      problem_col += 1
    if problems:
      raise ValueError(
          f"Unknown residues in MSA: {', '.join(problems)}. "
          f'target_sequence: {msa_sequences[0]}'
      )
    if upper_count != num_cols:
      raise ValueError(
          'Invalid shape all strings must have the same number '
          'of non-lowercase characters; First string has '
          f"{num_cols} non-lowercase characters but '{msa_sequence}' has "
          f'{upper_count}. target_sequence: {msa_sequences[0]}'
      )

  return msa_arr, deletions_arr


# UniProtKB SwissProt/TrEMBL dbs have the following description format:
# `db|UniqueIdentifier|EntryName`, e.g. `sp|P0C2L1|A3X1_LOXLA` or
# `tr|A0A146SKV9|A0A146SKV9_FUNHE`.
_UNIPROT_ENTRY_NAME_REGEX = re.compile(
    # UniProtKB TrEMBL or SwissProt database.
    r'(?:tr|sp)\|'
    # A primary accession number of the UniProtKB entry.
    r'(?:[A-Z0-9]{6,10})'
    # Occasionally there is an isoform suffix (e.g. _1 or _10) which we ignore.
    r'(?:_\d+)?\|'
    # TrEMBL: Same as AccessionId (6-10 characters).
    # SwissProt: A mnemonic protein identification code (1-5 characters).
    r'(?:[A-Z0-9]{1,10}_)'
    # A mnemonic species identification code.
    r'(?P<SpeciesId>[A-Z0-9]{1,5})'
)


def extract_species_ids(msa_descriptions: Sequence[str]) -> Sequence[str]:
  """Extracts species ID from MSA UniProtKB sequence identifiers.

  Args:
    msa_descriptions: The descriptions (the FASTA/A3M comment line) for each of
      the sequences.

  Returns:
    Extracted UniProtKB species IDs if there is a regex match for each
    description line, blank if the regex doesn't match.
  """
  species_ids = []
  for msa_description in msa_descriptions:
    msa_description = msa_description.strip()
    match = _UNIPROT_ENTRY_NAME_REGEX.match(msa_description)
    if match:
      species_ids.append(match.group('SpeciesId'))
    else:
      # Handle cases where the regex doesn't match
      # (e.g., append None or raise an error depending on your needs)
      species_ids.append('')
  return species_ids
