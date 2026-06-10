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

"""Functions for getting MSA and calculating alignment features."""

from collections.abc import MutableMapping, Sequence
import string
from typing import Self

from absl import logging
from alphafold3.constants import mmcif_names
from alphafold3.data import msa_config
from alphafold3.data import msa_features
from alphafold3.data import parsers
from alphafold3.data.tools import jackhmmer
from alphafold3.data.tools import msa_tool
from alphafold3.data.tools import nhmmer
import numpy as np


class Error(Exception):
  """Error indicating a problem with MSA Search."""


def _featurize(seq: str, chain_poly_type: str) -> str | list[int]:
  if mmcif_names.is_standard_polymer_type(chain_poly_type):
    featurized_seqs, _ = msa_features.extract_msa_features(
        msa_sequences=[seq], chain_poly_type=chain_poly_type
    )
    return featurized_seqs[0].tolist()
  # For anything else simply require an identical match.
  return seq


def sequences_are_feature_equivalent(
    sequence1: str,
    sequence2: str,
    chain_poly_type: str,
) -> bool:
  feat1 = _featurize(sequence1, chain_poly_type)
  feat2 = _featurize(sequence2, chain_poly_type)
  return feat1 == feat2


class Msa:
  """Multiple Sequence Alignment container with methods for manipulating it."""

  def __init__(
      self,
      query_sequence: str,
      chain_poly_type: str,
      sequences: Sequence[str],
      descriptions: Sequence[str],
      deduplicate: bool = True,
  ):
    """Raw constructor, prefer using the from_{a3m,multiple_msas} class methods.

    The first sequence must be equal (in featurised form) to the query sequence.
    If sequences/descriptions are empty, they will be initialised to the query.

    Args:
      query_sequence: The sequence that was used to search for MSA.
      chain_poly_type: Polymer type of the query sequence, see mmcif_names.
      sequences: The sequences returned by the MSA search tool.
      descriptions: Metadata for the sequences returned by the MSA search tool.
      deduplicate: If True, the MSA sequences will be deduplicated in the input
        order. Lowercase letters (insertions) are ignored when deduplicating.
    """
    if len(sequences) != len(descriptions):
      raise ValueError('The number of sequences and descriptions must match.')

    self.query_sequence = query_sequence
    self.chain_poly_type = chain_poly_type

    if not deduplicate:
      self.sequences = sequences
      self.descriptions = descriptions
    else:
      self.sequences = []
      self.descriptions = []
      # A replacement table that removes all lowercase characters.
      deletion_table = str.maketrans('', '', string.ascii_lowercase)
      unique_sequences = set()
      for seq, desc in zip(sequences, descriptions, strict=True):
        # Using string.translate is faster than re.sub('[a-z]+', '').
        sequence_no_deletions = seq.translate(deletion_table)
        if sequence_no_deletions not in unique_sequences:
          unique_sequences.add(sequence_no_deletions)
          self.sequences.append(seq)
          self.descriptions.append(desc)

    # Make sure the MSA always has at least the query.
    self.sequences = self.sequences or [query_sequence]
    self.descriptions = self.descriptions or ['Original query']

    # Check if the 1st MSA sequence matches the query sequence. Since it may be
    # mutated by the search tool (jackhmmer) check using the featurized version.
    if not sequences_are_feature_equivalent(
        self.sequences[0], query_sequence, chain_poly_type
    ):
      raise ValueError(
          f'First MSA sequence {self.sequences[0]} is not the {query_sequence=}'
      )

  @classmethod
  def from_multiple_msas(
      cls, msas: Sequence[Self], deduplicate: bool = True
  ) -> Self:
    """Initializes the MSA from multiple MSAs.

    Args:
      msas: A sequence of Msa objects representing individual MSAs produced by
        different tools/dbs.
      deduplicate: If True, the MSA sequences will be deduplicated in the input
        order. Lowercase letters (insertions) are ignored when deduplicating.

    Returns:
      An Msa object created by merging multiple MSAs.
    """
    if not msas:
      raise ValueError('At least one MSA must be provided.')

    query_sequence = msas[0].query_sequence
    chain_poly_type = msas[0].chain_poly_type
    sequences = []
    descriptions = []

    for msa in msas:
      if msa.query_sequence != query_sequence:
        raise ValueError(
            f'Query sequences must match: {[m.query_sequence for m in msas]}'
        )
      if msa.chain_poly_type != chain_poly_type:
        raise ValueError(
            f'Chain poly types must match: {[m.chain_poly_type for m in msas]}'
        )
      sequences.extend(msa.sequences)
      descriptions.extend(msa.descriptions)

    return cls(
        query_sequence=query_sequence,
        chain_poly_type=chain_poly_type,
        sequences=sequences,
        descriptions=descriptions,
        deduplicate=deduplicate,
    )

  @classmethod
  def from_multiple_a3ms(
      cls, a3ms: Sequence[str], chain_poly_type: str, deduplicate: bool = True
  ) -> Self:
    """Initializes the MSA from multiple A3M strings.

    Args:
      a3ms: A sequence of A3M strings representing individual MSAs produced by
        different tools/dbs.
      chain_poly_type: Polymer type of the query sequence, see mmcif_names.
      deduplicate: If True, the MSA sequences will be deduplicated in the input
        order. Lowercase letters (insertions) are ignored when deduplicating.

    Returns:
      An Msa object created by merging multiple A3Ms.
    """
    if not a3ms:
      raise ValueError('At least one A3M must be provided.')

    query_sequence = None
    all_sequences = []
    all_descriptions = []

    for a3m in a3ms:
      sequences, descriptions = parsers.parse_fasta(a3m)
      if query_sequence is None:
        query_sequence = sequences[0]

      if sequences[0] != query_sequence:
        raise ValueError(
            f'Query sequences must match: {sequences[0]=} != {query_sequence=}'
        )
      all_sequences.extend(sequences)
      all_descriptions.extend(descriptions)

    return cls(
        query_sequence=query_sequence,
        chain_poly_type=chain_poly_type,
        sequences=all_sequences,
        descriptions=all_descriptions,
        deduplicate=deduplicate,
    )

  @classmethod
  def from_a3m(
      cls,
      query_sequence: str,
      chain_poly_type: str,
      a3m: str,
      max_depth: int | None = None,
      deduplicate: bool = True,
  ) -> Self:
    """Parses the single A3M and builds the Msa object."""
    sequences, descriptions = parsers.parse_fasta(a3m)

    if max_depth is not None and 0 < max_depth < len(sequences):
      logging.info(
          'MSA cropped from depth of %d to %d for %s.',
          len(sequences),
          max_depth,
          query_sequence,
      )
      sequences = sequences[:max_depth]
      descriptions = descriptions[:max_depth]

    return cls(
        query_sequence=query_sequence,
        chain_poly_type=chain_poly_type,
        sequences=sequences,
        descriptions=descriptions,
        deduplicate=deduplicate,
    )

  @classmethod
  def from_empty(cls, query_sequence: str, chain_poly_type: str) -> Self:
    """Creates an empty Msa containing just the query sequence."""
    return cls(
        query_sequence=query_sequence,
        chain_poly_type=chain_poly_type,
        sequences=[],
        descriptions=[],
        deduplicate=False,
    )

  @property
  def depth(self) -> int:
    return len(self.sequences)

  def __repr__(self) -> str:
    return f'Msa({self.depth} sequences, {self.chain_poly_type})'

  def to_a3m(self) -> str:
    """Returns the MSA in the A3M format."""
    a3m_lines = []
    for desc, seq in zip(self.descriptions, self.sequences, strict=True):
      a3m_lines.append(f'>{desc}')
      a3m_lines.append(seq)
    return '\n'.join(a3m_lines) + '\n'

  def featurize(self) -> MutableMapping[str, np.ndarray]:
    """Featurises the MSA and returns a map of feature names to features.

    Returns:
      A dictionary mapping feature names to values.

    Raises:
      msa.Error:
        * If the sequences in the MSA don't have the same length after deletions
          (lower case letters) are removed.
        * If the MSA contains an unknown amino acid code.
        * If there are no sequences after aligning.
    """
    try:
      msa, deletion_matrix = msa_features.extract_msa_features(
          msa_sequences=self.sequences, chain_poly_type=self.chain_poly_type
      )
    except ValueError as e:
      raise Error(f'Error extracting MSA or deletion features: {e}') from e

    if msa.shape == (0, 0):
      raise Error(f'Empty MSA feature for {self}')

    species_ids = msa_features.extract_species_ids(self.descriptions)

    return {
        'msa_species_identifiers': np.array(species_ids, dtype=object),
        'num_alignments': np.array(self.depth, dtype=np.int32),
        'msa': msa,
        'deletion_matrix': deletion_matrix,
    }


def get_msa_tool(
    msa_tool_config: msa_config.JackhmmerConfig | msa_config.NhmmerConfig,
) -> msa_tool.MsaTool:
  """Returns the requested MSA tool."""

  match msa_tool_config:
    case msa_config.JackhmmerConfig():
      return jackhmmer.Jackhmmer(
          binary_path=msa_tool_config.binary_path,
          database_path=msa_tool_config.database_config.path,
          n_cpu=msa_tool_config.n_cpu,
          n_iter=msa_tool_config.n_iter,
          e_value=msa_tool_config.e_value,
          z_value=msa_tool_config.z_value,
          dom_z_value=msa_tool_config.dom_z_value,
          max_sequences=msa_tool_config.max_sequences,
          max_threads=msa_tool_config.max_parallel_shards,
      )
    case msa_config.NhmmerConfig():
      return nhmmer.Nhmmer(
          binary_path=msa_tool_config.binary_path,
          hmmalign_binary_path=msa_tool_config.hmmalign_binary_path,
          hmmbuild_binary_path=msa_tool_config.hmmbuild_binary_path,
          database_path=msa_tool_config.database_config.path,
          n_cpu=msa_tool_config.n_cpu,
          e_value=msa_tool_config.e_value,
          z_value=msa_tool_config.z_value,
          max_sequences=msa_tool_config.max_sequences,
          alphabet=msa_tool_config.alphabet,
          max_threads=msa_tool_config.max_parallel_shards,
      )
    case _:
      raise ValueError(f'Unknown MSA tool: {msa_tool_config}.')


def get_msa(
    target_sequence: str,
    run_config: msa_config.RunConfig,
    chain_poly_type: str,
    deduplicate: bool = False,
) -> Msa:
  """Computes the MSA for a given query sequence.

  Args:
    target_sequence: The target amino-acid sequence.
    run_config: MSA run configuration.
    chain_poly_type: The type of chain for which to get an MSA.
    deduplicate: If True, the MSA sequences will be deduplicated in the input
      order. Lowercase letters (insertions) are ignored when deduplicating.

  Returns:
    Aligned MSA sequences.
  """

  return Msa.from_a3m(
      query_sequence=target_sequence,
      chain_poly_type=chain_poly_type,
      a3m=get_msa_tool(run_config.config).query(target_sequence).a3m,
      max_depth=run_config.crop_size,
      deduplicate=deduplicate,
  )
