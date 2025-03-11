# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Functions for parsing various file formats."""

from collections.abc import Iterable, Sequence
from typing import IO, TypeAlias

from alphafold3.cpp import fasta_iterator
from alphafold3.cpp import msa_conversion


DeletionMatrix: TypeAlias = Sequence[Sequence[int]]


def lazy_parse_fasta_string(fasta_string: str) -> Iterable[tuple[str, str]]:
  """Lazily parses a FASTA/A3M string and yields (sequence, description) tuples.

  This implementation is more memory friendly than `fasta_sequence` while
  offering comparable performance. The underlying implementation is in C++ and
  is therefore faster than a pure Python implementation.

  Use this method when parsing FASTA files where you already have the FASTA
  string, but need to control how far you iterate through its sequences.

  Arguments:
    fasta_string: A string with the contents of FASTA/A3M file.

  Returns:
    Iterator of (sequence, description). In the description, the leading ">" is
    stripped.

  Raises:
    ValueError if the FASTA/A3M file is invalid, e.g. empty.
  """

  # The lifetime of the FastaStringIterator is tied to the lifetime of
  # fasta_string - fasta_string must be kept while the iterator is in use.
  return fasta_iterator.FastaStringIterator(fasta_string)


def parse_fasta(fasta_string: str) -> tuple[Sequence[str], Sequence[str]]:
  """Parses FASTA string and returns list of strings with amino-acid sequences.

  Arguments:
    fasta_string: The string contents of a FASTA file.

  Returns:
    A tuple of two lists:
    * A list of sequences.
    * A list of sequence descriptions taken from the comment lines. In the
      same order as the sequences.
  """
  return fasta_iterator.parse_fasta_include_descriptions(fasta_string)


def convert_a3m_to_stockholm(a3m: str, max_seqs: int | None = None) -> str:
  """Converts MSA in the A3M format to the Stockholm format."""
  sequences, descriptions = parse_fasta(a3m)
  if max_seqs is not None:
    sequences = sequences[:max_seqs]
    descriptions = descriptions[:max_seqs]

  stockholm = ['# STOCKHOLM 1.0', '']

  # Add the Stockholm header with the sequence metadata.
  names = []
  for i, description in enumerate(descriptions):
    name, _, rest = description.replace('\t', ' ').partition(' ')
    # Ensure that the names are unique - stockholm format requires that
    # the sequence names are unique.
    name = f'{name}_{i}'
    names.append(name)
    # Avoid zero-length description due to historic hmmbuild parsing bug.
    desc = rest.strip() or '<EMPTY>'
    stockholm.append(f'#=GS {name.strip()} DE {desc}')
  stockholm.append('')

  # Convert insertions in a sequence into gaps in all other sequences that don't
  # have an insertion in that column as well.
  sequences = msa_conversion.convert_a3m_to_stockholm(sequences)

  # Add the MSA data.
  max_name_width = max(len(name) for name in names)
  for name, sequence in zip(names, sequences, strict=True):
    # Align the names to the left and pad with spaces to the maximum length.
    stockholm.append(f'{name:<{max_name_width}s} {sequence}')

  # Add the reference annotation for the query (the first sequence).
  ref_annotation = ''.join('.' if c == '-' else 'x' for c in sequences[0])
  stockholm.append(f'{"#=GC RF":<{max_name_width}s} {ref_annotation}')
  stockholm.append('//')

  return '\n'.join(stockholm)


def convert_stockholm_to_a3m(
    stockholm: IO[str],
    max_sequences: int | None = None,
    remove_first_row_gaps: bool = True,
    linewidth: int | None = None,
) -> str:
  """Converts MSA in Stockholm format to the A3M format."""
  descriptions = {}
  sequences = {}
  reached_max_sequences = False

  if linewidth is not None and linewidth <= 0:
    raise ValueError('linewidth must be > 0 or None')

  for line in stockholm:
    reached_max_sequences = max_sequences and len(sequences) >= max_sequences
    line = line.strip()
    # Ignore blank lines, markup and end symbols - remainder are alignment
    # sequence parts.
    if not line or line.startswith(('#', '//')):
      continue
    seqname, aligned_seq = line.split(maxsplit=1)
    if seqname not in sequences:
      if reached_max_sequences:
        continue
      sequences[seqname] = ''
    sequences[seqname] += aligned_seq

  if not sequences:
    return ''

  stockholm.seek(0)
  for line in stockholm:
    line = line.strip()
    if line[:4] == '#=GS':
      # Description row - example format is:
      # #=GS UniRef90_Q9H5Z4/4-78            DE [subseq from] cDNA: FLJ22755 ...
      columns = line.split(maxsplit=3)
      seqname, feature = columns[1:3]
      value = columns[3] if len(columns) == 4 else ''
      if feature != 'DE':
        continue
      if reached_max_sequences and seqname not in sequences:
        continue
      descriptions[seqname] = value
      if len(descriptions) == len(sequences):
        break

  assert len(descriptions) <= len(sequences)

  # Convert sto format to a3m line by line
  a3m_sequences = {}
  # query_sequence is assumed to be the first sequence
  query_sequence = next(iter(sequences.values()))
  for seqname, sto_sequence in sequences.items():
    if remove_first_row_gaps:
      a3m_sequences[seqname] = msa_conversion.align_sequence_to_gapless_query(
          sequence=sto_sequence, query_sequence=query_sequence
      ).replace('.', '')
    else:
      a3m_sequences[seqname] = sto_sequence.replace('.', '')

  fasta_chunks = []

  for seqname, seq in a3m_sequences.items():
    fasta_chunks.append(f'>{seqname} {descriptions.get(seqname, "")}')

    if linewidth:
      fasta_chunks.extend(
          seq[i : linewidth + i] for i in range(0, len(seq), linewidth)
      )
    else:
      fasta_chunks.append(seq)

  return '\n'.join(fasta_chunks) + '\n'  # Include terminating newline.
