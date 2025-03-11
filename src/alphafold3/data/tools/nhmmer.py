# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Library to run Nhmmer from Python."""

import os
import pathlib
import tempfile
from typing import Final

from absl import logging
from alphafold3.data import parsers
from alphafold3.data.tools import hmmalign
from alphafold3.data.tools import hmmbuild
from alphafold3.data.tools import msa_tool
from alphafold3.data.tools import subprocess_utils

_SHORT_SEQUENCE_CUTOFF: Final[int] = 50


class Nhmmer(msa_tool.MsaTool):
  """Python wrapper of the Nhmmer binary."""

  def __init__(
      self,
      binary_path: str,
      hmmalign_binary_path: str,
      hmmbuild_binary_path: str,
      database_path: str,
      n_cpu: int = 8,
      e_value: float = 1e-3,
      max_sequences: int = 5000,
      filter_f3: float = 1e-5,
      alphabet: str | None = None,
      strand: str | None = None,
  ):
    """Initializes the Python Nhmmer wrapper.

    Args:
      binary_path: Path to the Nhmmer binary.
      hmmalign_binary_path: Path to the Hmmalign binary.
      hmmbuild_binary_path: Path to the Hmmbuild binary.
      database_path: MSA database path to search against. This can be either a
        FASTA (slow) or HMMERDB produced from the FASTA using the makehmmerdb
        binary. The HMMERDB is ~10x faster but experimental.
      n_cpu: The number of CPUs to give Nhmmer.
      e_value: The E-value, see Nhmmer docs for more details. Will be
        overwritten if bit_score is set.
      max_sequences: Maximum number of sequences to return in the MSA.
      filter_f3: Forward pre-filter, set to >1.0 to turn off.
      alphabet: The alphabet to assert when building a profile with hmmbuild.
        This must be 'rna', 'dna', or None.
      strand: "watson" searches query sequence, "crick" searches
        reverse-compliment and default is None which means searching for both.

    Raises:
      RuntimeError: If Nhmmer binary not found within the path.
    """
    self._binary_path = binary_path
    self._hmmalign_binary_path = hmmalign_binary_path
    self._hmmbuild_binary_path = hmmbuild_binary_path
    self._db_path = database_path

    subprocess_utils.check_binary_exists(path=self._binary_path, name='Nhmmer')

    if strand and strand not in {'watson', 'crick'}:
      raise ValueError(f'Invalid {strand=}. only "watson" or "crick" supported')

    if alphabet and alphabet not in {'rna', 'dna'}:
      raise ValueError(f'Invalid {alphabet=}, only "rna" or "dna" supported')

    self._e_value = e_value
    self._n_cpu = n_cpu
    self._max_sequences = max_sequences
    self._filter_f3 = filter_f3
    self._alphabet = alphabet
    self._strand = strand

  def query(self, target_sequence: str) -> msa_tool.MsaToolResult:
    """Query the database using Nhmmer."""
    logging.info('Query sequence: %s', target_sequence)

    with tempfile.TemporaryDirectory() as query_tmp_dir:
      input_a3m_path = os.path.join(query_tmp_dir, 'query.a3m')
      output_sto_path = os.path.join(query_tmp_dir, 'output.sto')
      pathlib.Path(output_sto_path).touch()
      subprocess_utils.create_query_fasta_file(
          sequence=target_sequence, path=input_a3m_path
      )

      cmd_flags = [
          *('-o', '/dev/null'),  # Don't pollute stdout with nhmmer output.
          '--noali',  # Don't include the alignment in stdout.
          *('--cpu', str(self._n_cpu)),
      ]

      cmd_flags.extend(['-E', str(self._e_value)])

      if self._alphabet:
        cmd_flags.extend([f'--{self._alphabet}'])

      if self._strand is not None:
        cmd_flags.extend([f'--{self._strand}'])

      cmd_flags.extend(['-A', output_sto_path])
      # As recommend by RNAcentral for short sequences.
      if (
          self._alphabet == 'rna'
          and len(target_sequence) < _SHORT_SEQUENCE_CUTOFF
      ):
        cmd_flags.extend(['--F3', str(0.02)])
      else:
        cmd_flags.extend(['--F3', str(self._filter_f3)])

      # The input A3M and the db are the last two arguments.
      cmd_flags.extend((input_a3m_path, self._db_path))

      cmd = [self._binary_path, *cmd_flags]

      subprocess_utils.run(
          cmd=cmd,
          cmd_name=f'Nhmmer ({os.path.basename(self._db_path)})',
          log_stdout=False,
          log_stderr=True,
          log_on_process_error=True,
      )

      if os.path.getsize(output_sto_path) > 0:
        with open(output_sto_path) as f:
          a3m_out = parsers.convert_stockholm_to_a3m(
              f, max_sequences=self._max_sequences - 1  # Query not included.
          )
        # Nhmmer hits are generally shorter than the query sequence. To get MSA
        # of width equal to the query sequence, align hits to the query profile.
        logging.info('Aligning output a3m of size %d bytes', len(a3m_out))

        aligner = hmmalign.Hmmalign(self._hmmalign_binary_path)
        target_sequence_fasta = f'>query\n{target_sequence}\n'
        profile_builder = hmmbuild.Hmmbuild(
            binary_path=self._hmmbuild_binary_path, alphabet=self._alphabet
        )
        profile = profile_builder.build_profile_from_a3m(target_sequence_fasta)
        a3m_out = aligner.align_sequences_to_profile(
            profile=profile, sequences_a3m=a3m_out
        )
        a3m_out = ''.join([target_sequence_fasta, a3m_out])

        # Parse the output a3m to remove line breaks.
        a3m = '\n'.join(
            [f'>{n}\n{s}' for s, n in parsers.lazy_parse_fasta_string(a3m_out)]
        )
      else:
        # Nhmmer returns an empty file if there are no hits.
        # In this case return only the query sequence.
        a3m = f'>query\n{target_sequence}'

    return msa_tool.MsaToolResult(
        target_sequence=target_sequence, e_value=self._e_value, a3m=a3m
    )
