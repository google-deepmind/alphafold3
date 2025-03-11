# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Library to run Jackhmmer from Python."""

import os
import tempfile

from absl import logging
from alphafold3.data import parsers
from alphafold3.data.tools import msa_tool
from alphafold3.data.tools import subprocess_utils


class Jackhmmer(msa_tool.MsaTool):
  """Python wrapper of the Jackhmmer binary."""

  def __init__(
      self,
      *,
      binary_path: str,
      database_path: str,
      n_cpu: int = 8,
      n_iter: int = 3,
      e_value: float | None = 1e-3,
      z_value: float | int | None = None,
      max_sequences: int = 5000,
      filter_f1: float = 5e-4,
      filter_f2: float = 5e-5,
      filter_f3: float = 5e-7,
  ):
    """Initializes the Python Jackhmmer wrapper.

    Args:
      binary_path: The path to the jackhmmer executable.
      database_path: The path to the jackhmmer database (FASTA format).
      n_cpu: The number of CPUs to give Jackhmmer.
      n_iter: The number of Jackhmmer iterations.
      e_value: The E-value, see Jackhmmer docs for more details.
      z_value: The Z-value representing the number of comparisons done (i.e
        correct database size) for E-value calculation.
      max_sequences: Maximum number of sequences to return in the MSA.
      filter_f1: MSV and biased composition pre-filter, set to >1.0 to turn off.
      filter_f2: Viterbi pre-filter, set to >1.0 to turn off.
      filter_f3: Forward pre-filter, set to >1.0 to turn off.

    Raises:
      RuntimeError: If Jackhmmer binary not found within the path.
    """
    self.binary_path = binary_path
    self.database_path = database_path

    subprocess_utils.check_binary_exists(
        path=self.binary_path, name='Jackhmmer'
    )

    if not os.path.exists(self.database_path):
      raise ValueError(f'Could not find Jackhmmer database {database_path}')

    self.n_cpu = n_cpu
    self.n_iter = n_iter
    self.e_value = e_value
    self.z_value = z_value
    self.max_sequences = max_sequences
    self.filter_f1 = filter_f1
    self.filter_f2 = filter_f2
    self.filter_f3 = filter_f3

  def query(self, target_sequence: str) -> msa_tool.MsaToolResult:
    """Queries the database using Jackhmmer."""
    logging.info('Query sequence: %s', target_sequence)
    with tempfile.TemporaryDirectory() as query_tmp_dir:
      input_fasta_path = os.path.join(query_tmp_dir, 'query.fasta')
      subprocess_utils.create_query_fasta_file(
          sequence=target_sequence, path=input_fasta_path
      )

      output_sto_path = os.path.join(query_tmp_dir, 'output.sto')

      # The F1/F2/F3 are the expected proportion to pass each of the filtering
      # stages (which get progressively more expensive), reducing these
      # speeds up the pipeline at the expensive of sensitivity.  They are
      # currently set very low to make querying Mgnify run in a reasonable
      # amount of time.
      cmd_flags = [
          *('-o', '/dev/null'),  # Don't pollute stdout with Jackhmmer output.
          *('-A', output_sto_path),
          '--noali',
          *('--F1', str(self.filter_f1)),
          *('--F2', str(self.filter_f2)),
          *('--F3', str(self.filter_f3)),
          *('--cpu', str(self.n_cpu)),
          *('-N', str(self.n_iter)),
      ]

      # Report only sequences with E-values <= x in per-sequence output.
      if self.e_value is not None:
        cmd_flags.extend(['-E', str(self.e_value)])

        # Use the same value as the reporting e-value (`-E` flag).
        cmd_flags.extend(['--incE', str(self.e_value)])

      if self.z_value is not None:
        cmd_flags.extend(['-Z', str(self.z_value)])

      cmd = (
          [self.binary_path]
          + cmd_flags
          + [input_fasta_path, self.database_path]
      )

      subprocess_utils.run(
          cmd=cmd,
          cmd_name=f'Jackhmmer ({os.path.basename(self.database_path)})',
          log_stdout=False,
          log_stderr=True,
          log_on_process_error=True,
      )

      with open(output_sto_path) as f:
        a3m = parsers.convert_stockholm_to_a3m(
            f, max_sequences=self.max_sequences
        )

    return msa_tool.MsaToolResult(
        target_sequence=target_sequence, a3m=a3m, e_value=self.e_value
    )
