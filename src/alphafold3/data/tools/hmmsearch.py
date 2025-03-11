# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""A Python wrapper for hmmsearch - search profile against a sequence db."""

import os
import tempfile

from absl import logging
from alphafold3.data import parsers
from alphafold3.data.tools import hmmbuild
from alphafold3.data.tools import subprocess_utils


class Hmmsearch(object):
  """Python wrapper of the hmmsearch binary."""

  def __init__(
      self,
      *,
      binary_path: str,
      hmmbuild_binary_path: str,
      database_path: str,
      alphabet: str = 'amino',
      filter_f1: float | None = None,
      filter_f2: float | None = None,
      filter_f3: float | None = None,
      e_value: float | None = None,
      inc_e: float | None = None,
      dom_e: float | None = None,
      incdom_e: float | None = None,
      filter_max: bool = False,
  ):
    """Initializes the Python hmmsearch wrapper.

    Args:
      binary_path: The path to the hmmsearch executable.
      hmmbuild_binary_path: The path to the hmmbuild executable. Used to build
        an hmm from an input a3m.
      database_path: The path to the hmmsearch database (FASTA format).
      alphabet: Chain type e.g. amino, rna, dna.
      filter_f1: MSV and biased composition pre-filter, set to >1.0 to turn off.
      filter_f2: Viterbi pre-filter, set to >1.0 to turn off.
      filter_f3: Forward pre-filter, set to >1.0 to turn off.
      e_value: E-value criteria for inclusion in tblout.
      inc_e: E-value criteria for inclusion in MSA/next round.
      dom_e: Domain e-value criteria for inclusion in tblout.
      incdom_e: Domain e-value criteria for inclusion of domains in MSA/next
        round.
      filter_max: Remove all filters, will ignore all filter_f* settings.

    Raises:
      RuntimeError: If hmmsearch binary not found within the path.
    """
    self.binary_path = binary_path
    self.hmmbuild_runner = hmmbuild.Hmmbuild(
        alphabet=alphabet, binary_path=hmmbuild_binary_path
    )
    self.database_path = database_path
    flags = []
    if filter_max:
      flags.append('--max')
    else:
      if filter_f1 is not None:
        flags.extend(('--F1', filter_f1))
      if filter_f2 is not None:
        flags.extend(('--F2', filter_f2))
      if filter_f3 is not None:
        flags.extend(('--F3', filter_f3))

    if e_value is not None:
      flags.extend(('-E', e_value))
    if inc_e is not None:
      flags.extend(('--incE', inc_e))
    if dom_e is not None:
      flags.extend(('--domE', dom_e))
    if incdom_e is not None:
      flags.extend(('--incdomE', incdom_e))

    self.flags = tuple(map(str, flags))

    subprocess_utils.check_binary_exists(
        path=self.binary_path, name='hmmsearch'
    )

    if not os.path.exists(self.database_path):
      logging.error('Could not find hmmsearch database %s', database_path)
      raise ValueError(f'Could not find hmmsearch database {database_path}')

  def query_with_hmm(self, hmm: str) -> str:
    """Queries the database using hmmsearch using a given hmm."""
    with tempfile.TemporaryDirectory() as query_tmp_dir:
      hmm_input_path = os.path.join(query_tmp_dir, 'query.hmm')
      sto_out_path = os.path.join(query_tmp_dir, 'output.sto')
      with open(hmm_input_path, 'w') as f:
        f.write(hmm)

      cmd = [
          self.binary_path,
          '--noali',  # Don't include the alignment in stdout.
          *('--cpu', '8'),
      ]
      # If adding flags, we have to do so before the output and input:
      if self.flags:
        cmd.extend(self.flags)
      cmd.extend([
          *('-A', sto_out_path),
          hmm_input_path,
          self.database_path,
      ])

      subprocess_utils.run(
          cmd=cmd,
          cmd_name=f'Hmmsearch ({os.path.basename(self.database_path)})',
          log_stdout=False,
          log_stderr=True,
          log_on_process_error=True,
      )

      with open(sto_out_path) as f:
        a3m_out = parsers.convert_stockholm_to_a3m(
            f, remove_first_row_gaps=False, linewidth=60
        )

    return a3m_out

  def query_with_a3m(self, a3m_in: str) -> str:
    """Query the database using hmmsearch using a given a3m."""

    # Only the "fast" model construction makes sense with A3M, as it doesn't
    # have any way to annotate reference columns.
    hmm = self.hmmbuild_runner.build_profile_from_a3m(a3m_in)
    return self.query_with_hmm(hmm)

  def query_with_sto(
      self, msa_sto: str, model_construction: str = 'fast'
  ) -> str:
    """Queries the database using hmmsearch using a given stockholm msa."""
    hmm = self.hmmbuild_runner.build_profile_from_sto(
        msa_sto, model_construction=model_construction
    )
    return self.query_with_hmm(hmm)
