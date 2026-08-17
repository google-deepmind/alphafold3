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

"""A Python wrapper for hmmbuild - construct HMM profiles from MSA."""

import os
import re
import tempfile
from typing import Literal

from alphafold3.data import parsers
from alphafold3.data.tools import subprocess_utils


class Hmmbuild(object):
  """Python wrapper of the hmmbuild binary."""

  def __init__(
      self,
      *,
      binary_path: str,
      singlemx: bool = False,
      alphabet: str | None = None,
  ):
    """Initializes the Python hmmbuild wrapper.

    Args:
      binary_path: The path to the hmmbuild executable.
      singlemx: Whether to use --singlemx flag. If True, it forces HMMBuild to
        just use a common substitution score matrix.
      alphabet: The alphabet to assert when building a profile. Useful when
        hmmbuild cannot guess the alphabet. If None, no alphabet is asserted.

    Raises:
      RuntimeError: If hmmbuild binary not found within the path.
    """
    self._binary_path = binary_path
    self._singlemx = singlemx
    self._alphabet = alphabet

    subprocess_utils.check_binary_exists(
        path=self._binary_path, name='hmmbuild'
    )

  def build_profile_from_sto(self, sto: str, model_construction='fast') -> str:
    """Builds a HHM for the aligned sequences given as an A3M string.

    Args:
      sto: A string with the aligned sequences in the Stockholm format.
      model_construction: Whether to use reference annotation in the msa to
        determine consensus columns ('hand') or default ('fast').

    Returns:
      A string with the profile in the HMM format.

    Raises:
      RuntimeError: If hmmbuild fails.
    """
    return self._build_profile(
        sto, informat='stockholm', model_construction=model_construction
    )

  def build_profile_from_a3m(self, a3m: str) -> str:
    """Builds a HHM for the aligned sequences given as an A3M string.

    Args:
      a3m: A string with the aligned sequences in the A3M format.

    Returns:
      A string with the profile in the HMM format.

    Raises:
      RuntimeError: If hmmbuild fails.
    """
    lines = []
    for sequence, description in parsers.lazy_parse_fasta_string(a3m):
      sequence = re.sub('[a-z]+', '', sequence)  # Remove inserted residues.
      lines.append(f'>{description}\n{sequence}\n')
    msa = ''.join(lines)
    return self._build_profile(msa, informat='afa')

  def _build_profile(
      self,
      msa: str,
      informat: Literal['afa', 'stockholm'],
      model_construction: str = 'fast',
  ) -> str:
    """Builds a HMM for the aligned sequences given as an MSA string.

    Args:
      msa: A string with the aligned sequences, in A3M or STO format.
      informat: One of 'afa' (aligned FASTA) or 'sto' (Stockholm).
      model_construction: Whether to use reference annotation in the msa to
        determine consensus columns ('hand') or default ('fast').

    Returns:
      A string with the profile in the HMM format.

    Raises:
      RuntimeError: If hmmbuild fails.
      ValueError: If unspecified arguments are provided.
    """
    if model_construction not in {'hand', 'fast'}:
      raise ValueError(f'Bad {model_construction=}. Only hand or fast allowed.')

    with tempfile.TemporaryDirectory() as query_tmp_dir:
      input_msa_path = os.path.join(query_tmp_dir, 'query.msa')
      output_hmm_path = os.path.join(query_tmp_dir, 'output.hmm')

      with open(input_msa_path, 'w') as f:
        f.write(msa)

      # Specify the format as we don't specify the input file extension. See
      # https://github.com/EddyRivasLab/hmmer/issues/321 for more details.
      cmd_flags = ['--informat', informat]
      # If adding flags, we have to do so before the output and input:
      if model_construction == 'hand':
        cmd_flags.append(f'--{model_construction}')
      if self._singlemx:
        cmd_flags.append('--singlemx')
      if self._alphabet:
        cmd_flags.append(f'--{self._alphabet}')

      cmd_flags.extend([output_hmm_path, input_msa_path])

      cmd = [self._binary_path, *cmd_flags]

      subprocess_utils.run(
          cmd=cmd,
          cmd_name='Hmmbuild',
          log_stdout=False,
          log_stderr=True,
          log_on_process_error=True,
      )

      with open(output_hmm_path) as f:
        hmm = f.read()

    return hmm
