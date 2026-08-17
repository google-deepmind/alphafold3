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

"""A Python wrapper for hmmalign from the HMMER Suite."""

from collections.abc import Mapping, Sequence
import os
import tempfile

from alphafold3.data import parsers
from alphafold3.data.tools import subprocess_utils


def _to_a3m(sequences: Sequence[str], name_prefix: str = 'sequence') -> str:
  a3m = ''
  for i, sequence in enumerate(sequences, 1):
    a3m += f'> {name_prefix} {i}\n{sequence}\n'
  return a3m


class Hmmalign:
  """Python wrapper of the hmmalign binary."""

  def __init__(self, binary_path: str):
    """Initializes the Python hmmalign wrapper.

    Args:
      binary_path: Path to the hmmalign binary.

    Raises:
      RuntimeError: If hmmalign binary not found within the path.
    """
    self._binary_path = binary_path

    subprocess_utils.check_binary_exists(
        path=self._binary_path, name='hmmalign'
    )

  def align_sequences(
      self,
      sequences: Sequence[str],
      profile: str,
      extra_flags: Mapping[str, str] | None = None,
  ) -> str:
    """Aligns sequence list to the profile and returns the alignment in A3M."""
    return self.align(
        a3m_str=_to_a3m(sequences, name_prefix='query'),
        profile=profile,
        extra_flags=extra_flags,
    )

  def align(
      self,
      a3m_str: str,
      profile: str,
      extra_flags: Mapping[str, str] | None = None,
  ) -> str:
    """Aligns sequences in A3M to the profile and returns the alignment in A3M.

    Args:
      a3m_str: A list of sequence strings.
      profile: A hmm file with the hmm profile to align the sequences to.
      extra_flags: Dictionary with extra flags, flag_name: flag_value, that are
        added to hmmalign.

    Returns:
      An A3M string with the aligned sequences.

    Raises:
      RuntimeError: If hmmalign fails.
    """
    with tempfile.TemporaryDirectory() as query_tmp_dir:
      input_profile = os.path.join(query_tmp_dir, 'profile.hmm')
      input_sequences = os.path.join(query_tmp_dir, 'sequences.a3m')
      output_a3m_path = os.path.join(query_tmp_dir, 'output.a3m')

      with open(input_profile, 'w') as f:
        f.write(profile)

      with open(input_sequences, 'w') as f:
        f.write(a3m_str)

      cmd = [
          self._binary_path,
          *('-o', output_a3m_path),
          *('--outformat', 'A2M'),  # A2M is A3M in the HMMER suite.
      ]
      if extra_flags:
        for flag_name, flag_value in extra_flags.items():
          cmd.extend([flag_name, flag_value])
      cmd.extend([input_profile, input_sequences])

      subprocess_utils.run(
          cmd=cmd,
          cmd_name='hmmalign',
          log_stdout=False,
          log_stderr=True,
          log_on_process_error=True,
      )

      with open(output_a3m_path, encoding='utf-8') as f:
        a3m = f.read()

    return a3m

  def align_sequences_to_profile(self, profile: str, sequences_a3m: str) -> str:
    """Aligns the sequences to profile and returns the alignment in A3M string.

    Uses hmmalign to align the sequences to the profile, then ouputs the
    sequence contatenated at the beginning of the sequences in the A3M format.
    As the sequences are represented by an alignment with possible gaps ('-')
    and insertions (lowercase characters), the method first removes the gaps,
    then uppercases the insertions to prepare the sequences for realignment.
    Sequences with gaps cannot be aligned, as '-'s are not a valid symbol to
    align; lowercase characters must be uppercased to preserve the original
    sequences before realignment.

    Args:
      profile: The Hmmbuild profile to align the sequences to.
      sequences_a3m: Sequences in A3M format to align to the profile.

    Returns:
      An A3M string with the aligned sequences.

    Raises:
      RuntimeError: If hmmalign fails.
    """
    deletion_table = str.maketrans('', '', '-')
    sequences_no_gaps_a3m = []
    for seq, desc in parsers.lazy_parse_fasta_string(sequences_a3m):
      sequences_no_gaps_a3m.append(f'>{desc}')
      sequences_no_gaps_a3m.append(seq.translate(deletion_table))
    sequences_no_gaps_a3m = '\n'.join(sequences_no_gaps_a3m)

    aligned_sequences = self.align(sequences_no_gaps_a3m, profile)

    return aligned_sequences
