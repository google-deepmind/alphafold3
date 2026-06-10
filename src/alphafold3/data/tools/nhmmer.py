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

"""Library to run Nhmmer from Python."""

from collections.abc import Iterable, Sequence
from concurrent import futures
import heapq
import os
import pathlib
import shutil
import tempfile
import time
from typing import Final

from absl import logging
from alphafold3.data import parsers
from alphafold3.data.tools import hmmalign
from alphafold3.data.tools import hmmbuild
from alphafold3.data.tools import msa_tool
from alphafold3.data.tools import shards
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
      z_value: float | int | None = None,
      max_sequences: int = 5000,
      filter_f3: float = 1e-5,
      alphabet: str | None = None,
      strand: str | None = None,
      max_threads: int | None = None,
  ):
    """Initializes the Python Nhmmer wrapper.

    NOTE: The MSA obtained by running against sharded dbs won't be always
    exactly the same as the MSA obtained by running against an unsharded db.
    This is because of Nhmmer deduplication logic, which won't spot duplicate
    hits across multiple shards. Usually this means that the sharded search
    finds more hits (likely bounded by the number of shards), but this should
    not pose an issue given how the results are used downstream. The problem is
    more pronounced with deep MSAs and lower in the hit list (higher e-values).

    Make sure to set the Z value when searching against a sharded database,
    otherwise the results won't match the normal unsharded search.

    Args:
      binary_path: Path to the Nhmmer binary.
      hmmalign_binary_path: Path to the Hmmalign binary.
      hmmbuild_binary_path: Path to the Hmmbuild binary.
      database_path: MSA database path to search against. This can be either a
        FASTA (slow) or HMMERDB produced from the FASTA using the makehmmerdb
        binary. The HMMERDB is ~10x faster but experimental.  Sharded file
        specs, e.g. <db_path>@<num_shards>, are supported.
      n_cpu: The number of CPUs to give Nhmmer.
      e_value: The E-value, see Nhmmer docs for more details. Will be
        overwritten if bit_score is set.
      z_value: The Z-value representing the number of comparisons done (i.e
        correct database size) for E-value calculation. Make sure to set this
        when searching against a sharded database, otherwise the e-values will
        be incorrectly scaled.
      max_sequences: Maximum number of sequences to return in the MSA.
      filter_f3: Forward pre-filter, set to >1.0 to turn off.
      alphabet: The alphabet to assert when building a profile with hmmbuild.
        This must be 'rna', 'dna', or None.
      strand: "watson" searches query sequence, "crick" searches
        reverse-compliment and default is None which means searching for both.
      max_threads: If given, the maximum number of threads used when running
        sharded databases.

    Raises:
      RuntimeError: If Nhmmer binary not found within the path.
      ValueError: If an invalid configuration is provided in the args.
    """
    self._database_path = database_path

    if shard_paths := shards.get_sharded_paths(self._database_path):
      if z_value is None:
        raise ValueError(
            'The Z-value must be set when searching against a sharded database '
            'to correctly scale e-values.'
        )
      if 'hmmerdb' in self._database_path:
        raise ValueError('HMMERDB is not supported in sharded mode.')

      if max_sequences <= 1:
        raise ValueError(
            'max_sequences must be greater than 1 when running in sharded '
            'mode, because each shard would return only the query sequence.'
        )

      self._shard_paths = shard_paths
      self._max_threads = len(self._shard_paths)
      if max_threads is not None:
        self._max_threads = min(max_threads, self._max_threads)
      logging.info('Nhmmer running with max_threads = %d', self._max_threads)
    else:
      self._shard_paths = None
      self._max_threads = None

    self._binary_path = binary_path
    self._hmmalign_binary_path = hmmalign_binary_path
    self._hmmbuild_binary_path = hmmbuild_binary_path
    subprocess_utils.check_binary_exists(path=self._binary_path, name='Nhmmer')

    if strand and strand not in {'watson', 'crick'}:
      raise ValueError(f'Invalid {strand=}. only "watson" or "crick" supported')

    if alphabet and alphabet not in {'rna', 'dna'}:
      raise ValueError(f'Invalid {alphabet=}, only "rna" or "dna" supported')

    self._e_value = e_value
    self._n_cpu = n_cpu
    self._z_value = z_value
    self._max_sequences = max_sequences
    self._filter_f3 = filter_f3
    self._alphabet = alphabet
    self._strand = strand

  def query(self, target_sequence: str) -> msa_tool.MsaToolResult:
    """Query the database (sharded or unsharded) using Nhmmer."""
    if self._shard_paths:
      # Sharded case, run the query against each database shard in parallel.
      logging.info(
          'Query sequence (sharded db): %s',
          target_sequence
          if len(target_sequence) <= 16
          else f'{target_sequence[:16]}... (len {len(target_sequence)})',
      )

      global_temp_dir = tempfile.mkdtemp()

      def _query_shard_fn(
          shard_path: str,
      ) -> tuple[msa_tool.MsaToolResult, float]:
        t_start = time.time()
        # Get tblout as it contains e-values we need for merging sequences.
        result = self._query_db_shard(
            target_sequence=target_sequence,
            db_shard_path=shard_path,
            get_tblout=True,  # Tblout contains e-values needed for merging.
            global_temp_dir=global_temp_dir,
        )
        return result, time.time() - t_start

      with futures.ThreadPoolExecutor(max_workers=self._max_threads) as ex:
        tool_outputs, timings = zip(*ex.map(_query_shard_fn, self._shard_paths))

      logging.info(
          'Finished query for %d shards, shard timings (seconds): %s',
          len(tool_outputs),
          ', '.join(f'{t:.1f}' for t in timings),
      )

      shutil.rmtree(global_temp_dir, ignore_errors=True)
      return _merge_nhmmer_results(tool_outputs, self._max_sequences)

    else:
      # Non-sharded case, run the query against the whole database.
      logging.info(
          'Query sequence (non-sharded db): %s',
          target_sequence
          if len(target_sequence) <= 16
          else f'{target_sequence[:16]}... (len {len(target_sequence)})',
      )
      return self._query_db_shard(
          target_sequence=target_sequence,
          db_shard_path=self._database_path,
          get_tblout=False,
      )

  def _query_db_shard(
      self,
      *,
      target_sequence: str,
      db_shard_path: str,
      get_tblout: bool,
      global_temp_dir: str | None = None,
  ) -> msa_tool.MsaToolResult:
    """Query the database shard using Nhmmer."""

    with tempfile.TemporaryDirectory(dir=global_temp_dir) as query_tmp_dir:
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

      if get_tblout:
        output_tblout_path = pathlib.Path(query_tmp_dir, 'tblout.txt')
        output_tblout_path.touch()
        cmd_flags.extend(['--tblout', str(output_tblout_path)])
      else:
        output_tblout_path = None

      cmd_flags.extend(['-E', str(self._e_value)])

      if self._z_value is not None:
        cmd_flags.extend(['-Z', str(self._z_value)])

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
      cmd_flags.extend((input_a3m_path, db_shard_path))

      cmd = [self._binary_path, *cmd_flags]
      subprocess_utils.run(
          cmd=cmd,
          cmd_name=f'Nhmmer ({os.path.basename(db_shard_path)})',
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

      # Get the tabular output which has e.g. e-value for each target.
      tbl = '' if output_tblout_path is None else output_tblout_path.read_text()

    return msa_tool.MsaToolResult(
        target_sequence=target_sequence,
        e_value=self._e_value,
        a3m=a3m,
        tblout=tbl,
    )


def _merge_nhmmer_results(
    nhmmer_results: Sequence[msa_tool.MsaToolResult],
    max_sequences: int,
) -> msa_tool.MsaToolResult:
  """Merges nhmmer result protos into a single one."""
  assert len(set(nh_res.target_sequence for nh_res in nhmmer_results)) == 1
  assert len(set(nh_res.e_value for nh_res in nhmmer_results)) == 1

  # Parse the TBL output, create a mapping from unique hit ID to TBL line.
  parsed_tbl = {}
  for nhmmer_result in nhmmer_results:
    assert nhmmer_result.tblout is not None
    for line in nhmmer_result.tblout.splitlines():
      if not line.startswith('#'):
        line_fields = line.split(maxsplit=15)
        accession = line_fields[0]
        alignment_from = line_fields[6]
        alignment_to = line_fields[7]
        # This is the unique ID that is used in the output A3M.
        unique_id = f'{accession}/{alignment_from}-{alignment_to}'
        parsed_tbl[unique_id] = line

  # Create an iterator and merge a3m info with tbl info.
  def _merged_a3m_tbl_iter(a3m: str) -> Iterable[tuple[str, str, str, str]]:
    # Don't parse the entire a3m, lazily parse only as many sequences as needed.
    iterator = iter(parsers.lazy_parse_fasta_string(a3m))
    next(iterator)  # Skip the query which isn't present in tblout.
    for sequence, description in iterator:
      name = description.partition(' ')[0]
      if tbl_info := parsed_tbl.get(name):
        # Skip sequences for which we don't have tbl information.
        yield sequence, description, tbl_info, name

  def sort_key(seq_data: tuple[str, str, str, str]) -> tuple[float, str]:
    unused_seq, unused_description, tbl_info, name = seq_data
    # Nucleic tblout has 16 space delimited columns. "-" used if no value
    # present. We want e-value in column 12, so do only 13 splits. Use the name
    # in case of an e-value tie.
    return float(tbl_info.split(maxsplit=13)[12]), name

  # A3M/TBL is sorted by e-value and name, hence we can merge them efficiently.
  merged_a3m_and_tblout = heapq.merge(
      *[_merged_a3m_tbl_iter(res.a3m) for res in nhmmer_results],
      key=sort_key,
  )

  # Truncate the a3m to max_sequences. Do not truncate the tblout.
  merged_tblout = []
  merged_a3m = [f'>query\n{nhmmer_results[0].target_sequence}']
  for seq, description, tbl_info, _ in merged_a3m_and_tblout:
    merged_tblout.append(tbl_info)
    if len(merged_a3m) < max_sequences:
      merged_a3m.append(f'>{description}\n{seq}')

  logging.info(
      'Limiting merged MSA depth from %d to %d',
      len(merged_tblout),
      max_sequences,
  )

  return msa_tool.MsaToolResult(
      target_sequence=nhmmer_results[0].target_sequence,
      a3m='\n'.join(merged_a3m),
      e_value=nhmmer_results[0].e_value,
      tblout=None,  # We no longer need the tblout.
  )
