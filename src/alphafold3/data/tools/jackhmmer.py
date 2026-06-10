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

"""Library to run Jackhmmer from Python."""

from collections.abc import Iterable, Sequence
from concurrent import futures
import heapq
import os
import pathlib
import shutil
import tempfile
import time

from absl import logging
from alphafold3.data import parsers
from alphafold3.data.tools import msa_tool
from alphafold3.data.tools import shards
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
      dom_e: float | None = None,
      dom_z_value: float | int | None = None,
      max_sequences: int = 5000,
      filter_f1: float = 5e-4,
      filter_f2: float = 5e-5,
      filter_f3: float = 5e-7,
      max_threads: int | None = None,
      **unused_kwargs,
  ):
    """Initializes the Python Jackhmmer wrapper.

    NOTE: The MSA obtained by running against sharded dbs won't be always
    exactly the same as the MSA obtained by running against an unsharded db.
    This is because of Jackhmmer deduplication logic, which won't spot duplicate
    hits across multiple shards. Usually this means that the sharded search
    finds more hits (likely bounded by the number of shards), but this should
    not pose an issue given how the results are used downstream. The problem is
    more pronounced with deep MSAs and lower in the hit list (higher e-values).

    Make sure to set the Z and domZ values when searching against a sharded
    database, otherwise the results won't match the normal unsharded search.

    Args:
      binary_path: The path to the jackhmmer executable.
      database_path: The path to the jackhmmer database (FASTA format). Sharded
        file specs, e.g. `<db_path>@<num_shards>`, are supported.
      n_cpu: The number of CPUs to give Jackhmmer.
      n_iter: The number of Jackhmmer iterations.
      e_value: The E-value, see Jackhmmer docs for more details.
      z_value: The Z-value representing the number of comparisons done (i.e
        correct database size) for E-value calculation. Make sure to set this
        when searching against a sharded database, otherwise the e-values will
        be incorrectly scaled.
      dom_e: Domain e-value criteria for inclusion in tblout.
      dom_z_value: Domain z-value representing the number of comparisons done
        (i.e correct database size) for domain E-value calculation. Make sure to
        set this when searching against a sharded database, otherwise the domain
        e-values will be incorrectly scaled.
      max_sequences: Maximum number of sequences to return in the MSA.
      filter_f1: MSV and biased composition pre-filter, set to >1.0 to turn off.
      filter_f2: Viterbi pre-filter, set to >1.0 to turn off.
      filter_f3: Forward pre-filter, set to >1.0 to turn off.
      max_threads: If given, the maximum number of threads used when running
        sharded databases.

    Raises:
      RuntimeError: If Jackhmmer binary not found within the path.
      ValueError: If an invalid configuration is provided in the args.
    """
    self._database_path = database_path

    if shard_paths := shards.get_sharded_paths(self._database_path):
      if n_iter != 1:
        raise ValueError('For a sharded db, only n_iter=1 is supported.')
      if z_value is None:
        raise ValueError(
            'The Z-value must be set when searching against a sharded database '
            'to correctly scale e-values.'
        )
      if max_sequences <= 1:
        raise ValueError(
            'max_sequences must be greater than 1 when running in sharded '
            'mode, because each shard would return only the query sequence.'
        )

      self._shard_paths = shard_paths
      self._max_threads = len(self._shard_paths)
      if max_threads is not None:
        self._max_threads = min(max_threads, self._max_threads)
      logging.info('Jackhmmer running with max_threads = %d', self._max_threads)
    else:
      self._shard_paths = None
      self._max_threads = None

    self._binary_path = binary_path
    subprocess_utils.check_binary_exists(
        path=self._binary_path, name='Jackhmmer'
    )

    self._n_cpu = n_cpu
    self._n_iter = n_iter
    self._e_value = e_value
    self._z_value = z_value
    self._dom_e = dom_e
    self._dom_z_value = dom_z_value
    self._max_sequences = max_sequences
    self._filter_f1 = filter_f1
    self._filter_f2 = filter_f2
    self._filter_f3 = filter_f3

    # If Jackhmmer supports the --seq_limit flag (via our patch), use it to
    # prevent writing out redundant sequences and increasing peak memory usage.
    # If not, the Jackhmmer will be run without the --seq_limit flag.
    self._supports_seq_limit = subprocess_utils.jackhmmer_seq_limit_supported(
        self._binary_path
    )

  def query(self, target_sequence: str) -> msa_tool.MsaToolResult:
    """Query the database (sharded or unsharded) using Jackhmmer."""
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
      return _merge_jackhmmer_results(tool_outputs, self._max_sequences)

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
    """Query the database shard using Jackhmmer."""

    with tempfile.TemporaryDirectory(dir=global_temp_dir) as query_tmp_dir:
      input_fasta_path = os.path.join(query_tmp_dir, 'query.fasta')
      subprocess_utils.create_query_fasta_file(
          sequence=target_sequence, path=input_fasta_path
      )
      output_sto_path = os.path.join(query_tmp_dir, 'output.sto')
      pathlib.Path(output_sto_path).touch()

      # The F1/F2/F3 are the expected proportion to pass each of the filtering
      # stages (which get progressively more expensive), reducing these
      # speeds up the pipeline at the expensive of sensitivity.  They are
      # currently set very low to make querying Mgnify run in a reasonable
      # amount of time.
      cmd_flags = [
          *('-o', '/dev/null'),  # Don't pollute stdout with Jackhmmer output.
          *('-A', output_sto_path),
          '--noali',
          *('--F1', str(self._filter_f1)),
          *('--F2', str(self._filter_f2)),
          *('--F3', str(self._filter_f3)),
          *('--cpu', str(self._n_cpu)),
          *('-N', str(self._n_iter)),
      ]

      if get_tblout:
        output_tblout_path = pathlib.Path(query_tmp_dir, 'tblout.txt')
        output_tblout_path.touch()
        cmd_flags.extend(['--tblout', str(output_tblout_path)])
      else:
        output_tblout_path = None

      # Report only sequences with E-values <= x in per-sequence output.
      if self._e_value is not None:
        cmd_flags.extend(['-E', str(self._e_value)])

        # Use the same value as the reporting e-value (`-E` flag).
        cmd_flags.extend(['--incE', str(self._e_value)])

      if self._z_value is not None:
        cmd_flags.extend(['-Z', str(self._z_value)])

      if self._dom_z_value is not None:
        cmd_flags.extend(['--domZ', str(self._dom_z_value)])

      if self._dom_e is not None:
        cmd_flags.extend(['--domE', str(self._dom_e)])

      if self._max_sequences is not None and self._supports_seq_limit:
        cmd_flags.extend(['--seq_limit', str(self._max_sequences)])

      # The input FASTA and the input db are the last two arguments.
      cmd = [self._binary_path] + cmd_flags + [input_fasta_path, db_shard_path]

      subprocess_utils.run(
          cmd=cmd,
          cmd_name=f'Jackhmmer ({os.path.basename(db_shard_path)})',
          log_stdout=False,
          log_stderr=True,
          log_on_process_error=True,
      )

      with open(output_sto_path) as f:
        a3m = parsers.convert_stockholm_to_a3m(
            f, max_sequences=self._max_sequences
        )

      # Get the tabular output which has e.g. e-value for each target.
      tbl = '' if output_tblout_path is None else output_tblout_path.read_text()

      return msa_tool.MsaToolResult(
          target_sequence=target_sequence,
          a3m=a3m,
          e_value=self._e_value,
          tblout=tbl,
      )


def _merge_jackhmmer_results(
    jh_results: Sequence[msa_tool.MsaToolResult], max_sequences: int
) -> msa_tool.MsaToolResult:
  """Merges Jackhmmer result protos into a single one."""
  assert len(set(jh_res.target_sequence for jh_res in jh_results)) == 1
  assert len(set(jh_res.e_value for jh_res in jh_results)) == 1

  # Parse the TBL output, create a mapping from hit name to TBL line.
  parsed_tbl = {}
  for jh_result in jh_results:
    assert jh_result.tblout is not None
    for line in jh_result.tblout.splitlines():
      if not line.startswith('#'):
        parsed_tbl[line.partition(' ')[0]] = line

  # Create an iterator and merge a3m info with tbl info.
  def _merged_a3m_tbl_iter(a3m: str) -> Iterable[tuple[str, str, str, str]]:
    # Don't parse the entire a3m, lazily parse only as many sequences as needed.
    iterator = iter(parsers.lazy_parse_fasta_string(a3m))
    next(iterator)  # Skip the query which isn't present in tblout.
    for sequence, description in iterator:
      name = description.partition(' ')[0].partition('/')[0]
      if tbl_info := parsed_tbl.get(name):
        # Skip sequences for which we don't have tbl information.
        yield sequence, description, tbl_info, name

  def sort_key(seq_data: tuple[str, str, str, str]) -> tuple[float, float, str]:
    unused_seq, unused_description, tbl_info, name = seq_data
    # Tblout lines have 19 whitespace delimited columns. "-" used if no value
    # present. We want e-value (column 5) and bit score (column 6), so do only 6
    # splits. E-value and bit score are equivalent, but bit score might have
    # higher resolution. Use the name in case of a tie.
    e_value, bit_score = tbl_info.split(maxsplit=6)[4:6]
    return float(e_value), -float(bit_score), name

  # A3M/TBL is sorted by e-value and name, hence we can merge them efficiently.
  merged_a3m_and_tblout = heapq.merge(
      *[_merged_a3m_tbl_iter(res.a3m) for res in jh_results],
      key=sort_key,
  )

  # Truncate the a3m to max_sequences. Do not truncate the tblout.
  merged_tblout = []
  merged_a3m = [f'>query\n{jh_results[0].target_sequence}']
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
      target_sequence=jh_results[0].target_sequence,
      a3m='\n'.join(merged_a3m),
      e_value=jh_results[0].e_value,
      tblout=None,  # We no longer need the tblout.
  )
