# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Helper functions for launching external tools."""

from collections.abc import Sequence
import os
import subprocess
import time
from typing import Any

from absl import logging


def create_query_fasta_file(sequence: str, path: str, linewidth: int = 80):
  """Creates a fasta file with the sequence with line width limit."""
  with open(path, 'w') as f:
    f.write('>query\n')

    i = 0
    while i < len(sequence):
      f.write(f'{sequence[i:(i + linewidth)]}\n')
      i += linewidth


def check_binary_exists(path: str, name: str) -> None:
  """Checks if a binary exists on the given path and raises otherwise."""
  if not os.path.exists(path):
    raise RuntimeError(f'{name} binary not found at {path}')


def jackhmmer_seq_limit_supported(jackhmmer_path: str) -> bool:
  """Checks if Jackhmmer supports the --seq-limit flag."""
  try:
    subprocess.run([jackhmmer_path, '-h', '--seq_limit', '1'], capture_output=True, check=True)
  except subprocess.CalledProcessError:
    return False
  return True


def run(
    cmd: Sequence[str],
    cmd_name: str,
    log_on_process_error: bool = False,
    log_stderr: bool = False,
    log_stdout: bool = False,
    max_out_streams_len: int | None = 500_000,
    **run_kwargs,
) -> subprocess.CompletedProcess[Any]:
  """Launches a subprocess, times it, and checks for errors.

  Args:
    cmd: Command to launch.
    cmd_name: Human-readable command name to be used in logs.
    log_on_process_error: Whether to use `logging.error` to log the process'
      stderr on failure.
    log_stderr: Whether to log the stderr of the command.
    log_stdout: Whether to log the stdout of the command.
    max_out_streams_len: Max length of prefix of stdout and stderr included in
      the exception message. Set to `None` to disable truncation.
    **run_kwargs: Any other kwargs for `subprocess.run`.

  Returns:
    The completed process object.

  Raises:
    RuntimeError: if the process completes with a non-zero return code.
  """

  logging.info('Launching subprocess "%s"', ' '.join(cmd))

  start_time = time.time()
  try:
    completed_process = subprocess.run(
        cmd,
        check=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True,
        **run_kwargs,
    )
  except subprocess.CalledProcessError as e:
    if log_on_process_error:
      # Logs have a 15k character limit, so log the error line by line.
      logging.error('%s failed. %s stderr begin:', cmd_name, cmd_name)
      for error_line in e.stderr.splitlines():
        if stripped_error_line := error_line.strip():
          logging.error(stripped_error_line)
      logging.error('%s stderr end.', cmd_name)

    error_msg = (
        f'{cmd_name} failed'
        f'\nstdout:\n{e.stdout[:max_out_streams_len]}\n'
        f'\nstderr:\n{e.stderr[:max_out_streams_len]}'
    )
    raise RuntimeError(error_msg) from e
  end_time = time.time()

  logging.info('Finished %s in %.3f seconds', cmd_name, end_time - start_time)
  stdout, stderr = completed_process.stdout, completed_process.stderr

  if log_stdout and stdout:
    logging.info('%s stdout:\n%s', cmd_name, stdout)

  if log_stderr and stderr:
    logging.info('%s stderr:\n%s', cmd_name, stderr)

  return completed_process
