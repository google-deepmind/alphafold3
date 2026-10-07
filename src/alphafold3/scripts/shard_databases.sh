#!/bin/bash
# Copyright 2026 DeepMind Technologies Limited
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
#
# Shards AlphaFold 3 genetic databases in parallel using SeqKit round-robin
# splitting.

set -euo pipefail

readonly DB_SOURCE_DIR="${1:-$HOME/public_databases}"
readonly SHARDED_DIR="${2:-$DB_SOURCE_DIR/sharded}"
readonly TMP_DIR="${3:-$SHARDED_DIR/.tmp}"

if ! command -v seqkit > /dev/null 2>&1; then
  echo "seqkit is not installed. Please install it." >&2
  exit 1
fi

mkdir -p "${SHARDED_DIR}" "${TMP_DIR}"
trap 'rm -rf "${TMP_DIR}"' EXIT
# Increase open file descriptor limit to accommodate splitting into hundreds of shards in parallel.
ulimit -n 65535 2>/dev/null || true

shard_db() {
  local filename="$1"
  local num_shards="$2"

  if [[ ! -f "${DB_SOURCE_DIR}/${filename}" ]]; then
    echo "Warning: ${DB_SOURCE_DIR}/${filename} not found, skipping."
    return 0
  fi

  echo "Splitting ${filename} into ${num_shards} shards (round-robin)"
  mkdir -p "${TMP_DIR}/split_${filename}"
  seqkit split2 --threads 8 --by-part "${num_shards}" \
    "${DB_SOURCE_DIR}/${filename}" -O "${TMP_DIR}/split_${filename}"

  local total_shards_fmt
  printf -v total_shards_fmt "%05d" "${num_shards}"
  local idx=0
  for part in "${TMP_DIR}/split_${filename}/"*.part_*; do
    [[ -e "${part}" ]] || continue
    local shard_idx_fmt
    printf -v shard_idx_fmt "%05d" "${idx}"
    mv "${part}" "${SHARDED_DIR}/${filename}-${shard_idx_fmt}-of-${total_shards_fmt}"
    idx=$((idx + 1))
  done

  rm -rf "${TMP_DIR}/split_${filename}"

  if (( idx != num_shards )); then
    echo "Error: Expected ${num_shards} shards for ${filename}, but found ${idx}." >&2
    exit 1
  fi

  echo "Finished splitting ${filename} (${idx} shards created)"
}

pids=()
shard_db "bfd-first_non_consensus_sequences.fasta" 60 & pids+=($!)
shard_db "uniref90_2022_05.fa" 144 & pids+=($!)
shard_db "uniprot_all_2021_04.fa" 204 & pids+=($!)
shard_db "mgy_clusters_2022_05.fa" 564 & pids+=($!)
shard_db "rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta" 16 & pids+=($!)
shard_db "rnacentral_active_seq_id_90_cov_80_linclust.fasta" 64 & pids+=($!)
shard_db "nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta" 384 & pids+=($!)

for _ in "${pids[@]}"; do
  wait -n
done
