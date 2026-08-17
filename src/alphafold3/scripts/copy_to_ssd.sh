#!/bin/bash
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

set -euo pipefail

readonly SOURCE_DIR=${1:-$HOME/public_databases}
readonly TARGET_DIR=${2:-/mnt/disks/ssd/public_databases}

mkdir -p "${TARGET_DIR}"

FILES=(pdb_seqres_2022_09_28.fasta \
      uniprot_all_2021_04.fa \
      mgy_clusters_2022_05.fa \
      uniref90_2022_05.fa \
      bfd-first_non_consensus_sequences.fasta \
      rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta \
      nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta \
      rnacentral_active_seq_id_90_cov_80_linclust.fasta)

NOT_COPIED_FILES=()

while (( ${#FILES[@]} )); do
  # Get total size of files to copy in bytes
  SOURCE_FILES=( "${FILES[@]/#/${SOURCE_DIR}/}" )
  TOTAL_SIZE=$(du -sbc "${SOURCE_FILES[@]}" | awk 'END{print $1}')

  # Get available space on target drive in bytes
  AVAILABLE_SPACE=$(df --portability --block-size=1 "$TARGET_DIR" | awk 'END{print $4}')

  # Compare sizes and copy if enough space
  if (( TOTAL_SIZE <= AVAILABLE_SPACE )); then
    printf 'Copying files... %s\n' "${FILES[@]}"
    echo "From ${SOURCE_DIR} -> ${TARGET_DIR}"

    for file in "${FILES[@]}"; do
      cp -r "${SOURCE_DIR}/${file}" "${TARGET_DIR}/" &
    done
    break
  else
    NOT_COPIED_FILES+=("${FILES[-1]}")
    unset 'FILES[-1]'
  fi
done

printf 'No room left on ssd for: %s\n' "${NOT_COPIED_FILES[@]}"
wait
