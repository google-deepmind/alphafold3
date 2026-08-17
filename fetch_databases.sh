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

readonly db_dir=${1:-$HOME/public_databases}

for cmd in wget tar zstd ; do
  if ! command -v "${cmd}" > /dev/null 2>&1; then
    echo "${cmd} is not installed. Please install it."
    exit 1
  fi
done

echo "Fetching databases to ${db_dir}"
mkdir -p "${db_dir}"

readonly SOURCE=https://storage.googleapis.com/alphafold-databases/v3.0

echo "Start Fetching and Untarring 'pdb_2022_09_28_mmcif_files.tar'"
wget --quiet --output-document=- \
    "${SOURCE}/pdb_2022_09_28_mmcif_files.tar.zst" | \
    tar --no-same-owner --no-same-permissions \
    --use-compress-program=zstd -xf - --directory="${db_dir}" &

for NAME in mgy_clusters_2022_05.fa \
            bfd-first_non_consensus_sequences.fasta \
            uniref90_2022_05.fa \
            uniprot_all_2021_04.fa \
            pdb_seqres_2022_09_28.fasta \
            rnacentral_active_seq_id_90_cov_80_linclust.fasta \
            nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta \
            rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta ; do
  echo "Start Fetching '${NAME}'"
  wget --quiet --output-document=- "${SOURCE}/${NAME}.zst" | \
      zstd --decompress > "${db_dir}/${NAME}" &
done

wait
echo "Complete"
