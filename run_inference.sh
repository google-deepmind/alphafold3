#!/usr/bin/bash
#

set -e

CWD=`realpath -s $0`
CWD=`dirname ${CWD}`

cd ${CWD}

help() {
  echo -e "usage: `basename $0` [options] -- [af3_opts] input_json1 input_json2 ...\n" \
          "options:\n" \
          "  -t, --timestamp <timestamp>\n" \
          "  -f, --force\n" \
    "\n" \
    "platforms:\n" \
    "  --msa_tool_home <msa_tool_home>\n" \
    "  --msa_tool_srv <msa_tool_srv>\n"
  exit $1
}

ARGS=$(getopt --options "t:fh" --longoptions "timestamp:,force_output_dir,help" -- "$@") || exit
eval "set -- ${ARGS}"
echo "ARGS: $*"
while true; do
  case "$1" in
    (-t | --timestamp) timestamp="$2"; shift 2;;
    (--force_output_dir) force_output_dir=1; shift 2;;
    (-h | --help) help 0 ;;
    (--) shift 1; break;;
    (*) help 1;
  esac
done

# To work around a known XLA issue causing the compilation time to greatly
# increase, the following environment variable setting XLA flags must be enabled
# when running AlphaFold 3. Note that if using CUDA capability 7 GPUs, it is
# necessary to set the following XLA_FLAGS value instead:
# export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
# (no need to disable gemm in that case as it is not supported for such GPU).
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
# Memory settings used for folding up to 5,120 tokens on A100 80 GB.
export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95

uv run python run_alphafold.py \
  --norun_data_pipeline \
  --force_output_dir \
  --jax_compilation_cache_dir ${CWD}/.cache/jax_compilation \
  --model_dir ${HOME}/public_databases/models \
  --input_dir ${af3_pred_data_dir}/msa \
  --output_dir ${af3_pred_data_dir}/msa
