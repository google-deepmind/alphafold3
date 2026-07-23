# Running AlphaFold 3 on Apple Silicon GPU (Metal / MPS)

> **Status: experimental feasibility path.** This mode runs AlphaFold 3
> inference natively on the Apple Silicon integrated GPU through the community
> [`jax-mps`](https://pypi.org/project/jax-mps/) Metal plugin for JAX. It is
> **not** an officially supported or numerically certified configuration. The
> NVIDIA CUDA path is unchanged and remains the reference. Use this for research,
> local development, and prototyping on Macs.

This guide describes a **native build with no Docker**. Docker Desktop on macOS
runs a Linux VM that cannot see the Mac GPU, so the containerized recipe in
[`installation.md`](installation.md) can only ever reach the CPU. To use the
Apple GPU you must build AlphaFold 3 natively, as described below.

If you only want CPU-only inference on a Mac (much slower), follow the "Running
AlphaFold 3 without a GPU" section of [`installation.md`](installation.md)
instead.

## What works and what to expect

Validated in July 2026 across a 52-target campaign (all 12 molecule-type classes
plus token-scaling targets up to ~2,200 tokens) on an **M2 Ultra Mac Studio
(76-GPU-core, 192 GB unified memory)**, macOS 15.

- **Correctness:** the full data pipeline (HMMER MSA + templates) and inference
  run end to end; predictions match the CUDA reference closely for the cases
  tested, and inference is **bit-reproducible run-to-run** with fixed seeds.
- **Attention:** only the portable `xla` implementation works on Metal;
  `triton` and `cudnn` require NVIDIA GPUs and are rejected early.
- **Memory is the main limit:** peak GPU (unified) memory reaches ~167 GB at
  ~2,200 tokens on a 192 GB machine — the practical single-structure ceiling on
  this hardware (see [Performance and memory](#performance-and-memory)).
- **Large targets are slow:** per-seed inference scales roughly cubically with
  tokens, with a steep Metal-specific penalty above ~1,500 tokens.

## Prerequisites

- Apple Silicon Mac (M1/M2/M3/M4 family). More unified memory = larger tractable
  structures; 64 GB is a reasonable floor, 128–192 GB recommended for large
  complexes.
- macOS 14 or newer.
- Xcode Command Line Tools (provides `clang`, `make`, `patch`):

  ```sh
  xcode-select --install
  ```

- [`uv`](https://docs.astral.sh/uv/) and a Python 3.11–3.13 toolchain (a
  Conda/Miniforge environment works too — see
  [Alternative: Conda](#alternative-conda-environment)).

## Step 1 — Install the HMMER suite (for the MSA data pipeline)

HMMER is only needed if you run the data pipeline (i.e. not with
`--norun_data_pipeline`). Prefer a user-local source build of HMMER 3.4 with the
sequence-truncation patch that ships in this repository
(`docker/jackhmmer_seq_limit.patch`):

```sh
export AF3_REPO="$PWD"                    # path to your alphafold3 clone
export HMMER_PREFIX="$HOME/af3-tools/hmmer-3.4"
mkdir -p "$HOME/af3-tools/hmmer-source" && cd "$HOME/af3-tools/hmmer-source"

curl -LO http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz
echo "ca70d94fd0cf271bd7063423aabb116d42de533117343a9b27a65c17ff06fbf3  hmmer-3.4.tar.gz" | shasum -a 256 -c
tar -xzf hmmer-3.4.tar.gz
cd hmmer-3.4
patch -p0 < "$AF3_REPO/docker/jackhmmer_seq_limit.patch"
./configure --prefix="$HMMER_PREFIX"
make -j"$(sysctl -n hw.perflevel0.logicalcpu)"
make install
"$HMMER_PREFIX/bin/jackhmmer" -h | head -1   # sanity check
```

Pass the binaries to `run_alphafold.py` in Step 4 via the
`--jackhmmer_binary_path`, `--nhmmer_binary_path`, `--hmmalign_binary_path`,
`--hmmsearch_binary_path`, and `--hmmbuild_binary_path` flags (or put
`$HMMER_PREFIX/bin` first on your `PATH`).

## Step 2 — Get the code, model weights, and databases

1. Clone the repository (if you have not already):

   ```sh
   git clone https://github.com/google-deepmind/alphafold3.git
   cd alphafold3
   ```

2. Obtain the model weights by following the process in the
   [main README](https://github.com/google-deepmind/alphafold3). Keep the
   weights **outside** the git checkout.
3. Download the genetic databases (see [`installation.md`](installation.md) and
   `fetch_databases.py`) onto fast local storage.

## Step 3 — Build and install AlphaFold 3 (with the MPS backend)

On `darwin arm64`, `pyproject.toml` adds the `jax-mps` Metal backend (the base
install already provides `jax`/`jaxlib` 0.10.2; the CUDA plugin is excluded on
macOS), so a normal install pulls in everything needed for GPU inference:

```sh
uv venv --python 3.12
source .venv/bin/activate
uv sync            # builds/installs alphafold3 AND the jax-mps Metal backend
uv run build_data  # builds the CCD / chemical-component data used at runtime
```

Verify that JAX sees the Metal device:

```sh
python -c "import jax; print([d.platform for d in jax.local_devices()])"
# -> ['mps']   (a one-time 'Platform mps is experimental' warning is expected)
```

## Step 4 — Run inference on the Apple GPU

Select the Metal backend with `--jax_backend mps` and the portable `xla`
attention implementation:

```sh
uv run python run_alphafold.py \
  --json_path=/path/to/fold_input.json \
  --model_dir=/path/to/weights \
  --db_dir=/path/to/databases \
  --output_dir=/path/to/output \
  --jax_backend=mps \
  --gpu_device=0 \
  --flash_attention_implementation=xla \
  --jackhmmer_binary_path="$HMMER_PREFIX/bin/jackhmmer" \
  --nhmmer_binary_path="$HMMER_PREFIX/bin/nhmmer" \
  --hmmalign_binary_path="$HMMER_PREFIX/bin/hmmalign" \
  --hmmsearch_binary_path="$HMMER_PREFIX/bin/hmmsearch" \
  --hmmbuild_binary_path="$HMMER_PREFIX/bin/hmmbuild"
```

For inference on a pre-computed input (skipping the MSA pipeline and HMMER), add
`--norun_data_pipeline` and drop the HMMER flags.

## Performance and memory

Measured on the M2 Ultra (192 GB). Inference time and peak memory are set by the
**padded token bucket**; MSA (CPU) is driven by chain count / RNA, not token
count. Treat as order-of-magnitude guidance; other Apple Silicon chips differ.

| Padded tokens | MSA time¹ | Inference / seed | Peak unified memory (GB) |
|--:|--:|--:|--:|
| 256 | ~7 min | ~29 s | 36 |
| 512 | ~13 min | ~90 s | 50 |
| 768 | ~13 min | ~3.9 min | 73 |
| 1,024 | ~7 min | ~8.6 min | 101 |
| 1,536 | ~19 min | ~26.6 min | 155 |
| ~2,180 | ~20 min | ~3.9 h | 167 |

¹ MSA / data pipeline is CPU-bound and driven by the number of unique chains and
RNA presence, not token count (observed range 7–35 min).

Key points:

- **Unified memory is the ceiling.** GPU tensors live in unified memory (macOS
  "Memory Used") and are **not** reflected in process RSS (~7 GB). On a 192 GB
  machine, ~2,200 tokens is the practical single-structure limit; larger
  structures will run out of memory. Scale expectations to your Mac's memory.
- **No persistent JAX compilation cache on Metal.** `--jax_compilation_cache_dir`
  does not persist across processes on MPS, so fold many inputs in one process
  (e.g. `--input_dir`) to amortize compile cost.
- Set HMMER thread counts so that (concurrent searches × threads-per-search) ≈
  the number of performance cores (`sysctl hw.perflevel0.logicalcpu`).

## Alternative: Conda environment

`pyproject.toml` still selects `jax-mps` on `darwin arm64`, so a plain install
pulls in the Metal backend:

```sh
conda create -n af3-mps python=3.13 cmake -y
conda activate af3-mps
pip install --no-build-isolation -e .   # installs alphafold3 + jax-mps
build_data
```

Then run `run_alphafold.py` with the same MPS flags as in Step 4.

## Troubleshooting

- **`No local MPS device was found`** — `jax-mps` is not installed. Re-check the
  install (`python -c "import jax; print(jax.local_devices())"`).
- **`--flash_attention_implementation must be set to "xla"`** — `triton`/`cudnn`
  are NVIDIA-only; use `xla` on Metal.
- **Out-of-memory / heavy swapping on large targets** — reduce the structure
  size or run on a Mac with more unified memory (see the table above).
