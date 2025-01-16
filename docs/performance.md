# Performance

## Data Pipeline

The runtime of the data pipeline (i.e. genetic sequence search and template
search) can vary significantly depending on the size of the input and the number
of homologous sequences found, as well as the available hardware (disk speed can
influence genetic search speed in particular). If you would like to improve
performance, it’s recommended to increase the disk speed (e.g. by leveraging a
RAM-backed filesystem), or increase the available CPU cores and add more
parallelisation. Also note that for sequences with deep MSAs, Jackhmmer or
Nhmmer may need a substantial amount of RAM beyond the recommended 64 GB of RAM.

## Model Inference

Table 8 in the Supplementary Information of the
[AlphaFold 3 paper](https://nature.com/articles/s41586-024-07487-w) provides
compile-free inference timings for AlphaFold 3 when configured to run on 16
NVIDIA A100s, with 40 GB of memory per device. In contrast, this repository
supports running AlphaFold 3 on a single NVIDIA A100 with 80 GB of memory in a
configuration optimised to maximise throughput.

We compare compile-free inference timings of these two setups in the table below
using GPU seconds (i.e. multiplying by 16 when using 16 A100s). The setup in
this repository is more efficient (by at least 2×) across all token sizes,
indicating its suitability for high-throughput applications.

Num Tokens | 1 A100 80 GB (GPU secs) | 16 A100 40 GB (GPU secs) | Improvement
:--------- | ----------------------: | -----------------------: | ----------:
1024       | 62                      | 352                      | 5.7×
2048       | 275                     | 1136                     | 4.1×
3072       | 703                     | 2016                     | 2.9×
4096       | 1434                    | 3648                     | 2.5×
5120       | 2547                    | 5552                     | 2.2×

## Running the Pipeline in Stages

The `run_alphafold.py` script can be executed in stages to optimise resource
utilisation. This can be useful for:

1.  Splitting the CPU-only data pipeline from model inference (which requires a
    GPU), to optimise cost and resource usage.
1.  Caching the results of MSA/template search, then reusing the augmented JSON
    for multiple different inferences across seeds or across variations of other
    features (e.g. a ligand).

### Data Pipeline Only

Launch `run_alphafold.py` with `--norun_inference` to generate Multiple Sequence
Alignments (MSAs) and templates, without running featurisation and model
inference. This stage can be quite costly in terms of runtime, CPU, and RAM use.
The output will be JSON files augmented with MSAs and templates that can then be
directly used as input for running inference.

### Featurisation and Model Inference Only

Launch `run_alphafold.py` with `--norun_data_pipeline` to skip the data pipeline
and run only featurisation and model inference. This stage requires the input
JSON file to contain pre-computed MSAs and templates.

## Accelerator Hardware Requirements

We officially support the following configurations, and have extensively tested
them for numerical accuracy and throughput efficiency:

-   1 NVIDIA A100 (80 GB)
-   1 NVIDIA H100 (80 GB)

We compare compile-free inference timings of both configurations in the
following table:

Num Tokens | 1 A100 80 GB (seconds) | 1 H100 80 GB (seconds)
:--------- | ---------------------: | ---------------------:
1024       | 62                     | 34
2048       | 275                    | 144
3072       | 703                    | 367
4096       | 1434                   | 774
5120       | 2547                   | 1416

### Other Hardware Configurations

#### NVIDIA A100 (40 GB)

AlphaFold 3 can run on inputs of size up to 4,352 tokens on a single NVIDIA A100
(40 GB) with the following configuration changes:

1.  Enabling [unified memory](#unified-memory).
1.  Adjusting `pair_transition_shard_spec` in `model_config.py`:

    ```py
      pair_transition_shard_spec: Sequence[_Shape2DType] = (
          (2048, None),  # for sequences up to 2048 tokens, do not shard
          (3072, 1024),  # for sequences up to 3072 tokens, shard in chunks of 1024
          (None, 512),   # for all other sequences, shard in chunks 512
      )
    ```

While numerically accurate, this configuration will have lower throughput
compared to the set up on the NVIDIA A100 (80 GB), due to less available memory.

#### NVIDIA V100

There are known numerical issues with CUDA Capability 7.x devices. To work
around the issue, set the ENV XLA_FLAGS to include
`--xla_disable_hlo_passes=custom-kernel-fusion-rewriter`.

With the above flag set, AlphaFold 3 can run on inputs of size up to 1,280
tokens on a single NVIDIA V100 using [unified memory](#unified-memory).

#### NVIDIA P100

AlphaFold 3 can run on inputs of size up to 1,024 tokens on a single NVIDIA P100
with no configuration changes needed.

#### Other devices

Large-scale numerical tests have not been performed on any other devices but
they are believed to be numerically accurate.

There are known numerical issues with CUDA Capability 7.x devices. To work
around the issue, set the environment variable `XLA_FLAGS` to include
`--xla_disable_hlo_passes=custom-kernel-fusion-rewriter`.

## Compilation Buckets

To avoid excessive re-compilation of the model, AlphaFold 3 implements
compilation buckets: ranges of input sizes using a single compilation of the
model.

When featurising an input, AlphaFold 3 determines the smallest bucket the input
fits into, then adds any necessary padding. This may avoid re-compiling the
model when running inference on the input if it belongs to the same bucket as a
previously processed input.

The configuration of bucket sizes involves a trade-off: more buckets leads to
more re-compilations of the model, but less padding.

By default, the largest bucket size is 5,120 tokens. Processing inputs larger
than this maximum bucket size triggers the creation of a new bucket for exactly
that input size, and a re-compilation of the model. In this case, you may wish
to redefine the compilation bucket sizes via the `--buckets` flag in
`run_alphafold.py` to add additional larger bucket sizes. For example, suppose
you are running inference on inputs with token sizes: `5132, 5280, 5342`. Using
the default bucket sizes configured in `run_alphafold.py` will trigger three
separate model compilations, one for each unique token size. If instead you pass
in the following flag to `run_alphafold.py`

```
--buckets 256,512,768,1024,1280,1536,2048,2560,3072,3584,4096,4608,5120,5376
```

when running inference on the above three input sizes, the model will be
compiled only once for the bucket size `5376`. **Note:** for this specific
example with input sizes `5132, 5280, 5342`, passing in `--buckets 5376` is
sufficient to achieve the desired compilation behaviour. The provided example
with multiple buckets illustrates a more general solution suitable for diverse
input sizes.

## Additional Flags

### Compilation Time Workaround with XLA Flags

To work around a known XLA issue causing the compilation time to greatly
increase, the following environment variable must be set (it is set by default
in the provided `Dockerfile`).

```sh
ENV XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
```

### CUDA Capability 7.x GPUs

For all CUDA Capability 7.x GPUs (e.g. V100) the environment variable
`XLA_FLAGS` must be changed to include
`--xla_disable_hlo_passes=custom-kernel-fusion-rewriter`. Disabling the Tritron
GEMM kernels is not necessary as they are not supported for such GPUs.

```sh
ENV XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
```

### GPU Memory

The following environment variables (set by default in the `Dockerfile`) enable
folding a single input of size up to 5,120 tokens on a single A100 (80 GB) or a
single H100 (80 GB):

```sh
ENV XLA_PYTHON_CLIENT_PREALLOCATE=true
ENV XLA_CLIENT_MEM_FRACTION=0.95
```

#### Unified Memory

If you would like to run AlphaFold 3 on inputs larger than 5,120 tokens, or on a
GPU with less memory (an A100 with 40 GB of memory, for instance), we recommend
enabling unified memory. Enabling unified memory allows the program to spill GPU
memory to host memory if there isn't enough space. This prevents an OOM, at the
cost of making the program slower by accessing host memory instead of device
memory. To learn more, check out the
[NVIDIA blog post](https://developer.nvidia.com/blog/unified-memory-cuda-beginners/).

You can enable unified memory by setting the following environment variables in
your `Dockerfile`:

```sh
ENV XLA_PYTHON_CLIENT_PREALLOCATE=false
ENV TF_FORCE_UNIFIED_MEMORY=true
ENV XLA_CLIENT_MEM_FRACTION=3.2
```

### JAX Persistent Compilation Cache

You may also want to make use of the JAX persistent compilation cache, to avoid
unnecessary recompilation of the model between runs. You can enable the
compilation cache with the `--jax_compilation_cache_dir <YOUR_DIRECTORY>` flag
in `run_alphafold.py`.

More detailed instructions are available in the
[JAX documentation](https://jax.readthedocs.io/en/latest/persistent_compilation_cache.html#persistent-compilation-cache),
and more specifically the instructions for use on
[Google Cloud](https://jax.readthedocs.io/en/latest/persistent_compilation_cache.html#persistent-compilation-cache).
In particular, note that if you would like to make use of a non-local
filesystem, such as Google Cloud Storage, you will need to install
[`etils`](https://github.com/google/etils) (this is not included by default in
the AlphaFold 3 Docker container).
