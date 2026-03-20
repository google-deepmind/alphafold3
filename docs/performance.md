# Performance

## Running the Pipeline in Stages

The `run_alphafold.py` script can be executed in stages to optimise resource
utilisation. This can be useful for:

1.  Splitting the CPU-only data pipeline from model inference (which requires a
    GPU), to optimise cost and resource usage.
1.  Generating the JSON output file from the data pipeline only run and then
    using it for multiple different inference only runs across seeds or across
    variations of other features (e.g. a ligand or a partner chain).
1.  Generating the JSON output for multiple individual monomer chains (e.g. for
    chains A, B, C, D), then running the inference on all possible chain pairs
    (AB, AC, AD, BC, BD, CD) by creating dimer JSONs by merging the monomer
    JSONs. By doing this, the MSA and template search need to be run just 4
    times (once for each chain), instead of 12 times.

### Data Pipeline Only

Launch `run_alphafold.py` with `--run_inference=false` to generate Multiple
Sequence Alignments (MSAs) and templates, without running featurisation and
model inference. This stage can be quite costly in terms of runtime, CPU, and
RAM use. The output will be JSON files augmented with MSAs and templates that
can then be directly used as input for running inference.

### Pre-computing and reusing MSA and templates

When folding multiple candidate chains with a set of fixed chains (i.e. chains
that are the same for all the runs), you can optimize the process by computing
the MSA and templates for the fixed chains only once. The computations for the
changing candidate chains will still be performed for each run:

1.  Run the AlphaFold 3 data pipeline for the fixed chains using the
    `--run_inference=false` flag. This step generates a JSON file containing the
    MSA and template data for these chains.
2.  When constructing your multimer input JSONs, populate the entries for the
    fixed chains using the data generated in the previous step.
    *   For the fixed chains: Specifically, copy the `unpairedMsa`, `pairedMsa`,
        and `templates` fields from the pre-computed JSON into the multimer
        input JSON. This prevents these fields from being recomputed.
    *   For the candidate chains: Leave these fields unset (or `null`) in the
        multimer input JSON. This will signal the pipeline to compute them
        dynamically for each run.

This technique can also be extended to efficiently process all combinations of
*n* first chains and *m* second chains. Instead of performing *n* × *m* full
computations, you can reduce this to *n* + *m* data pipeline runs.

In this scenario:

1.  Run the data pipeline (step 1 above, with `--run_inference=false`) for all
    *n* individual first chains and all *m* individual second chains.
2.  Assemble the dimer input JSONs for each desired pair by combining their
    respective pre-computed monomer JSONs.
3.  Run only the inference step on these assembled JSONs using the
    `--run_data_pipeline=false` flag.

This approach has been discussed in multiple GitHub issues, such as:
https://github.com/google-deepmind/alphafold3/issues/171 (which links to other
similar issues).

### Featurisation and Model Inference Only

Launch `run_alphafold.py` with `--run_data_pipeline=false` to skip the data
pipeline and run only featurisation and model inference. This stage requires the
input JSON file to contain pre-computed MSAs and templates (or they must be
explicitly set to empty if you want to run MSA and template free).

## Data Pipeline

The runtime of the data pipeline (i.e. genetic sequence search and template
search) can vary significantly depending on the size of the input and the number
of homologous sequences found, as well as the available hardware – the disk
speed can influence genetic search speed in particular.

If you would like to improve performance, it's recommended to increase the disk
speed (e.g. by leveraging a RAM-backed filesystem), or increase the available
CPU cores and add more parallelisation. This can help because AlphaFold 3 runs
genetic search against 4 databases in parallel, so the optimal number of cores
is the number of cores used for each Jackhmmer process times 4. Also note that
for sequences with deep MSAs, Jackhmmer or Nhmmer may need a substantial amount
of RAM beyond the recommended 64 GB of RAM.

### Sharded genetic databases

The run time of the genetic database search can be *significantly* sped up by
splitting the genetic databases if a machine with many CPU cores is used and the
databases are on very fast SSD or in a RAM-backed filesystem. With this
technique you can make Jackhmmer/Nhmmer genetic search fully utilize your
hardware and take advantage of multi-core systems.

Each genetic database with *n* sequences is split into *s* shards, each
containing roughly *n* / *s* sequences. We recommend splitting the sequences
between shards randomly to make sure each shard has similar sequence length
distribution. This could be achieved using standard tools:

1.  Shuffle the sequences in the fasta. This can be done for example by running:
    `seqkit shuffle --two-pass <db.fasta>`
2.  Split the shuffled fasta in *s* shards. This can be done for example by
    running: `seqkit split2 --by-part <s> <db.fasta>`

Make sure the shards names follow this pattern:
`prefix-<shard_index>-of-<total_shards>`, both `shard_index` and `total_shards`
having always 5 digits, with leading zeros as needed. The `shard_index` goes
from 0 to `total_shards - 1`. A file "path" (spec) for a sharded file is
`prefix@<total_shards>`.

E.g. for a file named `uniprot.fasta` split into 3 shards, the names of the
shards should be:

*   `uniprot.fasta-00000-of-00003`
*   `uniprot.fasta-00001-of-00003`
*   `uniprot.fasta-00002-of-00003`

The file spec for these files is `uniprot.fasta@3`.

Save the total number of sequences in the protein databases, and the total
number of nucleic bases in the RNA databases – these will be needed later as a
flag to Jackhmmer/Nhmmer to correctly scale e-values across all shards.

Save the sharded databases on a fast SSD or in a RAM-backed filesystem, then
launch AlphaFold with the sharded paths instead of normal paths and set the
Z-values.

For instance with each database sharded into 16 shards:

```bash
python run_alphafold.py \
    --small_bfd_database_path="bfd-first_non_consensus_sequences.fasta@64" \
    --small_bfd_z_value=65984053 \
    --mgnify_database_path="mgy_clusters_2022_05.fa@512" \
    --mgnify_z_value=623796864 \
    --uniprot_cluster_annot_database_path="uniprot_cluster_annot_2021_04.fasta@256" \
    --uniprot_cluster_annot_z_value=225619586 \
    --uniref90_database_path="uniref90_2022_05.fasta@128" \
    --uniref90_z_value=153742194 \
    --ntrna_database_path="nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta@256" \
    --ntrna_z_value=76752.808514 \
    --rfam_database_path="rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta@16" \
    --rfam_z_value=138.115553 \
    --rna_central_database_path="rnacentral_active_seq_id_90_cov_80_linclust.fasta@64" \
    --rna_central_z_value=13271.415730
    --jackhmmer_n_cpu=2 \
    --jackhmmer_max_parallel_shards=16 \
    --nhmmer_n_cpu=2 \
    --nhmmer_max_parallel_shards=16
```

This run will utilize (2 CPUs) × (16 max parallel shards) × (4 protein dbs
searched in parallel) = 128 cores for each protein chain, and (2 CPUs) × (16 max
parallel shards) × (3 RNA dbs searched in parallel) = 96 cores for each RNA
chain. Make sure to tune:

*   the Jackhmmer/Nhmmer number of CPUs,
*   the maximum number of shards searched in parallel,
*   and the number of shards for each database

so that the memory bandwidth and CPUs on your machine are optimally utilized.
You should aim for consistent shard sizes across all databases (so e.g. if
database A is split into 16 shards and is 3× smaller than database B, database B
should be split into 3 × 16 = 48 shards).

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
          (2048, None),
          (3072, 1024),
          (None, 512),
      )
    ```

The format of entries in `pair_transition_shard_spec` is
`(num_tokens_upper_bound, shard_size)`. Setting `shard_size=None` means there is
no upper bound.

For the example above:

*   `(2048, None)`: for sequences up to 2,048 tokens, do not shard
*   `(3072, 1024)`: for sequences up to 3,072 tokens, shard in chunks of 1,024
*   `(None, 512)`: for all other sequences, shard in chunks of 512

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
