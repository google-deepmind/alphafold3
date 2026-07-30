# Frequently Asked Questions

## Will two runs with the same input and settings give the same result?

Yes for the algorithm, no for the bytes. Every stochastic step in AlphaFold 3 is
derived from the seed in the input JSON, so the same input, the same flags, the
same machine and the same databases produce the same structure. However, nothing
in the code guarantees bit-identical numerics, and the output files are never
byte-identical because a timestamp is written into the mmCIF.

The details are broken down below by layer, since each layer can be pinned or
broken independently.

### 1. All randomness derives from `modelSeeds`

Seeds are read from the input JSON
([`folding_input.py`](../src/alphafold3/common/folding_input.py)):

```python
rng_seeds=[int(seed) for seed in raw_json['modelSeeds']],
```

Featurisation creates one `RandomState` per seed
([`featurisation.py`](../src/alphafold3/data/featurisation.py)):

```python
batch = data_pipeline.process_item(
    fold_input=fold_input,
    ccd=ccd,
    random_state=np.random.RandomState(rng_seed),
    random_seed=rng_seed,
)
```

The featurisation pipeline then re-seeds from the integer rather than continuing
to consume the passed state, so the result does not depend on how many draws
happened earlier
([`model/pipeline/pipeline.py`](../src/alphafold3/model/pipeline/pipeline.py)):

```python
if random_seed is None:
  random_seed = random_state.randint(2**31)
random_state = np.random.RandomState(seed=random_seed)
```

That chain reaches the genuinely stochastic parts of featurisation. RDKit ligand
conformer generation (a stochastic ETKDG embedding) takes its seed from it
([`model/features.py`](../src/alphafold3/model/features.py),
[`rdkit_utils.py`](../src/alphafold3/data/tools/rdkit_utils.py)):

```python
conformer_random_seed = int(random_state.randint(1, 1 << 31))
...
params = rd_all_chem.ETKDGv3()
params.randomSeed = random_seed
```

Reference frame construction is pinned to a hard-coded constant:

```python
_DETERMINISTIC_FRAMES_RANDOM_SEED = 12312837
```

Inference gets one PRNG key per seed
([`run_alphafold.py`](../run_alphafold.py)):

```python
rng_key = jax.random.PRNGKey(seed)
result = model_runner.run_inference(example, rng_key)
```

Everything downstream splits off that key: diffusion noise, `random_augmentation`
in `diffusion_head.py`, and MSA row sampling in `featurization.py`. There is no
unseeded RNG in the model path.

`--num_seeds` is deterministic as well; the extra seeds are consecutive
integers, not samples:

```python
rng_seeds=list(range(self.rng_seeds[0], self.rng_seeds[0] + num_seeds)),
```

### 2. When seeds are sampled for you, results legitimately differ

With the `alphafoldserver` JSON format and `"modelSeeds": []`, a fresh seed is
drawn on every run:

```python
def _sample_rng_seed() -> int:
  """Sample a random seed for AlphaFoldServer job."""
  return random.randint(0, 2**32 - 1)
```

The same applies to inputs converted from mmCIF. In these cases the two runs are
not actually the same job. The seeds that were used are recorded in the
`*_data.json` file AlphaFold 3 writes to the output directory.

### 3. The data pipeline is deterministic relative to fixed databases

The MSA path is written to be order-stable:

*   Deduplication preserves input order (`msa.py`: "deduplicated in the input
    order").
*   MSAs are combined in a fixed list order (`data/pipeline.py`).
*   The sharded search merge sorts on a total order with a name tiebreak, so
    thread completion order cannot leak into the result
    ([`jackhmmer.py`](../src/alphafold3/data/tools/jackhmmer.py)):

    ```python
    e_value, bit_score = tbl_info.split(maxsplit=6)[4:6]
    return float(e_value), -float(bit_score), name
    ```

Template filtering uses a fixed date default rather than the wall clock:

```python
_MAX_TEMPLATE_DATE = flags.DEFINE_string(
    'max_template_date',
    '2021-09-30',  # By default, use the date from the AlphaFold 3 paper.
```

Two things break this in practice:

*   If the sequence databases or the PDB structure store changed between runs,
    the MSA and templates change.
*   Sharded and unsharded genetic search do not produce the same MSA. See
    [Known Issues](known_issues.md#msa-discrepancy-between-alphafold-3-and-alphafold-server).
    Same flags, different database layout, different answer.

Note also that in the unsharded path the hit ordering produced by
Jackhmmer/Nhmmer is consumed as-is; that ordering is outside this repository's
control.

To remove this layer entirely, run the data pipeline once with
`--norun_inference` and feed the resulting JSON to both inference runs. This is
the split workflow described in the
[performance documentation](performance.md).

### 4. Numerical bit-identity is not guaranteed

The forward pass is a jitted, device-pinned function:

```python
return functools.partial(
    jax.jit(forward_fn.apply, device=self._device), self.model_params
)
```

Given the same compiled executable on the same GPU, this is deterministic. It
stops being deterministic when the executable differs:

*   a different `--flash_attention_implementation` (`triton`, `cudnn`, `xla`),
*   a different GPU architecture,
*   different `XLA_FLAGS`,
*   a different JAX or tokamax version.

The model runs in bfloat16 by default (`model_config.py`, `bfloat16 = 'all'`),
so small reduction-order differences are amplified through 10 recycles and the
diffusion sampler. [Known Issues](known_issues.md) shows both ends of this
sensitivity: CUDA Capability 7.x GPUs produce unusable output without a specific
XLA flag, and tokamax silently substitutes a different kernel when one is not
available on the device.

One caveat that cannot be verified from this repository: XLA's GPU autotuner
selects kernel configurations by timing them, which can in principle vary
between runs on a loaded machine. Nothing here pins it. Treat repeated-run
bit-identity as likely rather than guaranteed.

If numerics do shift slightly, the ranking can visibly flip, because best-sample
selection is a strict comparison over the sample order:

```python
if max_ranking_score is None or ranking_score > max_ranking_score:
```

### 5. The output files are never byte-identical

Wall-clock time is written into the mmCIF metadata
([`post_processing.py`](../src/alphafold3/model/post_processing.py)):

```python
timestamp = datetime.datetime.now().isoformat(sep=' ', timespec='seconds')
```

The second run also will not write to the same directory: a non-empty output
directory gets a timestamped sibling unless `--force_output_dir` is passed.

```python
new_output_dir = (
    output_dir.parent
    / f'{output_dir.name}_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}'
)
```

So `diff` on the mmCIF will always report a difference. Compare coordinates and
the confidence JSON files instead.

### Summary

Reproducible when: seeds are explicit in the input JSON, the same
container/binary is used, the same GPU model, the same databases (including the
same sharding), and the same flags.

Not reproducible when: seeds are auto-sampled, the databases were updated, or
the GPU, attention implementation, or XLA flags changed.
