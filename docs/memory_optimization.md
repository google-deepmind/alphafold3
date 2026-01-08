# Memory Optimization Guide

This guide explains how to handle Out of Memory (OOM) errors when running AlphaFold 3, particularly for large homomer complexes.

## Understanding Memory Requirements

AlphaFold 3's memory usage scales with:

1. **Number of tokens** (O(n²) for attention operations)
2. **Number of chains** (homomers require more memory)
3. **Number of recycles** (default: 10)
4. **Number of diffusion samples** (default: 5)
5. **Flash attention implementation** (significant impact)

### Memory Scaling for Homomers

Homomer complexes (multiple copies of the same protein) have significantly higher memory requirements:

| Tokens | Single Chain | 2-mer | 4-mer | 8-mer |
|--------|--------------|-------|-------|-------|
| 1024   | ~15 GB       | ~22 GB| ~35 GB| ~65 GB|
| 2048   | ~28 GB       | ~42 GB| ~68 GB| OOM   |
| 4608   | ~58 GB       | ~78 GB| ~95 GB| OOM   |

## Automatic Memory Optimization

AlphaFold 3 now includes automatic memory optimization (enabled by default):

```bash
python run_alphafold.py \
    --json_path=input.json \
    --model_dir=/root/models \
    --output_dir=/root/output \
    --auto_memory_optimization=true  # Default
```

When enabled, the system will:
1. Estimate memory requirements before inference
2. Detect if OOM is likely
3. Automatically reduce `num_recycles` and `num_diffusion_samples`
4. Print what optimizations were applied

### Example Output

```
Memory Analysis for WP_147157570.1_copies_4:
  Estimated tokens (with padding): 4608
  Number of chains: 4
  Available GPU memory: 80.0 GB

Memory Estimate:
  Total: 86.42 GB
  Model parameters: 5.00 GB
  Embeddings: 5.38 GB
  Attention: 65.04 GB
  Diffusion: 11.00 GB
  Tokens: 4608
  Chains: 4
  Homomer: True

Warning: Estimated memory (86.4 GB) exceeds available memory (80.0 GB)

Optimization Suggestion:
  num_recycles: 5
  num_diffusion_samples: 3
  flash_attention_implementation: triton
  Estimated memory: 72.15 GB
  Will fit: Yes

Applying automatic memory optimizations:
  num_recycles: 10 -> 5
  num_diffusion_samples: 5 -> 3
```

## Manual Optimization

### 1. Reduce Recycles and Samples

```bash
python run_alphafold.py \
    --json_path=input.json \
    --model_dir=/root/models \
    --output_dir=/root/output \
    --num_recycles=5 \
    --num_diffusion_samples=3
```

**Impact:**
- Reduces memory by ~30-40%
- Minimal impact on prediction quality
- Faster inference time

### 2. Enable Flash Attention

```bash
python run_alphafold.py \
    --json_path=input.json \
    --model_dir=/root/models \
    --output_dir=/root/output \
    --flash_attention_implementation=triton  # or cudnn
```

**Requirements:**
- NVIDIA Ampere GPUs or later (A100, H100, RTX 3090+)
- Reduces attention memory by ~70%

For older GPUs (V100, P100):
```bash
--flash_attention_implementation=xla
```

### 3. Memory Estimation Only

Estimate memory without running inference:

```bash
python run_alphafold.py \
    --json_path=input.json \
    --model_dir=/root/models \
    --output_dir=/root/output \
    --estimate_memory_only=true
```

This will print memory estimates and suggestions, then exit.

### 4. Override GPU Memory Limit

If you want to enable CPU memory spillover:

```bash
python run_alphafold.py \
    --json_path=input.json \
    --model_dir=/root/models \
    --output_dir=/root/output \
    --max_gpu_memory_gb=160  # Allow spillover to system RAM
```

**Warning:** This will be significantly slower due to PCIe transfer overhead.

## Docker Configuration

### Recommended Environment Variables

```bash
docker run --rm \
    --volume ./input:/root/af_input \
    --volume ./output:/root/af_output \
    --volume $HOME/AF3_model:/root/models \
    --volume $HOME/AF3_db:/root/public_databases \
    --gpus all \
    -e XLA_PYTHON_CLIENT_PREALLOCATE=false \
    -e XLA_CLIENT_MEM_FRACTION=0.90 \
    -e XLA_FLAGS="--xla_gpu_enable_async_all_reduce=true" \
    alphafold3 python run_alphafold.py \
        --input_dir=/root/af_input \
        --model_dir=/root/models \
        --output_dir=/root/af_output \
        --flash_attention_implementation=triton
```

### Common Mistakes to Avoid

❌ **Don't do this:**
```bash
-e XLA_PYTHON_CLIENT_PREALLOCATE=true  # Preallocates all memory
-e XLA_CLIENT_MEM_FRACTION=0.98        # Too aggressive
-e XLA_CLIENT_MEM_FRACTION=3.2         # Invalid (>1.0)
```

✅ **Do this instead:**
```bash
-e XLA_PYTHON_CLIENT_PREALLOCATE=false
-e XLA_CLIENT_MEM_FRACTION=0.90
```

## Strategies for Large Homomers

### Strategy 1: Fold Single Copy First

For a 4-copy homomer:

1. Create input JSON with single copy
2. Fold the single copy
3. Use the result as a template for the full complex

### Strategy 2: Fold in Stages

For an 8-copy homomer:

1. Fold 2-copy version
2. Use result as template for 4-copy
3. Use 4-copy result as template for 8-copy

### Strategy 3: Reduce Bucket Size

If your input has 4184 tokens but gets padded to 4608:

```bash
--buckets=256,512,768,1024,1280,1536,2048,2560,3072,3584,4096
# Removed 4608 and 5120
```

This forces padding to 4096 instead, saving memory.

## Troubleshooting

### Error: "Out of memory trying to allocate X bytes"

**Solution:**
1. Enable auto-optimization (should be on by default)
2. Manually reduce `--num_recycles=3 --num_diffusion_samples=1`
3. Enable flash attention
4. Check for background processes using GPU memory

### Error: "Cannot fit even with aggressive optimization"

**Solutions:**
1. Use a GPU with more memory (H100 80GB, A100 80GB)
2. Split homomer into smaller chunks
3. Enable CPU memory spillover with `--max_gpu_memory_gb=160`
4. Reduce bucket size to avoid padding

### Slow Performance with Memory Spillover

If using `XLA_CLIENT_MEM_FRACTION=3.2` or similar:

**Problem:** GPU is swapping data to/from system RAM over PCIe

**Solution:**
- Reduce memory usage instead of relying on spillover
- Use automatic optimization
- Consider smaller bucket sizes

### GPU Memory Already in Use

Check for background processes:

```bash
nvidia-smi
# Look for other processes using GPU memory

# Kill specific process
kill -9 <PID>

# Or reset GPU
sudo nvidia-smi --gpu-reset
```

## Performance vs. Quality Trade-offs

| Setting | Memory Savings | Quality Impact | Speed Impact |
|---------|----------------|----------------|--------------|
| `num_recycles: 10 -> 5` | ~30% | Minimal | 2x faster |
| `num_diffusion_samples: 5 -> 3` | ~15% | Minimal | 1.7x faster |
| `flash_attention: xla -> triton` | ~70% | None | 1.5x faster |
| `num_recycles: 10 -> 3` | ~50% | Moderate | 3x faster |
| `num_diffusion_samples: 5 -> 1` | ~35% | Moderate | 5x faster |

## Hardware Recommendations

### For Large Homomers (>4000 tokens)

**Minimum:**
- NVIDIA H100 80GB or A100 80GB
- Flash attention support (Ampere or later)

**Recommended:**
- NVIDIA H100 80GB with NVLink
- 128+ GB system RAM for data pipeline

### For Medium Complexes (2000-4000 tokens)

**Minimum:**
- NVIDIA A100 40GB or RTX A6000
- Flash attention support

**Recommended:**
- NVIDIA A100 80GB or H100

### For Small Complexes (<2000 tokens)

**Minimum:**
- NVIDIA V100 32GB or RTX 3090 24GB

**Recommended:**
- NVIDIA A100 40GB

## Additional Resources

- [AlphaFold 3 Paper](https://www.nature.com/articles/s41586-024-07487-w)
- [Performance Documentation](performance.md)
- [Input Format Guide](input.md)
- [GitHub Issues](https://github.com/google-deepmind/alphafold3/issues)

## Getting Help

If you continue to experience OOM errors:

1. Run with `--estimate_memory_only=true` and share the output
2. Include your GPU model and available memory
3. Share the input JSON (number of chains, tokens)
4. Report on GitHub Issues with the "OOM" label
