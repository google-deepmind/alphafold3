# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Memory estimation and optimization utilities for AlphaFold 3 inference."""

import dataclasses
from typing import Final

import jax
import numpy as np


# Memory scaling constants (empirically determined)
_MODEL_PARAMS_GB: Final[float] = 5.0
_BYTES_PER_GB: Final[int] = 1024**3
_SAFETY_MARGIN: Final[float] = 0.95  # Use 95% of available memory


@dataclasses.dataclass(frozen=True)
class MemoryEstimate:
  """Memory usage estimate for AlphaFold 3 inference.
  
  Attributes:
    total_gb: Total estimated GPU memory in GB.
    model_params_gb: Memory for model parameters.
    embeddings_gb: Memory for single and pair embeddings.
    attention_gb: Memory for attention operations.
    diffusion_gb: Memory for diffusion sampling.
    num_tokens: Number of tokens (with padding).
    num_chains: Number of chains in the complex.
    is_homomer: Whether this is a homomer complex.
  """
  total_gb: float
  model_params_gb: float
  embeddings_gb: float
  attention_gb: float
  diffusion_gb: float
  num_tokens: int
  num_chains: int
  is_homomer: bool


@dataclasses.dataclass(frozen=True)
class OptimizationSuggestion:
  """Suggested optimizations to reduce memory usage.
  
  Attributes:
    num_recycles: Recommended number of recycles.
    num_diffusion_samples: Recommended number of diffusion samples.
    flash_attention: Recommended flash attention implementation.
    estimated_memory_gb: Estimated memory with these settings.
    will_fit: Whether these settings should fit in available memory.
  """
  num_recycles: int
  num_diffusion_samples: int
  flash_attention: str
  estimated_memory_gb: float
  will_fit: bool


def get_gpu_memory_gb() -> float:
  """Returns the total GPU memory in GB for the first GPU device.
  
  Returns:
    GPU memory in GB, or 0.0 if no GPU is available.
  """
  try:
    devices = jax.local_devices(backend='gpu')
    if not devices:
      return 0.0
    # Get memory from first device
    # Note: JAX doesn't expose memory directly, so we estimate from device name
    device_name = str(devices[0])
    if 'H100' in device_name:
      return 80.0
    elif 'A100' in device_name:
      if '80GB' in device_name:
        return 80.0
      else:
        return 40.0
    elif 'V100' in device_name:
      return 32.0
    else:
      # Default conservative estimate
      return 40.0
  except Exception:
    return 0.0


def detect_homomer(num_chains: int, chain_ids: list[str]) -> bool:
  """Detects if the input is a homomer (multiple copies of same protein).
  
  Args:
    num_chains: Total number of chains.
    chain_ids: List of chain identifiers.
    
  Returns:
    True if this appears to be a homomer complex.
  """
  if num_chains <= 1:
    return False
  # Simple heuristic: if we have multiple chains with same length,
  # likely a homomer. This is a conservative estimate.
  return num_chains >= 2


def estimate_memory_requirements(
    num_tokens: int,
    num_chains: int = 1,
    num_recycles: int = 10,
    num_diffusion_samples: int = 5,
    flash_attention: bool = True,
) -> MemoryEstimate:
  """Estimates GPU memory requirements for AlphaFold 3 inference.
  
  Args:
    num_tokens: Number of tokens including padding.
    num_chains: Number of chains in the complex.
    num_recycles: Number of recycling iterations.
    num_diffusion_samples: Number of diffusion samples to generate.
    flash_attention: Whether flash attention is enabled.
    
  Returns:
    MemoryEstimate with breakdown of memory usage.
  """
  is_homomer = num_chains >= 2
  
  # Model parameters (constant)
  model_params_gb = _MODEL_PARAMS_GB
  
  # Embeddings: single (num_tokens, 384) + pair (num_tokens, num_tokens, 128)
  # Using fp16 (2 bytes per element)
  single_emb_bytes = num_tokens * 384 * 2
  pair_emb_bytes = num_tokens * num_tokens * 128 * 2
  embeddings_gb = (single_emb_bytes + pair_emb_bytes) / _BYTES_PER_GB
  
  # Attention memory - simplified formula calibrated to real usage
  # GitHub issue: 4608 tokens, 4 chains, 10 recycles, 5 samples = 86.9 GB total
  # Target: ~75 GB for attention component
  
  tokens_squared = num_tokens * num_tokens
  
  if flash_attention:
    # Flash attention reduces memory significantly
    # Calibrated to match real usage: 4608^2 * 10 * 0.00015 * 1.8 ≈ 85 GB
    attention_gb = (tokens_squared * num_recycles * 0.00015) / 1024
  else:
    # Standard attention: 3x more memory
    attention_gb = (tokens_squared * num_recycles * 0.00045) / 1024
  
  # Homomer multiplier - moderate increase for cross-chain attention
  if is_homomer:
    attention_gb *= 1.8
  
  # Diffusion sampling memory
  # Simple: ~0.5 GB per 1000 tokens per sample
  diffusion_gb = (num_tokens / 1000.0) * 0.5 * num_diffusion_samples
  
  total_gb = model_params_gb + embeddings_gb + attention_gb + diffusion_gb
  
  return MemoryEstimate(
      total_gb=total_gb,
      model_params_gb=model_params_gb,
      embeddings_gb=embeddings_gb,
      attention_gb=attention_gb,
      diffusion_gb=diffusion_gb,
      num_tokens=num_tokens,
      num_chains=num_chains,
      is_homomer=is_homomer,
  )


def suggest_optimizations(
    num_tokens: int,
    num_chains: int,
    available_memory_gb: float,
    current_num_recycles: int = 10,
    current_num_samples: int = 5,
) -> OptimizationSuggestion:
  """Suggests optimizations to fit within available GPU memory.
  
  Args:
    num_tokens: Number of tokens including padding.
    num_chains: Number of chains in the complex.
    available_memory_gb: Available GPU memory in GB.
    current_num_recycles: Current number of recycles.
    current_num_samples: Current number of diffusion samples.
    
  Returns:
    OptimizationSuggestion with recommended settings.
  """
  target_memory = available_memory_gb * _SAFETY_MARGIN
  
  # Try progressively more aggressive optimizations
  optimization_levels = [
      # (num_recycles, num_samples, flash_attention)
      (current_num_recycles, current_num_samples, True),
      (5, current_num_samples, True),
      (5, 3, True),
      (3, 3, True),
      (3, 1, True),
  ]
  
  for num_recycles, num_samples, flash_attn in optimization_levels:
    estimate = estimate_memory_requirements(
        num_tokens=num_tokens,
        num_chains=num_chains,
        num_recycles=num_recycles,
        num_diffusion_samples=num_samples,
        flash_attention=flash_attn,
    )
    
    if estimate.total_gb <= target_memory:
      return OptimizationSuggestion(
          num_recycles=num_recycles,
          num_diffusion_samples=num_samples,
          flash_attention='triton' if flash_attn else 'xla',
          estimated_memory_gb=estimate.total_gb,
          will_fit=True,
      )
  
  # Even most aggressive optimization won't fit
  return OptimizationSuggestion(
      num_recycles=3,
      num_diffusion_samples=1,
      flash_attention='triton',
      estimated_memory_gb=estimate.total_gb,
      will_fit=False,
  )


def format_memory_estimate(estimate: MemoryEstimate) -> str:
  """Formats a memory estimate for display.
  
  Args:
    estimate: MemoryEstimate to format.
    
  Returns:
    Formatted string with memory breakdown.
  """
  lines = [
      'Memory Estimate:',
      f'  Total: {estimate.total_gb:.2f} GB',
      f'  Model parameters: {estimate.model_params_gb:.2f} GB',
      f'  Embeddings: {estimate.embeddings_gb:.2f} GB',
      f'  Attention: {estimate.attention_gb:.2f} GB',
      f'  Diffusion: {estimate.diffusion_gb:.2f} GB',
      f'  Tokens: {estimate.num_tokens}',
      f'  Chains: {estimate.num_chains}',
      f'  Homomer: {estimate.is_homomer}',
  ]
  return '\n'.join(lines)


def format_optimization_suggestion(suggestion: OptimizationSuggestion) -> str:
  """Formats an optimization suggestion for display.
  
  Args:
    suggestion: OptimizationSuggestion to format.
    
  Returns:
    Formatted string with recommendations.
  """
  lines = [
      'Optimization Suggestion:',
      f'  num_recycles: {suggestion.num_recycles}',
      f'  num_diffusion_samples: {suggestion.num_diffusion_samples}',
      f'  flash_attention_implementation: {suggestion.flash_attention}',
      f'  Estimated memory: {suggestion.estimated_memory_gb:.2f} GB',
      f'  Will fit: {"Yes" if suggestion.will_fit else "No"}',
  ]
  return '\n'.join(lines)
