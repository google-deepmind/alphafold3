#!/usr/bin/env python3
"""
Enhanced AlphaFold3 runner script with threading fixes and diagnostics.
Addresses Issue #83: Docker runtime failures with OpenBLAS threading.

This script provides an enhanced wrapper around the AlphaFold3 inference
with automatic threading configuration and error diagnostics.
"""

import os
import sys
import subprocess
import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

def setup_threading_environment():
    """Configure threading environment variables for stable execution."""
    
    # Get system CPU count
    try:
        import multiprocessing
        cpu_count = multiprocessing.cpu_count()
    except:
        cpu_count = 4  # Safe fallback
    
    # Calculate optimal thread counts
    if cpu_count <= 4:
        optimal_threads = 1
    elif cpu_count <= 8:
        optimal_threads = 2
    elif cpu_count <= 16:
        optimal_threads = 4
    else:
        optimal_threads = 8
    
    # Set threading environment variables if not already set
    threading_config = {
        'OPENBLAS_NUM_THREADS': str(optimal_threads),
        'MKL_NUM_THREADS': str(optimal_threads),
        'OMP_NUM_THREADS': str(optimal_threads), 
        'NUMBA_NUM_THREADS': str(optimal_threads),
        'XLA_FLAGS': '--xla_cpu_multi_thread_eigen=false --xla_disable_hlo_passes=custom-kernel-fusion-rewriter'
    }
    
    applied_configs = []
    for key, value in threading_config.items():
        if key not in os.environ:
            os.environ[key] = value
            applied_configs.append(f"{key}={value}")
    
    if applied_configs:
        print(f"Applied threading configuration: {', '.join(applied_configs)}")
    
    return optimal_threads

def test_dependencies():
    """Test that required dependencies can be imported."""
    
    print("Testing dependencies...")
    
    # Test critical imports
    dependencies = [
        ('numpy', 'NumPy'),
        ('jax', 'JAX'), 
        ('haiku', 'Haiku'),
        ('ml_collections', 'ML Collections')
    ]
    
    failed_imports = []
    
    for module_name, display_name in dependencies:
        try:
            __import__(module_name)
            print(f"✓ {display_name} import successful")
        except ImportError as e:
            print(f"✗ {display_name} import failed: {e}")
            failed_imports.append(display_name)
        except Exception as e:
            print(f"⚠ {display_name} import had issues: {e}")
    
    if failed_imports:
        print(f"Failed to import: {', '.join(failed_imports)}")
        return False
    
    return True

def check_gpu_availability():
    """Check GPU availability and CUDA setup."""
    
    print("Checking GPU availability...")
    
    try:
        import jax
        devices = jax.devices()
        gpu_devices = [d for d in devices if d.device_kind == 'gpu']
        
        if gpu_devices:
            print(f"✓ Found {len(gpu_devices)} GPU device(s)")
            for i, device in enumerate(gpu_devices):
                print(f"  GPU {i}: {device}")
        else:
            print("⚠ No GPU devices found, will use CPU")
            
        return len(gpu_devices) > 0
        
    except Exception as e:
        print(f"⚠ Could not check GPU availability: {e}")
        return False

def validate_input_format(input_path: str) -> bool:
    """Validate that the input file is in the correct format."""
    
    if not os.path.exists(input_path):
        print(f"✗ Input file not found: {input_path}")
        return False
    
    try:
        with open(input_path, 'r') as f:
            if input_path.endswith('.json'):
                data = json.load(f)
                if 'sequences' not in data:
                    print("⚠ Input JSON missing 'sequences' field")
                    return False
            elif input_path.endswith('.fasta'):
                content = f.read()
                if not content.strip().startswith('>'):
                    print("⚠ Input file doesn't appear to be valid FASTA format")
                    return False
            else:
                print("⚠ Input file should be .json or .fasta format")
                
        print(f"✓ Input file format appears valid: {input_path}")
        return True
        
    except Exception as e:
        print(f"✗ Error validating input file: {e}")
        return False

def run_alphafold3_inference(
    input_path: str,
    output_dir: str,
    model_paths: str,
    databases_path: str,
    additional_args: List[str] = None
) -> Tuple[bool, str]:
    """
    Run AlphaFold3 inference with enhanced error handling.
    
    Returns:
        Tuple of (success: bool, output: str)
    """
    
    # Build command
    cmd = [
        sys.executable, 
        'run_alphafold.py',
        f'--json_path={input_path}',
        f'--output_dir={output_dir}',
        f'--model_paths={model_paths}',
        f'--db_dir={databases_path}'
    ]
    
    if additional_args:
        cmd.extend(additional_args)
    
    print(f"Running command: {' '.join(cmd)}")
    print("=" * 60)
    
    try:
        # Run with real-time output
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            bufsize=1
        )
        
        output_lines = []
        
        # Stream output in real-time
        for line in iter(process.stdout.readline, ''):
            print(line, end='')
            output_lines.append(line)
        
        process.stdout.close()
        return_code = process.wait()
        
        full_output = ''.join(output_lines)
        
        if return_code == 0:
            print("=" * 60)
            print("✓ AlphaFold3 inference completed successfully!")
            return True, full_output
        else:
            print("=" * 60)
            print(f"✗ AlphaFold3 inference failed with return code: {return_code}")
            return False, full_output
            
    except Exception as e:
        error_msg = f"Error running AlphaFold3: {e}"
        print(error_msg)
        return False, error_msg

def diagnose_threading_issues(error_output: str) -> List[str]:
    """Analyze error output for threading-related issues and provide suggestions."""
    
    suggestions = []
    
    # Common threading error patterns
    threading_errors = [
        ('pthread_create failed', 'Thread creation failure detected. Try reducing OPENBLAS_NUM_THREADS to 1.'),
        ('OpenBLAS', 'OpenBLAS threading issue detected. Set OPENBLAS_NUM_THREADS=1.'),
        ('MKL', 'Intel MKL threading issue detected. Set MKL_NUM_THREADS=1.'),
        ('numpy', 'NumPy import issue detected. Check threading configuration.'),
        ('out of memory', 'Memory issue detected. Consider reducing batch size or thread count.'),
        ('CUDA', 'CUDA issue detected. Check GPU memory and drivers.')
    ]
    
    for pattern, suggestion in threading_errors:
        if pattern.lower() in error_output.lower():
            suggestions.append(suggestion)
    
    return suggestions

def main():
    """Main function to run AlphaFold3 with enhanced configuration."""
    
    parser = argparse.ArgumentParser(
        description='Enhanced AlphaFold3 runner with threading fixes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 run_alphafold3_fixed.py --input input.json --output /tmp/af3_output --models /path/to/models --databases /path/to/db
  python3 run_alphafold3_fixed.py --input protein.fasta --output ./results --models ./models --databases ./databases --max_sequences 256
        """
    )
    
    parser.add_argument('--input', required=True, help='Input file (JSON or FASTA)')
    parser.add_argument('--output', required=True, help='Output directory')  
    parser.add_argument('--models', required=True, help='Path to model files')
    parser.add_argument('--databases', required=True, help='Path to databases')
    parser.add_argument('--max_sequences', type=int, help='Maximum sequences to use')
    parser.add_argument('--num_diffn_timesteps', type=int, help='Number of diffusion timesteps')
    parser.add_argument('--skip_tests', action='store_true', help='Skip dependency tests')
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    
    args = parser.parse_args()
    
    print("AlphaFold3 Enhanced Runner")
    print("=" * 50)
    
    # Step 1: Setup threading environment
    optimal_threads = setup_threading_environment()
    
    # Step 2: Run tests unless skipped
    if not args.skip_tests:
        print("\nRunning pre-flight checks...")
        
        if not test_dependencies():
            print("✗ Dependency tests failed. Use --skip_tests to bypass.")
            return 1
            
        check_gpu_availability()
        
        if not validate_input_format(args.input):
            print("✗ Input validation failed.")
            return 1
    
    # Step 3: Prepare additional arguments
    additional_args = []
    
    if args.max_sequences:
        additional_args.append(f'--max_sequences={args.max_sequences}')
        
    if args.num_diffn_timesteps:
        additional_args.append(f'--num_diffn_timesteps={args.num_diffn_timesteps}')
    
    if args.debug:
        additional_args.extend(['--alsologtostderr', '--verbosity=1'])
    
    # Step 4: Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Step 5: Run AlphaFold3
    print(f"\nStarting AlphaFold3 inference...")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print(f"Models: {args.models}")
    print(f"Databases: {args.databases}")
    print(f"Threading: {optimal_threads} threads")
    
    success, output = run_alphafold3_inference(
        args.input,
        args.output, 
        args.models,
        args.databases,
        additional_args
    )
    
    # Step 6: Handle results
    if success:
        print(f"\n✓ Results saved to: {args.output}")
        return 0
    else:
        print(f"\n✗ AlphaFold3 inference failed")
        
        # Provide diagnostic suggestions
        suggestions = diagnose_threading_issues(output)
        if suggestions:
            print("\nSuggested fixes:")
            for suggestion in suggestions:
                print(f"  • {suggestion}")
        
        print("\nFor additional help:")
        print("  • Check system resources with 'ulimit -a'")
        print("  • Try running with --debug for more information")
        print("  • Consider reducing --max_sequences")
        
        return 1

if __name__ == "__main__":
    sys.exit(main())