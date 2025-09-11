#!/bin/bash
# Enhanced Docker entrypoint script for AlphaFold3
# Addresses threading and OpenBLAS configuration issues
# Resolves Issue #83: Docker runtime failures with threading problems

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}AlphaFold3 Docker Container - Enhanced Startup${NC}"
echo "=========================================="

# Function to detect available CPU cores and memory
detect_system_resources() {
    local cpu_cores=$(nproc)
    local total_memory=$(grep MemTotal /proc/meminfo | awk '{print int($2/1024/1024)}')
    
    echo -e "${BLUE}System Resources Detected:${NC}"
    echo "  CPU Cores: $cpu_cores"
    echo "  Total Memory: ${total_memory}GB"
    
    # Set conservative thread limits based on available resources
    if [ $cpu_cores -le 4 ]; then
        export RECOMMENDED_THREADS=1
    elif [ $cpu_cores -le 8 ]; then
        export RECOMMENDED_THREADS=2
    elif [ $cpu_cores -le 16 ]; then
        export RECOMMENDED_THREADS=4
    else
        export RECOMMENDED_THREADS=8
    fi
    
    echo "  Recommended Threads: $RECOMMENDED_THREADS"
}

# Function to configure OpenBLAS and threading
configure_threading() {
    echo -e "${YELLOW}Configuring Threading and OpenBLAS...${NC}"
    
    # Set OpenBLAS thread count if not already set
    if [ -z "$OPENBLAS_NUM_THREADS" ]; then
        export OPENBLAS_NUM_THREADS=$RECOMMENDED_THREADS
        echo "  Set OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"
    else
        echo "  Using existing OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"
    fi
    
    # Set MKL thread count if Intel MKL is present
    if [ -z "$MKL_NUM_THREADS" ]; then
        export MKL_NUM_THREADS=$RECOMMENDED_THREADS
        echo "  Set MKL_NUM_THREADS=$MKL_NUM_THREADS"
    fi
    
    # Set OpenMP thread count
    if [ -z "$OMP_NUM_THREADS" ]; then
        export OMP_NUM_THREADS=$RECOMMENDED_THREADS
        echo "  Set OMP_NUM_THREADS=$OMP_NUM_THREADS"
    fi
    
    # Set NumPy threading
    if [ -z "$NUMBA_NUM_THREADS" ]; then
        export NUMBA_NUM_THREADS=$RECOMMENDED_THREADS
        echo "  Set NUMBA_NUM_THREADS=$NUMBA_NUM_THREADS"
    fi
    
    # JAX threading configuration
    if [ -z "$XLA_FLAGS" ]; then
        export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false --xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
        echo "  Set XLA_FLAGS for CUDA 7.x compatibility"
    fi
}

# Function to check system limits and adjust if needed
check_system_limits() {
    echo -e "${YELLOW}Checking System Limits...${NC}"
    
    # Check current limits
    local max_processes=$(ulimit -u)
    local max_file_desc=$(ulimit -n)
    local stack_size=$(ulimit -s)
    
    echo "  Max Processes: $max_processes"
    echo "  Max File Descriptors: $max_file_desc"
    echo "  Stack Size: $stack_size KB"
    
    # Set minimum required limits if possible
    if [ "$max_processes" -lt 4096 ]; then
        echo -e "  ${YELLOW}Warning: Low process limit ($max_processes). Consider increasing with 'ulimit -u 8192'${NC}"
    fi
    
    if [ "$max_file_desc" -lt 1024 ]; then
        echo -e "  ${YELLOW}Warning: Low file descriptor limit ($max_file_desc). Consider increasing with 'ulimit -n 4096'${NC}"
    fi
}

# Function to validate Python and NumPy installation
validate_python_environment() {
    echo -e "${YELLOW}Validating Python Environment...${NC}"
    
    # Check Python version
    python_version=$(python3 --version 2>&1 | cut -d' ' -f2)
    echo "  Python Version: $python_version"
    
    # Test NumPy import with threading configuration
    if python3 -c "import numpy as np; print(f'NumPy Version: {np.__version__}'); print(f'NumPy Config: {np.show_config()}')" > /tmp/numpy_test.log 2>&1; then
        echo -e "  ${GREEN}✓ NumPy import successful${NC}"
    else
        echo -e "  ${RED}✗ NumPy import failed. See details:${NC}"
        cat /tmp/numpy_test.log
        echo -e "${YELLOW}  Attempting to fix NumPy threading issues...${NC}"
        
        # Try to fix common NumPy threading issues
        export OPENBLAS_NUM_THREADS=1
        export MKL_NUM_THREADS=1
        export OMP_NUM_THREADS=1
        
        if python3 -c "import numpy as np; print('NumPy import successful after threading fix')" > /tmp/numpy_fix_test.log 2>&1; then
            echo -e "  ${GREEN}✓ NumPy import fixed with single-threading${NC}"
        else
            echo -e "  ${RED}✗ NumPy import still failing:${NC}"
            cat /tmp/numpy_fix_test.log
            exit 1
        fi
    fi
}

# Function to test JAX functionality
test_jax_functionality() {
    echo -e "${YELLOW}Testing JAX Functionality...${NC}"
    
    if python3 -c "import jax; import jax.numpy as jnp; print(f'JAX Version: {jax.__version__}'); print(f'JAX Devices: {jax.devices()}'); x = jnp.array([1, 2, 3]); print(f'JAX Test: {x.sum()}')" > /tmp/jax_test.log 2>&1; then
        echo -e "  ${GREEN}✓ JAX functionality test passed${NC}"
    else
        echo -e "  ${YELLOW}! JAX test had issues (may still work for AF3):${NC}"
        cat /tmp/jax_test.log
    fi
}

# Function to display final configuration
display_configuration() {
    echo -e "${BLUE}Final Configuration:${NC}"
    echo "  OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"
    echo "  MKL_NUM_THREADS=${MKL_NUM_THREADS:-not set}"
    echo "  OMP_NUM_THREADS=$OMP_NUM_THREADS"
    echo "  NUMBA_NUM_THREADS=${NUMBA_NUM_THREADS:-not set}"
    echo "  XLA_FLAGS=$XLA_FLAGS"
}

# Main execution flow
main() {
    echo -e "${GREEN}Starting AlphaFold3 Docker container with enhanced configuration...${NC}"
    
    # Step 1: Detect system resources
    detect_system_resources
    
    # Step 2: Configure threading
    configure_threading
    
    # Step 3: Check system limits
    check_system_limits
    
    # Step 4: Validate Python environment
    validate_python_environment
    
    # Step 5: Test JAX functionality
    test_jax_functionality
    
    # Step 6: Display final configuration
    display_configuration
    
    echo -e "${GREEN}Container initialization completed successfully!${NC}"
    echo "=========================================="
    
    # Execute the original command or start AlphaFold3
    if [ $# -eq 0 ]; then
        echo -e "${BLUE}Starting interactive shell...${NC}"
        exec /bin/bash
    else
        echo -e "${BLUE}Executing command: $@${NC}"
        exec "$@"
    fi
}

# Handle script interruption
trap 'echo -e "\n${RED}Container startup interrupted${NC}"; exit 1' INT TERM

# Run main function
main "$@"