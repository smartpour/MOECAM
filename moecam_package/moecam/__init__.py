"""
MOECAM - Multi-Objective Evolutionary Comparison and Analysis Module
===================================================================

Main module providing CFFI interfaces to C/C++ optimization algorithms.
"""

try:
    from moecam._moecam_cffi import lib, ffi
    CFFI_AVAILABLE = True
except ImportError:
    CFFI_AVAILABLE = False
    lib = None
    ffi = None

import numpy as np
from typing import List, Tuple, Callable, Optional

def hello_moecam():
    """Test function to verify CFFI is working."""
    if not CFFI_AVAILABLE:
        return "CFFI not available"

    lib.hello_world()
    result = lib.add_numbers(5, 3)
    return f"CFFI working! 5 + 3 = {result}"

def extract_pareto_front(objectives: np.ndarray, return_indices: bool = True) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """
    Extract Pareto front from objective values using CFFI.

    Args:
        objectives: Array of shape (n_points, n_objectives) with objective values
        return_indices: Whether to return indices of Pareto points

    Returns:
        Tuple of (pareto_objectives, pareto_indices) if return_indices=True
        Otherwise just pareto_objectives
    """
    if not CFFI_AVAILABLE:
        raise RuntimeError("CFFI not available. Install package with 'pip install -e .'")

    objectives = np.asarray(objectives, dtype=np.float64)
    n_points, n_objectives = objectives.shape

    # Allocate memory for results
    max_pareto = n_points  # Upper bound
    pareto_points = np.zeros((max_pareto, n_objectives), dtype=np.float64)
    pareto_indices = np.zeros(max_pareto, dtype=np.int32)

    # Create CFFI arrays
    c_points = ffi.cast("double*", ffi.from_buffer(objectives))
    c_pareto_points = ffi.cast("double*", ffi.from_buffer(pareto_points))
    c_pareto_indices = ffi.cast("int*", ffi.from_buffer(pareto_indices))

    # Call C function
    pareto_count = lib.extract_pareto_front(
        c_points, n_points, n_objectives,
        c_pareto_points, c_pareto_indices
    )

    # Extract results
    result_pareto = pareto_points[:pareto_count].copy()
    result_indices = pareto_indices[:pareto_count].copy() if return_indices else None

    return (result_pareto, result_indices) if return_indices else result_pareto

# Mock functions for compatibility with comprehensive test
def minimize_ecam_cffi(objective_func: Callable, bounds: List[Tuple[float, float]],
                      num_objectives: int = 2, max_iterations: int = 100):
    """Mock ECAM optimization - replace with real implementation."""
    raise NotImplementedError("ECAM algorithm not yet implemented in CFFI")

def minimize_random_start_cffi(objective_func: Callable, bounds: List[Tuple[float, float]],
                              num_objectives: int = 2, max_iterations: int = 100):
    """Mock Random Start optimization - replace with real implementation."""
    raise NotImplementedError("Random Start algorithm not yet implemented in CFFI")

def calculate_hypervolume_cffi(points: np.ndarray, reference_point: np.ndarray):
    """Mock hypervolume calculation - replace with real implementation."""
    raise NotImplementedError("Hypervolume calculation not yet implemented in CFFI")

__all__ = [
    'hello_moecam',
    'extract_pareto_front',
    'minimize_ecam_cffi',
    'minimize_random_start_cffi',
    'calculate_hypervolume_cffi',
    'CFFI_AVAILABLE'
]
