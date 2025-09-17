"""
Pareto Front Extraction CFFI Interface
====================================

Proper Python-to-C++ interface for Pareto front extraction
using CFFI like the toy example.
"""

import cffi
import os
import numpy as np
from pathlib import Path

# Setup CFFI
ffi = cffi.FFI()

# Find the compiled library
_this_dir = Path(__file__).parent
_lib_path = _this_dir / 'libpareto.so'

# Define the C interface
ffi.cdef("""
// Pareto front extraction interface
int extract_pareto_front(double* all_points, int num_points, int num_objectives,
                        double* pareto_points, int max_pareto_points);
""")

# Load the library (will be None if not built yet)
try:
    lib = ffi.dlopen(str(_lib_path))
except OSError:
    lib = None

class ParetoExtractorCFFI:
    """CFFI-based Pareto front extractor."""

    def __init__(self):
        if lib is None:
            raise RuntimeError(
                f"Pareto library not found at {_lib_path}. "
                f"Run 'python build_libraries.py' to build it."
            )

    def extract_pareto_front(self, objectives):
        """Extract Pareto front using CFFI interface.

        Args:
            objectives: numpy array of shape (n_points, n_objectives)

        Returns:
            numpy array: Pareto front points, shape (n_pareto, n_objectives)
        """
        objectives = np.asarray(objectives, dtype=np.float64)

        if objectives.ndim != 2:
            raise ValueError("objectives must be 2D array")

        n_points, n_objectives = objectives.shape

        # Prepare data for C interface
        objectives_flat = np.ascontiguousarray(objectives.flatten(), dtype=np.float64)

        # Allocate space for Pareto front (worst case: all points are Pareto)
        max_pareto = n_points
        pareto_flat = np.zeros(max_pareto * n_objectives, dtype=np.float64)

        # Convert to CFFI pointers
        objectives_ptr = ffi.cast("double *", objectives_flat.ctypes.data)
        pareto_ptr = ffi.cast("double *", pareto_flat.ctypes.data)

        # Call C function
        num_pareto = lib.extract_pareto_front(objectives_ptr, n_points, n_objectives,
                                             pareto_ptr, max_pareto)

        if num_pareto < 0:
            raise RuntimeError("Pareto front extraction failed")

        if num_pareto == 0:
            return np.empty((0, n_objectives), dtype=np.float64)

        # Reshape result
        pareto_front = pareto_flat[:num_pareto * n_objectives].reshape(num_pareto, n_objectives)

        return pareto_front

    def extract_pareto_mask(self, objectives):
        """Extract Pareto front as boolean mask.

        Args:
            objectives: numpy array of shape (n_points, n_objectives)

        Returns:
            numpy array: Boolean mask indicating which points are Pareto optimal
        """
        objectives = np.asarray(objectives, dtype=np.float64)
        pareto_front = self.extract_pareto_front(objectives)

        # Create boolean mask by finding which original points match Pareto front
        mask = np.zeros(len(objectives), dtype=bool)

        # For each Pareto point, find its index in the original array
        for pareto_point in pareto_front:
            # Find matching point in original objectives (within tolerance)
            distances = np.sum((objectives - pareto_point)**2, axis=1)
            closest_idx = np.argmin(distances)
            if distances[closest_idx] < 1e-10:  # Very close match
                mask[closest_idx] = True

        return mask

# Convenience function
def extract_pareto_front_cffi(objectives):
    """Extract Pareto front using CFFI interface."""
    extractor = ParetoExtractorCFFI()
    return extractor.extract_pareto_front(objectives)
