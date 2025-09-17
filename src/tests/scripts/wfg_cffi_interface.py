"""
WFG Hypervolume CFFI Interface
=============================

Python-to-C++ interface for WFG hypervolume calculation
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
_lib_path = _this_dir / 'libwfg.so'

# Define the C interface
ffi.cdef("""
// WFG hypervolume interface
int wfg_init(int num_objectives, int max_points);
double wfg_calculate_hypervolume(double* points, double* reference_point, int num_points, int num_objectives);
void wfg_cleanup();
""")

# Load the library (will be None if not built yet)
try:
    lib = ffi.dlopen(str(_lib_path))
except OSError:
    lib = None

class WFGHypervolumeCFFI:
    """CFFI-based WFG hypervolume calculator."""

    def __init__(self):
        if lib is None:
            raise RuntimeError(
                f"WFG library not found at {_lib_path}. "
                f"Run 'python build_libraries.py' to build it."
            )
        self.initialized = False

    def calculate_hypervolume(self, points, reference_point):
        """Calculate hypervolume using CFFI interface.

        Args:
            points: numpy array of shape (n_points, n_objectives)
            reference_point: numpy array of shape (n_objectives,)

        Returns:
            float: hypervolume value
        """
        points = np.asarray(points, dtype=np.float64)
        reference_point = np.asarray(reference_point, dtype=np.float64)

        if points.ndim != 2:
            raise ValueError("points must be 2D array")

        n_points, n_objectives = points.shape

        if reference_point.shape[0] != n_objectives:
            raise ValueError("reference_point must match number of objectives")

        # Initialize if needed
        if not self.initialized:
            result = lib.wfg_init(n_objectives, n_points * 2)  # Some buffer
            if result != 0:
                raise RuntimeError("Failed to initialize WFG")
            self.initialized = True

        # Prepare data for C interface
        points_flat = np.ascontiguousarray(points.flatten(), dtype=np.float64)
        reference_flat = np.ascontiguousarray(reference_point, dtype=np.float64)

        # Convert to CFFI pointers
        points_ptr = ffi.cast("double *", points_flat.ctypes.data)
        ref_ptr = ffi.cast("double *", reference_flat.ctypes.data)

        # Call C function
        hypervolume = lib.wfg_calculate_hypervolume(points_ptr, ref_ptr, n_points, n_objectives)

        if hypervolume < 0:
            raise RuntimeError("WFG hypervolume calculation failed")

        return float(hypervolume)

    def __del__(self):
        """Cleanup WFG resources."""
        if lib is not None and self.initialized:
            lib.wfg_cleanup()

# Convenience function
def calculate_hypervolume_cffi(points, reference_point):
    """Calculate hypervolume using CFFI interface."""
    calculator = WFGHypervolumeCFFI()
    return calculator.calculate_hypervolume(points, reference_point)
