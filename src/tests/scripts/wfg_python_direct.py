"""
Python binding for WFG hypervolume calculation using direct C library integration.
This provides a fast, direct interface to the WFG algorithm without subprocess overhead.
"""

import ctypes
import numpy as np
import os

class WFGDirectInterface:
    """Direct C library interface for WFG hypervolume calculation."""

    def __init__(self, library_path=None):
        """Initialize the WFG interface.

        Args:
            library_path: Path to libwfg_simple.so. If None, looks in current directory.
        """
        if library_path is None:
            # Look for library in the same directory as this script
            script_dir = os.path.dirname(os.path.abspath(__file__))
            library_path = os.path.join(script_dir, 'libwfg_simple.so')

        if not os.path.exists(library_path):
            raise FileNotFoundError(f"WFG library not found at {library_path}")

        # Load the shared library
        self.lib = ctypes.CDLL(library_path)

        # Define function signatures
        self.lib.wfg_simple_init.argtypes = [ctypes.c_int, ctypes.c_int]
        self.lib.wfg_simple_init.restype = ctypes.c_int

        self.lib.wfg_simple_calculate.argtypes = [
            ctypes.POINTER(ctypes.c_double),  # pts
            ctypes.POINTER(ctypes.c_double),  # nadir
            ctypes.c_int,                     # nPoints
            ctypes.c_int                      # dimensions
        ]
        self.lib.wfg_simple_calculate.restype = ctypes.c_double

        self.lib.wfg_simple_cleanup.argtypes = []
        self.lib.wfg_simple_cleanup.restype = None

        self.initialized = False

    def initialize(self, dimensions, max_points=1000):
        """Initialize the WFG algorithm for given dimensions.

        Args:
            dimensions: Number of objectives
            max_points: Maximum number of points to handle

        Returns:
            bool: True if initialization successful
        """
        result = self.lib.wfg_simple_init(dimensions, max_points)
        self.initialized = (result == 1)
        return self.initialized

    def calculate_hypervolume(self, points, reference_point):
        """Calculate hypervolume for given points and reference point.

        Args:
            points: 2D numpy array of shape (n_points, n_objectives)
            reference_point: 1D numpy array of shape (n_objectives,)

        Returns:
            float: Hypervolume value
        """
        # Convert inputs to numpy arrays
        points = np.asarray(points, dtype=np.float64)
        reference_point = np.asarray(reference_point, dtype=np.float64)

        # Validate dimensions
        if points.ndim != 2:
            raise ValueError("Points must be a 2D array")

        n_points, n_objectives = points.shape

        if len(reference_point) != n_objectives:
            raise ValueError("Reference point dimension must match number of objectives")

        # Initialize if not done yet
        if not self.initialized:
            if not self.initialize(n_objectives, n_points + 100):
                raise RuntimeError("Failed to initialize WFG algorithm")

        # Flatten points array for C interface
        points_flat = points.flatten()

        # Create ctypes arrays
        points_c = (ctypes.c_double * len(points_flat))(*points_flat)
        ref_c = (ctypes.c_double * len(reference_point))(*reference_point)

        # Call C function
        result = self.lib.wfg_simple_calculate(
            points_c,
            ref_c,
            n_points,
            n_objectives
        )

        if result < 0:
            raise RuntimeError("WFG calculation failed")

        return result

    def cleanup(self):
        """Clean up WFG resources."""
        if hasattr(self, 'lib'):
            self.lib.wfg_simple_cleanup()
        self.initialized = False

    def __del__(self):
        """Cleanup on destruction."""
        self.cleanup()

# Convenience function for one-shot hypervolume calculation
def calculate_hypervolume(points, reference_point, library_path=None):
    """Calculate hypervolume using WFG algorithm.

    Args:
        points: 2D array-like of shape (n_points, n_objectives)
        reference_point: 1D array-like of shape (n_objectives,)
        library_path: Optional path to WFG library

    Returns:
        float: Hypervolume value
    """
    wfg = WFGDirectInterface(library_path)
    try:
        return wfg.calculate_hypervolume(points, reference_point)
    finally:
        wfg.cleanup()

if __name__ == "__main__":
    # Test the interface
    print("Testing WFG Direct Python Interface")

    # Test 1: Simple 2D case
    points_2d = np.array([
        [2.0, 8.0],
        [4.0, 6.0],
        [6.0, 4.0],
        [8.0, 2.0]
    ])
    ref_2d = np.array([10.0, 10.0])

    hv_2d = calculate_hypervolume(points_2d, ref_2d)
    print(f"2D Hypervolume: {hv_2d}")

    # Test 2: 3D case
    points_3d = np.array([
        [1.0, 8.0, 9.0],
        [2.0, 7.0, 8.0],
        [3.0, 6.0, 7.0],
        [4.0, 5.0, 6.0]
    ])
    ref_3d = np.array([10.0, 10.0, 10.0])

    hv_3d = calculate_hypervolume(points_3d, ref_3d)
    print(f"3D Hypervolume: {hv_3d}")

    # Test 3: Comparison with Kursawe-like data
    points_kur = np.array([
        [-19.0, 2.0],
        [-18.5, 3.0],
        [-17.8, 4.5],
        [-16.2, 7.0],
        [-14.5, 10.5]
    ])
    ref_kur = np.array([0.0, 15.0])

    hv_kur = calculate_hypervolume(points_kur, ref_kur)
    print(f"Kursawe-like Hypervolume: {hv_kur}")

    print("All tests completed successfully!")
