"""
MOECAM Direct C++ Library Interface
===================================

This module provides high-performance direct C++ library interfaces for:
1. Pareto front extraction
2. Hypervolume calculation using WFG algorithm

Key advantages over subprocess-based approach:
- ~10x faster execution (milliseconds vs seconds)
- No subprocess overhead
- Direct memory management
- Seamless integration with Python

Example usage:
    from moecam_direct import MOECAMDirect

    moecam = MOECAMDirect()

    # Your optimization algorithm generates these
    objectives = [[f1_1, f2_1], [f1_2, f2_2], ...]

    # Extract Pareto front and calculate hypervolume
    result = moecam.analyze(objectives)

    print(f"Pareto points: {len(result['pareto_front'])}")
    print(f"Hypervolume: {result['hypervolume']}")
"""

import sys
import os
import ctypes
import tempfile
import numpy as np
from pathlib import Path

# Add local modules to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from pareto_interface import ParetoFrontExtractor

class DirectWFGInterface:
    """High-performance direct C++ interface to WFG hypervolume algorithm."""

    def __init__(self, library_path=None):
        """Initialize WFG interface.

        Args:
            library_path: Path to libwfg_simple.so (auto-detected if None)
        """
        if library_path is None:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            library_path = os.path.join(script_dir, 'libwfg_simple.so')

        if not os.path.exists(library_path):
            raise FileNotFoundError(
                f"WFG library not found at {library_path}. "
                f"Run 'make libwfg_simple.so' to build it."
            )

        self.lib = ctypes.CDLL(library_path)
        self._setup_function_signatures()
        self.initialized = False

    def _setup_function_signatures(self):
        """Define C function signatures for ctypes."""
        self.lib.wfg_simple_init.argtypes = [ctypes.c_int, ctypes.c_int]
        self.lib.wfg_simple_init.restype = ctypes.c_int

        self.lib.wfg_simple_calculate.argtypes = [
            ctypes.POINTER(ctypes.c_double),  # points array
            ctypes.POINTER(ctypes.c_double),  # reference point
            ctypes.c_int,                     # number of points
            ctypes.c_int                      # number of objectives
        ]
        self.lib.wfg_simple_calculate.restype = ctypes.c_double

        self.lib.wfg_simple_cleanup.argtypes = []
        self.lib.wfg_simple_cleanup.restype = None

    def calculate_hypervolume(self, points, reference_point=None):
        """Calculate hypervolume for given points.

        Args:
            points: List of lists or numpy array, shape (n_points, n_objectives)
            reference_point: Reference point for hypervolume calculation.
                           If None, uses nadir point + small offset.

        Returns:
            float: Hypervolume value
        """
        # Convert to list of lists for consistency
        if isinstance(points, np.ndarray):
            points = points.tolist()

        if not points:
            return 0.0

        n_points = len(points)
        n_objectives = len(points[0])

        # Auto-generate reference point if not provided
        if reference_point is None:
            reference_point = self._calculate_nadir_reference(points)

        # Initialize WFG if needed
        if not self.initialized:
            result = self.lib.wfg_simple_init(n_objectives, n_points + 100)
            if result != 1:
                raise RuntimeError("Failed to initialize WFG algorithm")
            self.initialized = True

        # Prepare data for C interface
        points_flat = []
        for point in points:
            points_flat.extend(point)

        points_c = (ctypes.c_double * len(points_flat))(*points_flat)
        ref_c = (ctypes.c_double * len(reference_point))(*reference_point)

        # Calculate hypervolume
        result = self.lib.wfg_simple_calculate(
            points_c, ref_c, n_points, n_objectives
        )

        if result < 0:
            raise RuntimeError("WFG calculation failed")

        return result

    def _calculate_nadir_reference(self, points):
        """Calculate nadir point + offset for reference point."""
        points_array = np.array(points)
        nadir = np.max(points_array, axis=0)

        # Add small offset to nadir point
        offset = 0.1 * np.abs(nadir) + 1.0
        return (nadir + offset).tolist()

    def cleanup(self):
        """Clean up WFG resources."""
        if hasattr(self, 'lib') and self.initialized:
            self.lib.wfg_simple_cleanup()
        self.initialized = False

    def __del__(self):
        """Cleanup on destruction."""
        self.cleanup()

class MOECAMDirect:
    """Complete MOECAM interface using direct C++ libraries.

    Provides integrated Pareto front extraction and hypervolume calculation
    with maximum performance through direct C++ library calls.
    """

    def __init__(self, wfg_library_path=None):
        """Initialize MOECAM direct interface.

        Args:
            wfg_library_path: Path to WFG shared library (auto-detected if None)
        """
        self.pareto_extractor = ParetoFrontExtractor()
        self.wfg_interface = DirectWFGInterface(wfg_library_path)

    def analyze(self, objectives, reference_point=None, include_details=False):
        """Complete analysis: Pareto extraction + hypervolume calculation.

        Args:
            objectives: Array-like of shape (n_points, n_objectives)
            reference_point: Reference point for hypervolume (auto if None)
            include_details: If True, include detailed timing and statistics

        Returns:
            dict: {
                'pareto_front': numpy array of Pareto optimal points,
                'hypervolume': float,
                'summary': dict with statistics (if include_details=True)
            }
        """
        import time

        start_total = time.time()

        # Step 1: Extract Pareto front
        start_pareto = time.time()
        pareto_result = self.pareto_extractor.extract_pareto_front(
            np.array(objectives), include_indices=True
        )
        pareto_time = time.time() - start_pareto

        pareto_front = pareto_result['pareto_front']

        # Step 2: Calculate hypervolume
        start_hv = time.time()
        if len(pareto_front) > 0:
            hypervolume = self.wfg_interface.calculate_hypervolume(
                pareto_front, reference_point
            )
        else:
            hypervolume = 0.0
        hv_time = time.time() - start_hv

        total_time = time.time() - start_total

        result = {
            'pareto_front': pareto_front,
            'hypervolume': hypervolume
        }

        if include_details:
            result['summary'] = {
                'original_points': len(objectives),
                'pareto_points': len(pareto_front),
                'reduction_ratio': (len(objectives) - len(pareto_front)) / len(objectives),
                'pareto_time': pareto_time,
                'hypervolume_time': hv_time,
                'total_time': total_time,
                'pareto_indices': pareto_result.get('indices', [])
            }

        return result

    def batch_analyze(self, objectives_list, reference_points=None):
        """Analyze multiple objective sets efficiently.

        Args:
            objectives_list: List of objective arrays
            reference_points: List of reference points (one per objective set)

        Returns:
            list: Results for each objective set
        """
        results = []

        for i, objectives in enumerate(objectives_list):
            ref_point = None
            if reference_points is not None:
                ref_point = reference_points[i]

            result = self.analyze(objectives, ref_point, include_details=True)
            results.append(result)

        return results

    def cleanup(self):
        """Clean up resources."""
        if hasattr(self, 'wfg_interface'):
            self.wfg_interface.cleanup()

# Convenience functions for quick usage
def extract_pareto_front(objectives):
    """Extract Pareto front from objectives.

    Args:
        objectives: Array-like of shape (n_points, n_objectives)

    Returns:
        numpy.ndarray: Pareto optimal points
    """
    extractor = ParetoFrontExtractor()
    result = extractor.extract_pareto_front(np.array(objectives))
    return result['pareto_front']

def calculate_hypervolume(objectives, reference_point=None):
    """Calculate hypervolume for objectives.

    Args:
        objectives: Array-like of shape (n_points, n_objectives)
        reference_point: Reference point (auto-calculated if None)

    Returns:
        float: Hypervolume value
    """
    wfg = DirectWFGInterface()
    try:
        return wfg.calculate_hypervolume(objectives, reference_point)
    finally:
        wfg.cleanup()

def complete_analysis(objectives, reference_point=None):
    """Complete Pareto + hypervolume analysis.

    Args:
        objectives: Array-like of shape (n_points, n_objectives)
        reference_point: Reference point (auto-calculated if None)

    Returns:
        dict: Complete analysis results
    """
    moecam = MOECAMDirect()
    try:
        return moecam.analyze(objectives, reference_point, include_details=True)
    finally:
        moecam.cleanup()

if __name__ == "__main__":
    # Demo and validation
    print("MOECAM Direct C++ Library Interface Demo")
    print("========================================")

    # Generate test data (Kursawe function)
    import math
    import random

    random.seed(42)

    def kursawe(x):
        f1 = sum(-10 * math.exp(-0.2 * math.sqrt(x[i]**2 + x[i+1]**2)) for i in range(len(x)-1))
        f2 = sum(abs(xi)**0.8 + 5 * math.sin(xi**3) for xi in x)
        return [f1, f2]

    # Generate 200 test points
    test_objectives = []
    for _ in range(200):
        x = [random.uniform(-5, 5) for _ in range(3)]
        test_objectives.append(kursawe(x))

    print(f"\nAnalyzing {len(test_objectives)} objective vectors...")

    # Complete analysis
    result = complete_analysis(test_objectives)

    print(f"\nResults:")
    print(f"  Original points: {result['summary']['original_points']}")
    print(f"  Pareto points: {result['summary']['pareto_points']}")
    print(f"  Reduction: {result['summary']['reduction_ratio']:.1%}")
    print(f"  Hypervolume: {result['hypervolume']:.6f}")
    print(f"  Total time: {result['summary']['total_time']:.4f}s")
    print(f"    - Pareto extraction: {result['summary']['pareto_time']:.4f}s")
    print(f"    - Hypervolume calc: {result['summary']['hypervolume_time']:.6f}s")

    print(f"\n✅ MOECAM Direct Interface ready for use!")
