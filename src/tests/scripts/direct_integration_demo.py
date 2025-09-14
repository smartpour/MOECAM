"""
Complete MOECAM integration demo using direct C++ library interfaces.
This demonstrates the full pipeline: objective recording -> Pareto extraction -> hypervolume calculation
using direct C++ library calls instead of subprocess overhead.
"""

import sys
import os
import ctypes
import time
from pathlib import Path

# Add the scripts directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Import our existing interfaces
from pareto_interface import ParetoFrontExtractor
from objective_recorder import ObjectiveRecorder

class DirectWFGInterface:
    """Direct C library interface for WFG hypervolume calculation."""

    def __init__(self, library_path=None):
        if library_path is None:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            library_path = os.path.join(script_dir, 'libwfg_simple.so')

        if not os.path.exists(library_path):
            raise FileNotFoundError(f"WFG library not found at {library_path}")

        self.lib = ctypes.CDLL(library_path)

        # Define function signatures
        self.lib.wfg_simple_init.argtypes = [ctypes.c_int, ctypes.c_int]
        self.lib.wfg_simple_init.restype = ctypes.c_int

        self.lib.wfg_simple_calculate.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_int,
            ctypes.c_int
        ]
        self.lib.wfg_simple_calculate.restype = ctypes.c_double

        self.lib.wfg_simple_cleanup.argtypes = []
        self.lib.wfg_simple_cleanup.restype = None

        self.initialized = False

    def calculate_hypervolume(self, points_2d_list, reference_point):
        """Calculate hypervolume for list of points.

        Args:
            points_2d_list: List of lists, each inner list is one point's objectives
            reference_point: List of reference point coordinates

        Returns:
            float: Hypervolume value
        """
        if not points_2d_list:
            return 0.0

        n_points = len(points_2d_list)
        n_objectives = len(points_2d_list[0])

        # Validate reference point
        if len(reference_point) != n_objectives:
            raise ValueError("Reference point dimension must match number of objectives")

        # Initialize if needed
        if not self.initialized:
            result = self.lib.wfg_simple_init(n_objectives, n_points + 100)
            if result != 1:
                raise RuntimeError("Failed to initialize WFG algorithm")
            self.initialized = True

        # Flatten points array for C interface
        points_flat = []
        for point in points_2d_list:
            points_flat.extend(point)

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
        """Clean up resources."""
        if hasattr(self, 'lib') and self.initialized:
            self.lib.wfg_simple_cleanup()
        self.initialized = False

    def __del__(self):
        self.cleanup()

def demonstrate_direct_pipeline():
    """Demonstrate the complete pipeline using direct C++ libraries."""

    print("=== MOECAM Direct C++ Library Integration Demo ===")
    print("This demo shows objective recording -> Pareto extraction -> hypervolume calculation")
    print("using direct C++ library calls for maximum performance.\n")

    # Step 1: Record objectives (simulating optimization algorithm)
    print("Step 1: Recording objective values from optimization run...")

    recorder = ObjectiveRecorder()

    # Simulate Kursawe function evaluations
    import math
    def kursawe(x):
        f1 = sum(-10 * math.exp(-0.2 * math.sqrt(x[i]**2 + x[i+1]**2)) for i in range(len(x)-1))
        f2 = sum(abs(xi)**0.8 + 5 * math.sin(xi**3) for xi in x)
        return [f1, f2]

    # Generate test points
    test_points = []
    import random
    random.seed(42)

    for _ in range(100):
        x = [random.uniform(-5, 5) for _ in range(3)]
        objectives = kursawe(x)
        recorder.record(x, objectives)
        test_points.append(objectives)

    print(f"   Recorded {len(test_points)} objective evaluations")

    # Step 2: Extract Pareto front using direct C++ interface
    print("\nStep 2: Extracting Pareto front using direct C++ library...")

    start_time = time.time()

    # Convert to numpy array for Pareto extraction
    import numpy as np
    objectives_array = np.array(test_points)

    extractor = ParetoFrontExtractor()
    pareto_result = extractor.extract_pareto_front(objectives_array)

    pareto_time = time.time() - start_time

    print(f"   Pareto extraction completed in {pareto_time:.4f}s")
    print(f"   Original points: {len(test_points)}")
    print(f"   Pareto-optimal points: {len(pareto_result['pareto_front'])}")
    print(f"   Reduction ratio: {(len(test_points) - len(pareto_result['pareto_front'])) / len(test_points):.1%}")

    # Step 3: Calculate hypervolume using direct WFG library
    print("\nStep 3: Calculating hypervolume using direct WFG C++ library...")

    start_time = time.time()

    # Calculate reference point (slightly worse than nadir point)
    pareto_points = pareto_result['pareto_front'].tolist()

    if not pareto_points:
        print("   Warning: No Pareto points found!")
        return

    # Calculate nadir point (worst values in each objective)
    obj1_values = [p[0] for p in pareto_points]
    obj2_values = [p[1] for p in pareto_points]

    nadir_point = [max(obj1_values) + 1.0, max(obj2_values) + 1.0]
    reference_point = nadir_point  # Use nadir as reference

    print(f"   Reference point: [{reference_point[0]:.3f}, {reference_point[1]:.3f}]")

    # Calculate hypervolume using direct interface
    wfg = DirectWFGInterface()
    try:
        hypervolume = wfg.calculate_hypervolume(pareto_points, reference_point)
        hv_time = time.time() - start_time

        print(f"   Hypervolume calculation completed in {hv_time:.6f}s")
        print(f"   Hypervolume: {hypervolume:.6f}")

    finally:
        wfg.cleanup()

    # Step 4: Performance comparison summary
    print("\n=== Performance Summary ===")
    print(f"Pareto extraction time: {pareto_time:.4f}s")
    print(f"Hypervolume time: {hv_time:.6f}s")
    print(f"Total pipeline time: {pareto_time + hv_time:.4f}s")

    print("\n=== Results Summary ===")
    print(f"Original evaluation points: {len(test_points)}")
    print(f"Pareto-optimal points: {len(pareto_points)}")
    print(f"Final hypervolume: {hypervolume:.6f}")

    return {
        'original_points': len(test_points),
        'pareto_points': len(pareto_points),
        'hypervolume': hypervolume,
        'pareto_time': pareto_time,
        'hv_time': hv_time,
        'total_time': pareto_time + hv_time
    }

if __name__ == "__main__":
    try:
        result = demonstrate_direct_pipeline()
        print(f"\n✅ Direct C++ library integration successful!")
        print(f"   Achieved {result['original_points']} -> {result['pareto_points']} point reduction")
        print(f"   Hypervolume: {result['hypervolume']:.6f}")
        print(f"   Total processing time: {result['total_time']:.4f}s")

    except Exception as e:
        print(f"\n❌ Error in direct integration: {e}")
        import traceback
        traceback.print_exc()
