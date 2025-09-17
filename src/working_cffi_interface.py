"""
Working MOECAM CFFI Interface
============================

Simple Python interface to the working CFFI extension.
"""

import numpy as np
from typing import List, Tuple, Dict, Any
import warnings

try:
    from moecam._moecam_cffi import ffi, lib
    CFFI_AVAILABLE = True
except ImportError:
    ffi = None
    lib = None
    CFFI_AVAILABLE = False
    warnings.warn("MOECAM CFFI library not found. Run build_working.py to build it.")


def test_basic_functionality():
    """Test basic CFFI functionality."""
    if not CFFI_AVAILABLE:
        raise RuntimeError("CFFI library not available")

    # Test toy functions
    result = lib.add(2.5, 3.7)
    print(f"✅ add(2.5, 3.7) = {result}")

    lib.print_message(b"Hello from CFFI!")
    return True


def extract_pareto_front_cffi(points: np.ndarray, strict_mode: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract Pareto front from a set of objective vectors.

    Args:
        points: Array of shape (n_points, n_objectives)
        strict_mode: Use strict Pareto dominance if True, weak if False

    Returns:
        Tuple of (pareto_points, pareto_indices)
    """
    if not CFFI_AVAILABLE:
        raise RuntimeError("CFFI library not available")

    points = np.asarray(points, dtype=np.float64)
    if points.ndim != 2:
        raise ValueError("Points must be a 2D array")

    n_points, n_objectives = points.shape
    max_pareto = n_points  # At most all points can be Pareto optimal

    # Flatten points for C interface
    points_flat = ffi.new("double[]", points.flatten().tolist())
    pareto_points_flat = ffi.new("double[]", max_pareto * n_objectives)
    pareto_indices = ffi.new("int[]", max_pareto)

    # Extract Pareto front
    n_pareto = lib.extract_pareto_front(points_flat, n_points, n_objectives,
                                      pareto_points_flat, max_pareto, int(strict_mode))

    if n_pareto < 0:
        raise RuntimeError("Pareto front extraction failed")

    # Extract results
    pareto_points = np.zeros((n_pareto, n_objectives))
    for i in range(n_pareto):
        for j in range(n_objectives):
            pareto_points[i, j] = pareto_points_flat[i * n_objectives + j]

    # Get indices too
    n_indices = lib.extract_pareto_indices(points_flat, n_points, n_objectives,
                                         pareto_indices, max_pareto, int(strict_mode))

    indices = np.array([pareto_indices[i] for i in range(n_indices)])

    return pareto_points, indices


def generate_test_data(func_name: str = "ZDT1", n_points: int = 100, n_dim: int = 2) -> np.ndarray:
    """Generate test data for multi-objective functions."""
    np.random.seed(42)  # For reproducible results

    if func_name == "ZDT1":
        # Random points in [0,1]^n_dim
        x = np.random.uniform(0, 1, (n_points, n_dim))
        objectives = np.zeros((n_points, 2))

        for i in range(n_points):
            f1 = x[i, 0]
            g = 1 + 9 * np.sum(x[i, 1:]) / (n_dim - 1) if n_dim > 1 else 1
            f2 = g * (1 - np.sqrt(f1 / g))
            objectives[i] = [f1, f2]

        return objectives

    elif func_name == "ZDT2":
        # Random points in [0,1]^n_dim
        x = np.random.uniform(0, 1, (n_points, n_dim))
        objectives = np.zeros((n_points, 2))

        for i in range(n_points):
            f1 = x[i, 0]
            g = 1 + 9 * np.sum(x[i, 1:]) / (n_dim - 1) if n_dim > 1 else 1
            f2 = g * (1 - (f1 / g)**2)
            objectives[i] = [f1, f2]

        return objectives

    else:
        # Random test data
        return np.random.uniform(0, 10, (n_points, 2))


def run_working_test():
    """Run a working test of the CFFI implementation."""
    print("🔧 MOECAM Working CFFI Test")
    print("=" * 30)

    if not CFFI_AVAILABLE:
        print("❌ CFFI library not available. Run build_working.py first.")
        return False

    try:
        # Test basic functionality
        print("📦 Testing basic functions...")
        test_basic_functionality()

        # Generate test data
        print("\n🧮 Generating test data...")
        test_points = generate_test_data("ZDT1", n_points=50, n_dim=3)
        print(f"   Generated {len(test_points)} points with {test_points.shape[1]} objectives")

        # Test Pareto front extraction
        print("\n📊 Testing Pareto front extraction...")
        pareto_points, pareto_indices = extract_pareto_front_cffi(test_points, strict_mode=False)

        print(f"✅ Pareto front extraction successful!")
        print(f"   Total points: {len(test_points)}")
        print(f"   Pareto points: {len(pareto_points)}")
        print(f"   Pareto ratio: {len(pareto_points)/len(test_points):.2%}")

        # Save results
        print("\n💾 Saving results...")
        output_dir = "../working_test_results"
        import os
        os.makedirs(output_dir, exist_ok=True)

        np.savetxt(f"{output_dir}/all_points.txt", test_points,
                   header="All test points (ZDT1)", fmt="%.6f")
        np.savetxt(f"{output_dir}/pareto_points.txt", pareto_points,
                   header="Pareto front points", fmt="%.6f")
        np.savetxt(f"{output_dir}/pareto_indices.txt", pareto_indices,
                   header="Pareto front indices", fmt="%d")

        # Create basic plot if matplotlib available
        try:
            import matplotlib.pyplot as plt

            plt.figure(figsize=(10, 6))
            plt.scatter(test_points[:, 0], test_points[:, 1],
                       alpha=0.6, s=20, c='lightblue', label='All Points')
            plt.scatter(pareto_points[:, 0], pareto_points[:, 1],
                       alpha=0.8, s=50, c='red', label='Pareto Front')

            # Add true Pareto front for ZDT1
            f1_true = np.linspace(0, 1, 100)
            f2_true = 1 - np.sqrt(f1_true)
            plt.plot(f1_true, f2_true, 'g-', linewidth=2, label='True Pareto Front')

            plt.xlabel('Objective 1')
            plt.ylabel('Objective 2')
            plt.title('MOECAM Working Test - ZDT1 Pareto Front')
            plt.legend()
            plt.grid(True, alpha=0.3)

            plt.savefig(f"{output_dir}/pareto_test.png", dpi=300, bbox_inches='tight')
            plt.close()

            print(f"📈 Plot saved to {output_dir}/pareto_test.png")

        except ImportError:
            print("📈 Matplotlib not available, skipping plot")

        print(f"\n📁 Results saved to: {output_dir}/")
        print("\n🎉 All tests passed! Working CFFI implementation is functional.")
        return True

    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = run_working_test()
    exit(0 if success else 1)
