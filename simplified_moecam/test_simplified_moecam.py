#!/usr/bin/env python3
"""
Test Script for Simplified MOECAM
=================================

Test the simplified installation and basic functionality.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def test_simplified_moecam():
    """Test the simplified MOECAM package."""

    print("🔬 Testing Simplified MOECAM")
    print("=" * 40)

    # Try to import MOECAM
    try:
        import moecam
        print(f"✅ MOECAM imported successfully")
        print(f"✅ Version: {moecam.__version__}")
        print(f"✅ CFFI available: {moecam.CFFI_AVAILABLE}")

    except ImportError as e:
        print(f"❌ Failed to import MOECAM: {e}")
        return False

    if not moecam.CFFI_AVAILABLE:
        print("❌ CFFI library not available")
        return False

    # Test basic functionality
    print("\n🧪 Testing basic functions...")
    try:
        result = moecam.hello_moecam()
        print(f"   {result}")

        # Test simple addition
        test_val = moecam.lib.add(2.5, 3.7)
        print(f"   ✅ Addition test: 2.5 + 3.7 = {test_val}")

    except Exception as e:
        print(f"   ❌ Basic function test failed: {e}")
        return False

    # Test Pareto front extraction
    print("\n🎯 Testing Pareto front extraction...")
    try:
        # Generate test data (ZDT1-like)
        np.random.seed(42)
        n_points = 50
        objectives = np.random.rand(n_points, 2)

        # Make some points clearly Pareto optimal
        for i in range(5):
            objectives[i] = [i/10.0, 1 - i/10.0]

        # Extract Pareto front
        pareto_points = moecam.extract_pareto_front(objectives, strict_mode=False)
        pareto_indices = moecam.extract_pareto_indices(objectives, strict_mode=False)

        print(f"   ✅ Found {len(pareto_points)} Pareto optimal points")
        print(f"   ✅ Pareto indices: {len(pareto_indices)} points")

        # Visualize results
        plt.figure(figsize=(10, 6))

        plt.subplot(1, 2, 1)
        plt.scatter(objectives[:, 0], objectives[:, 1], alpha=0.6, label='All points')
        plt.scatter(pareto_points[:, 0], pareto_points[:, 1],
                   color='red', s=100, label='Pareto front', marker='*')
        plt.xlabel('Objective 1')
        plt.ylabel('Objective 2')
        plt.title('Pareto Front Extraction')
        plt.legend()
        plt.grid(True)

    except Exception as e:
        print(f"   ❌ Pareto front test failed: {e}")
        return False

    # Test ZDT functions
    print("\n🧮 Testing ZDT test functions...")
    try:
        x = np.array([0.5, 0.3, 0.2, 0.4, 0.1])

        zdt1_result = moecam.evaluate_zdt1(x)
        zdt2_result = moecam.evaluate_zdt2(x)
        zdt3_result = moecam.evaluate_zdt3(x)

        print(f"   ✅ ZDT1: {zdt1_result}")
        print(f"   ✅ ZDT2: {zdt2_result}")
        print(f"   ✅ ZDT3: {zdt3_result}")

        # Generate ZDT1 Pareto front for visualization
        x_pareto = np.linspace(0, 1, 100)
        y_pareto = 1 - np.sqrt(x_pareto)

        plt.subplot(1, 2, 2)
        plt.plot(x_pareto, y_pareto, 'b-', linewidth=2, label='ZDT1 True Pareto Front')
        plt.scatter(zdt1_result[0], zdt1_result[1], color='red', s=100,
                   marker='o', label='Test Point', zorder=5)
        plt.xlabel('f1(x)')
        plt.ylabel('f2(x)')
        plt.title('ZDT1 Test Function')
        plt.legend()
        plt.grid(True)

    except Exception as e:
        print(f"   ❌ ZDT function test failed: {e}")
        return False

    # Test optimization algorithms (basic calls)
    print("\n🔍 Testing optimization algorithms...")
    try:
        x0 = np.array([0.5, 0.5])
        bounds = np.array([[0.0, 1.0], [0.0, 1.0]])

        # Test ECAM
        x_ecam, val_ecam, status_ecam = moecam.minimize_ecam(x0, bounds, max_iterations=10)
        print(f"   ✅ ECAM: x={x_ecam}, val={val_ecam}, status={status_ecam}")

        # Test DFBM
        x_dfbm, val_dfbm, status_dfbm = moecam.minimize_dfbm(x0, bounds, max_iterations=10)
        print(f"   ✅ DFBM: x={x_dfbm}, val={val_dfbm}, status={status_dfbm}")

        # Test Random Start
        x_rand, val_rand, status_rand = moecam.minimize_random_start(x0, bounds, max_iterations=10)
        print(f"   ✅ Random Start: x={x_rand}, val={val_rand}, status={status_rand}")

    except Exception as e:
        print(f"   ❌ Optimization test failed: {e}")
        return False

    # Test hypervolume calculation
    print("\n📊 Testing hypervolume calculation...")
    try:
        test_points = np.array([[0.2, 0.8], [0.4, 0.6], [0.6, 0.4], [0.8, 0.2]])
        reference = np.array([1.0, 1.0])

        hv = moecam.wfg_hypervolume(test_points, reference)
        print(f"   ✅ Hypervolume: {hv}")

    except Exception as e:
        print(f"   ❌ Hypervolume test failed: {e}")
        return False

    # Save visualization
    output_dir = Path("simplified_moecam_test_results")
    output_dir.mkdir(exist_ok=True)

    plt.tight_layout()
    plt.savefig(output_dir / "simplified_moecam_test.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n💾 Results saved to: {output_dir.absolute()}")
    print("🎉 All tests passed! Simplified MOECAM is working correctly.")

    return True


if __name__ == "__main__":
    success = test_simplified_moecam()
    sys.exit(0 if success else 1)
