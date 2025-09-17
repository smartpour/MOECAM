#!/usr/bin/env python3
"""
Simple MOECAM Installation Test
==============================

Quick verification test for MOECAM CFFI installation.
This test verifies all core functionality is working properly.
"""

import sys
import numpy as np
from pathlib import Path

# Add source to path for development testing
sys.path.insert(0, str(Path(__file__).parent / 'src'))

def test_zdt1(x):
    """Simple ZDT1 test function."""
    x = np.asarray(x)
    f1 = x[0]
    g = 1 + 9 * sum(x[1:]) / (len(x) - 1) if len(x) > 1 else 1
    f2 = g * (1 - (f1 / g)**0.5)
    return [f1, f2]

def run_installation_test():
    """Run basic installation verification test."""
    print("🔧 MOECAM Installation Test")
    print("=" * 30)

    try:
        # Try to import the CFFI interface
        print("📦 Testing CFFI import...")
        from moecam.cffi_interface import (
            minimize_ecam_cffi,
            extract_pareto_front_cffi,
            calculate_hypervolume_cffi
        )
        print("✅ CFFI interface imported successfully")

        # Test ECAM optimization
        print("\n🧮 Testing ECAM optimization...")
        bounds = [(0, 1), (0, 1)]

        result = minimize_ecam_cffi(
            test_zdt1,
            bounds,
            num_objectives=2,
            max_iterations=20
        )

        if result['success']:
            print(f"✅ ECAM optimization successful!")
            print(f"   Function evaluations: {result['function_evaluations']}")
            print(f"   Best objectives: {result['best_objectives']}")
        else:
            print("❌ ECAM optimization failed")
            return False

        # Test Pareto front extraction
        print("\n📊 Testing Pareto front extraction...")
        all_objectives = result['all_objectives']

        pareto_points, pareto_indices = extract_pareto_front_cffi(
            all_objectives,
            strict_mode=False
        )

        print(f"✅ Pareto front extraction successful!")
        print(f"   Total points: {len(all_objectives)}")
        print(f"   Pareto points: {len(pareto_points)}")

        # Test hypervolume calculation
        print("\n📏 Testing hypervolume calculation...")
        reference_point = np.array([1.1, 1.1])

        hypervolume = calculate_hypervolume_cffi(pareto_points, reference_point)

        print(f"✅ Hypervolume calculation successful!")
        print(f"   Hypervolume: {hypervolume:.6f}")

        print("\n🎉 All tests passed! MOECAM is properly installed.")
        return True

    except ImportError as e:
        print(f"❌ Import error: {e}")
        print("\n💡 To fix this issue:")
        print("1. Navigate to the project directory")
        print("2. Activate your virtual environment")
        print("3. Run: pip install -e .")
        return False

    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = run_installation_test()
    sys.exit(0 if success else 1)
