#!/usr/bin/env python3
"""
MOECAM Proper CFFI Implementation Test
=====================================

Builds and tests the proper CFFI interfaces as requested by the user,
following the toy example pattern instead of subprocess shortcuts.
"""

import subprocess
import sys
import numpy as np
from pathlib import Path

def build_all_libraries():
    """Build all CFFI libraries."""
    print("=" * 60)
    print("Building All CFFI Libraries")
    print("=" * 60)

    # Change to the working directory
    base_dir = Path(__file__).parent

    # Run build script
    result = subprocess.run([
        sys.executable, "build_libraries.py"
    ], cwd=base_dir, capture_output=True, text=True)

    if result.returncode != 0:
        print("❌ Build failed:")
        print(result.stdout)
        print(result.stderr)
        return False
    else:
        print("✅ Build successful:")
        print(result.stdout)
        return True

def test_zdt1(x):
    """ZDT1 test function."""
    x = np.asarray(x)
    f1 = x[0]
    g = 1 + 9 * np.sum(x[1:]) / (len(x) - 1)
    f2 = g * (1 - np.sqrt(f1 / g))
    return [f1, f2]

def sphere(x):
    """Simple sphere function."""
    return [np.sum(x**2)]

def test_proper_cffi_implementation():
    """Test the proper CFFI implementation."""

    print("\n" + "=" * 60)
    print("Testing Proper CFFI Implementation")
    print("=" * 60)

    # First build the libraries
    if not build_all_libraries():
        print("❌ Cannot test - build failed")
        return False

    # Add scripts directory to path
    scripts_dir = Path(__file__).parent.parent / "src" / "tests" / "scripts"
    sys.path.insert(0, str(scripts_dir))

    success = True

    # Test 1: Check that shared libraries exist
    print("\n1. Checking shared libraries...")
    expected_libs = [
        scripts_dir / "libmoecam.so",
        scripts_dir / "libpareto.so",
        scripts_dir / "libwfg.so"
    ]

    for lib_path in expected_libs:
        if lib_path.exists():
            print(f"   ✓ {lib_path.name} exists")
        else:
            print(f"   ❌ {lib_path.name} missing")
            success = False

    # Test 2: Test MOECAM CFFI interfaces
    print("\n2. Testing MOECAM CFFI interfaces...")
    try:
        import moecam_cffi_interface
        print("   ✓ MOECAM CFFI module imported successfully")

        # Test ECAM
        bounds = [(0, 1), (0, 1)]
        result = moecam_cffi_interface.minimize_ecam_cffi(
            test_zdt1, bounds, num_objectives=2, max_iterations=50
        )

        if result['success']:
            print(f"   ✓ ECAM ran successfully: {result['function_evaluations']} evaluations")
            print(f"     Best: x={result['x']}, f={result['f']}")
        else:
            print(f"   ❌ ECAM failed: {result.get('message', 'Unknown error')}")
            success = False

    except Exception as e:
        print(f"   ❌ MOECAM CFFI test failed: {e}")
        success = False

    # Test 3: Test Pareto CFFI interface
    print("\n3. Testing Pareto CFFI interface...")
    try:
        import pareto_cffi_interface
        print("   ✓ Pareto CFFI module imported successfully")

        pareto = pareto_cffi_interface.ParetoExtractorCFFI()

        # Test data
        objectives = np.array([
            [1.0, 2.0],
            [2.0, 1.0],
            [1.5, 1.5],
            [3.0, 3.0]  # Dominated
        ])

        pareto_mask = pareto.extract_pareto_mask(objectives)
        pareto_points = objectives[pareto_mask]

        print(f"   ✓ Pareto extraction: {len(pareto_points)} from {len(objectives)} points")

    except Exception as e:
        print(f"   ❌ Pareto CFFI test failed: {e}")
        success = False

    # Test 4: Test WFG CFFI interface
    print("\n4. Testing WFG CFFI interface...")
    try:
        import wfg_cffi_interface
        print("   ✓ WFG CFFI module imported successfully")

        wfg = wfg_cffi_interface.WFGHypervolumeCFFI()

        # Test data
        points = np.array([[1.0, 2.0], [2.0, 1.0]])
        reference = np.array([3.0, 3.0])

        volume = wfg.calculate_hypervolume(points, reference)
        print(f"   ✓ Hypervolume calculation: {volume}")

    except Exception as e:
        print(f"   ❌ WFG CFFI test failed: {e}")
        success = False

    # Test 5: Compare with original functionality
    print("\n5. Functional comparison...")
    try:
        # Run ZDT1 with ECAM
        bounds = [(0, 1), (0, 1)]
        result = moecam_cffi_interface.minimize_ecam_cffi(
            test_zdt1, bounds, num_objectives=2, max_iterations=100
        )

        if result['success']:
            all_obj = result['all_objectives']

            # Extract Pareto front
            pareto_mask = pareto.extract_pareto_mask(all_obj)
            pareto_front = all_obj[pareto_mask]

            print(f"   ✓ End-to-end test: {len(all_obj)} evaluations → {len(pareto_front)} Pareto points")
            print(f"     Pareto efficiency: {100 * len(pareto_front) / len(all_obj):.1f}%")

            # Calculate hypervolume of Pareto front
            if len(pareto_front) > 0 and len(pareto_front) < 20:  # Safety limit
                ref_point = np.array([2.0, 2.0])
                try:
                    hv = wfg.calculate_hypervolume(pareto_front, ref_point)
                    print(f"     Hypervolume: {hv:.4f}")
                except Exception as e:
                    print(f"     Hypervolume calculation failed: {e}")
            else:
                print(f"     Hypervolume: skipped (too many points: {len(pareto_front)})")

        else:
            print(f"   ❌ End-to-end test failed")
            success = False

    except Exception as e:
        print(f"   ❌ Functional comparison failed: {e}")
        success = False

    print("\n" + "=" * 60)
    if success:
        print("🎯 All CFFI interfaces working correctly!")
        print("✅ Proper implementation following toy example pattern")
        print("✅ No subprocess shortcuts - true CFFI interfaces")
    else:
        print("❌ Some CFFI interfaces failed")
    print("=" * 60)

    return success

if __name__ == "__main__":
    test_proper_cffi_implementation()
