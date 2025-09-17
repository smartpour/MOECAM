#!/usr/bin/env python3
"""
Test Proper CFFI Interfaces
===========================

Tests the proper CFFI interfaces that follow the toy example pattern,
replacing the subprocess shortcuts.
"""

import sys
import numpy as np
from pathlib import Path

# Add the scripts directory to the path
scripts_dir = Path(__file__).parent.parent / "src" / "tests" / "scripts"
sys.path.insert(0, str(scripts_dir))

def test_zdt1(x):
    """ZDT1 test function."""
    x = np.asarray(x)
    f1 = x[0]
    g = 1 + 9 * np.sum(x[1:]) / (len(x) - 1)
    f2 = g * (1 - np.sqrt(f1 / g))
    return [f1, f2]

def sphere(x):
    """Simple sphere function for single-objective testing."""
    return [np.sum(x**2)]

def test_cffi_interfaces():
    """Test all CFFI interfaces."""

    print("=" * 60)
    print("Testing Proper CFFI Interfaces")
    print("=" * 60)

    # Test 1: MOECAM Algorithms CFFI
    print("\n1. Testing MOECAM ECAM Algorithm CFFI...")
    try:
        from moecam_cffi_interface import minimize_ecam_cffi

        bounds = [(-1, 1), (-1, 1)]
        result = minimize_ecam_cffi(test_zdt1, bounds, num_objectives=2, max_iterations=50)

        print(f"   ✓ ECAM completed with {result['function_evaluations']} evaluations")
        print(f"   ✓ Best solution: {result['x']}")
        print(f"   ✓ Best objectives: {result['f']}")
        print(f"   ✓ Evaluated {len(result['all_solutions'])} points")

    except ImportError as e:
        print(f"   ❌ MOECAM CFFI interface not available: {e}")
    except Exception as e:
        print(f"   ❌ MOECAM CFFI test failed: {e}")

    # Test 2: WFG Hypervolume CFFI
    print("\n2. Testing WFG Hypervolume CFFI...")
    try:
        from wfg_cffi_interface import WFGHypervolumeCFFI

        wfg = WFGHypervolumeCFFI()

        # Test data
        points = np.array([[1.0, 2.0], [2.0, 1.0], [1.5, 1.5]])
        reference = np.array([3.0, 3.0])

        volume = wfg.calculate_hypervolume(points, reference)
        print(f"   ✓ Hypervolume calculated: {volume}")

    except ImportError as e:
        print(f"   ❌ WFG CFFI interface not available: {e}")
    except Exception as e:
        print(f"   ❌ WFG CFFI test failed: {e}")

    # Test 3: Pareto Front CFFI
    print("\n3. Testing Pareto Front Extraction CFFI...")
    try:
        from pareto_cffi_interface import ParetoFrontCFFI

        pareto = ParetoFrontCFFI()

        # Test data
        objectives = np.array([
            [1.0, 2.0],
            [2.0, 1.0],
            [1.5, 1.5],
            [3.0, 3.0]  # Dominated
        ])

        pareto_mask = pareto.extract_pareto_front(objectives)
        pareto_points = objectives[pareto_mask]

        print(f"   ✓ Found {len(pareto_points)} Pareto points from {len(objectives)} total")
        print(f"   ✓ Pareto points: {pareto_points}")

    except ImportError as e:
        print(f"   ❌ Pareto CFFI interface not available: {e}")
    except Exception as e:
        print(f"   ❌ Pareto CFFI test failed: {e}")

    # Test 4: Random Start Algorithm CFFI
    print("\n4. Testing Random Start Algorithm CFFI...")
    try:
        from moecam_cffi_interface import minimize_random_start_cffi

        bounds = [(-2, 2), (-2, 2)]
        result = minimize_random_start_cffi(test_zdt1, bounds, num_objectives=2, max_iterations=30)

        print(f"   ✓ Random Start completed with {result['function_evaluations']} evaluations")
        print(f"   ✓ Best solution: {result['x']}")
        print(f"   ✓ Best objectives: {result['f']}")

    except ImportError as e:
        print(f"   ❌ Random Start CFFI interface not available: {e}")
    except Exception as e:
        print(f"   ❌ Random Start CFFI test failed: {e}")

    # Test 5: DIRECT Algorithm CFFI
    print("\n5. Testing DIRECT Algorithm CFFI...")
    try:
        from moecam_cffi_interface import direct_optimize_cffi

        bounds = [(-2, 2), (-2, 2)]
        result = direct_optimize_cffi(sphere, bounds, max_iterations=25)

        print(f"   ✓ DIRECT completed with {result['function_evaluations']} evaluations")
        print(f"   ✓ Best solution: {result['x']}")
        print(f"   ✓ Best objective: {result['fun']}")

    except ImportError as e:
        print(f"   ❌ DIRECT CFFI interface not available: {e}")
    except Exception as e:
        print(f"   ❌ DIRECT CFFI test failed: {e}")

    print("\n" + "=" * 60)
    print("CFFI Interface Testing Complete")
    print("=" * 60)

if __name__ == "__main__":
    test_cffi_interfaces()
