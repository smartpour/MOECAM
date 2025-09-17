#!/usr/bin/env python3
"""
Simple CFFI Algorithm Test
==========================

Test with a very simple objective function to isolate callback issues.
"""

import sys
from pathlib import Path

scripts_dir = Path(__file__).parent.parent / "src" / "tests" / "scripts"
sys.path.insert(0, str(scripts_dir))

def simple_objective(x):
    """Simple objective function without numpy dependencies."""
    # Convert to plain Python list if needed
    if hasattr(x, '__iter__'):
        x_list = list(x)
    else:
        x_list = [x]

    # Simple sphere function for 2 objectives
    f1 = sum(xi**2 for xi in x_list)
    f2 = sum((xi - 1)**2 for xi in x_list)

    return [f1, f2]

print("Testing simple CFFI algorithm call...")

try:
    import moecam_cffi_interface

    # Test with very simple bounds and low iterations
    bounds = [(0, 1), (0, 1)]

    print("Creating MOECAM instance...")
    moecam = moecam_cffi_interface.MOECAMAlgorithmsCFFI()

    print("Running ECAM with 5 iterations...")
    result = moecam.minimize_ecam(
        objective_func=simple_objective,
        bounds=bounds,
        num_objectives=2,
        max_iterations=5
    )

    print(f"✓ ECAM completed successfully!")
    print(f"  Success: {result['success']}")
    print(f"  Evaluations: {result['function_evaluations']}")
    print(f"  Best solution: {result['x']}")
    print(f"  Best objectives: {result['f']}")

except Exception as e:
    print(f"❌ Test failed: {e}")
    import traceback
    traceback.print_exc()
