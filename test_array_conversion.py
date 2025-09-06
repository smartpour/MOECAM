#!/usr/bin/env python
"""Test script for array conversion utilities and C++ interface robustness.

Tests various input types:
- Python lists
- NumPy arrays (different dtypes)
- Non-contiguous arrays
- Empty arrays
- Scalar values
"""

import numpy as np
import sys
import os

# Add package to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
from moecam.core.cffi_interface import CppOptimizer, convert_py_float_to_cffi, convert_py_int_to_cffi

def test_array_conversions():
    """Test the array conversion utilities."""
    print("=== Testing Array Conversion Utilities ===")

    # Test cases for float conversion
    test_cases_float = [
        ("Python list", [1.0, 2.0, 3.0]),
        ("NumPy float64", np.array([1.0, 2.0, 3.0], dtype=np.float64)),
        ("NumPy float32", np.array([1.0, 2.0, 3.0], dtype=np.float32)),
        ("NumPy int", np.array([1, 2, 3], dtype=np.int32)),
        ("Non-contiguous", np.array([1.0, 2.0, 3.0, 4.0, 5.0])[::2]),  # Every other element
        ("Empty array", []),
        ("None", None),
    ]

    for name, data in test_cases_float:
        try:
            np_arr, cffi_ptr = convert_py_float_to_cffi(data)
            print(f"✓ {name:15s}: shape={np_arr.shape}, dtype={np_arr.dtype}, contiguous={np_arr.flags.c_contiguous}")
        except Exception as e:
            print(f"✗ {name:15s}: ERROR: {e}")

    # Test cases for int conversion
    test_cases_int = [
        ("Python list", [1, 2, 3]),
        ("NumPy int32", np.array([1, 2, 3], dtype=np.int32)),
        ("NumPy int64", np.array([1, 2, 3], dtype=np.int64)),
        ("NumPy float", np.array([1.0, 2.0, 3.0], dtype=np.float64)),
    ]

    print("\n--- Integer Conversion Tests ---")
    for name, data in test_cases_int:
        try:
            np_arr, cffi_ptr = convert_py_int_to_cffi(data)
            print(f"✓ {name:15s}: shape={np_arr.shape}, dtype={np_arr.dtype}, contiguous={np_arr.flags.c_contiguous}")
        except Exception as e:
            print(f"✗ {name:15s}: ERROR: {e}")

def test_optimizer_with_various_inputs():
    """Test CppOptimizer with various input array types."""
    print("\n=== Testing CppOptimizer with Various Input Types ===")

    def simple_obj(x):
        return np.array([np.sum(x**2), np.sum((x-1)**2)])

    # Test different bound input types
    bound_test_cases = [
        ("Python lists", [-2, -2, -2], [2, 2, 2]),
        ("NumPy float64", np.array([-2, -2, -2], dtype=np.float64), np.array([2, 2, 2], dtype=np.float64)),
        ("NumPy float32", np.array([-2, -2, -2], dtype=np.float32), np.array([2, 2, 2], dtype=np.float32)),
        ("Mixed types", [-2.0, -2, -2], np.array([2, 2, 2])),
        ("Non-contiguous", np.array([-3, -2, -1, 0, 1, 2, 3])[1:4], np.array([0, 1, 2, 3, 4, 5, 6])[4:7]),
    ]

    for name, lb, ub in bound_test_cases:
        try:
            opt = CppOptimizer(3, 2, lb, ub, simple_obj)
            result = opt.evaluate([0, 0, 0])
            print(f"✓ {name:15s}: Created optimizer, eval([0,0,0]) = {result}")
        except Exception as e:
            print(f"✗ {name:15s}: ERROR: {e}")

    # Test evaluation with different input types
    opt = CppOptimizer(3, 2, [-2, -2, -2], [2, 2, 2], simple_obj)

    eval_test_cases = [
        ("Python list", [1.0, 0.5, -0.5]),
        ("NumPy float64", np.array([1.0, 0.5, -0.5], dtype=np.float64)),
        ("NumPy float32", np.array([1.0, 0.5, -0.5], dtype=np.float32)),
        ("NumPy int", np.array([1, 0, -1], dtype=np.int32)),
        ("Non-contiguous", np.array([0, 1.0, 2, 0.5, 4, -0.5, 6])[1::2]),
    ]

    print("\n--- Evaluation Input Tests ---")
    for name, x in eval_test_cases:
        try:
            result = opt.evaluate(x)
            print(f"✓ {name:15s}: eval({x if len(str(x)) < 20 else '...'}) = {result}")
        except Exception as e:
            print(f"✗ {name:15s}: ERROR: {e}")

def test_crossover_mutate_inputs():
    """Test crossover and mutation with various input types."""
    print("\n=== Testing Crossover/Mutation with Various Inputs ===")

    def simple_obj(x):
        return np.array([np.sum(x**2), np.sum((x-1)**2)])

    opt = CppOptimizer(3, 2, [-2, -2, -2], [2, 2, 2], simple_obj)

    # Test crossover inputs
    parent_test_cases = [
        ("Lists", [1.0, 0.0, -1.0], [0.5, 1.0, 0.5]),
        ("NumPy arrays", np.array([1.0, 0.0, -1.0]), np.array([0.5, 1.0, 0.5])),
        ("Mixed", [1.0, 0.0, -1.0], np.array([0.5, 1.0, 0.5])),
        ("Non-contiguous",
         np.array([0, 1.0, 2, 0.0, 4, -1.0])[1::2],
         np.array([0, 0.5, 2, 1.0, 4, 0.5])[1::2]),
    ]

    print("\n--- Crossover Input Tests ---")
    for name, p1, p2 in parent_test_cases:
        try:
            c1, c2 = opt.crossover(p1, p2)
            print(f"✓ {name:15s}: crossover successful, child1 shape={c1.shape}")
        except Exception as e:
            print(f"✗ {name:15s}: ERROR: {e}")

    print("\n--- Mutation Input Tests ---")
    for name, individual, _ in parent_test_cases[:3]:  # Skip non-contiguous for mutation
        try:
            mutated = opt.mutate(individual, 0.1)
            print(f"✓ {name:15s}: mutation successful, result shape={mutated.shape}")
        except Exception as e:
            print(f"✗ {name:15s}: ERROR: {e}")

if __name__ == "__main__":
    test_array_conversions()
    test_optimizer_with_various_inputs()
    test_crossover_mutate_inputs()
    print("\n=== Array Conversion Tests Complete ===")
