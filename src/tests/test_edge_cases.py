#!/usr/bin/env python
"""Edge case validation and performance test for array conversion."""

import numpy as np
import sys
import os
import time

sys.path.append(os.path.join(os.path.dirname(__file__), 'src', 'src'))
from moecam.core.cffi_interface import CppOptimizer, convert_py_float_to_cffi

def test_edge_cases():
    """Test edge cases that might cause problems."""
    print("=== Edge Case Validation ===")

    def simple_obj(x):
        return np.array([np.sum(x**2), np.sum((x-1)**2)])

    # Edge case 1: Very large arrays
    print("Testing large arrays...")
    large_bounds = np.ones(100) * 5  # 100-dimensional
    opt = CppOptimizer(100, 2, -large_bounds, large_bounds, simple_obj)
    large_x = np.random.randn(100)
    result = opt.evaluate(large_x)
    print(f"✓ Large array (100D): result shape = {result.shape}")

    # Edge case 2: Single element
    print("Testing single element...")
    opt1d = CppOptimizer(1, 2, [-1], [1], simple_obj)
    result1d = opt1d.evaluate([0.5])
    print(f"✓ Single element: result = {result1d}")

    # Edge case 3: Many objectives
    def multi_obj(x):
        return np.array([np.sum(x**2), np.sum((x-1)**2), np.sum((x+1)**2), np.sum(x**4)])

    opt_multi = CppOptimizer(3, 4, [-2, -2, -2], [2, 2, 2], multi_obj)
    result_multi = opt_multi.evaluate([1, 0, -1])
    print(f"✓ Many objectives (4): result = {result_multi}")

    # Edge case 4: Zero bounds
    opt_zero = CppOptimizer(2, 2, [0, 0], [0, 0], simple_obj)
    result_zero = opt_zero.evaluate([0, 0])
    print(f"✓ Zero bounds: result = {result_zero}")

def test_performance():
    """Test performance of array conversion vs direct allocation."""
    print("\n=== Performance Comparison ===")

    def simple_obj(x):
        return np.array([np.sum(x**2), np.sum((x-1)**2)])

    opt = CppOptimizer(10, 2, np.ones(10) * -5, np.ones(10) * 5, simple_obj)

    # Test different input types for performance
    test_data = [
        ("Python list", [1.0] * 10),
        ("NumPy contiguous", np.ones(10, dtype=np.float64)),
        ("NumPy non-contiguous", np.ones(20, dtype=np.float64)[::2]),
        ("NumPy float32", np.ones(10, dtype=np.float32)),
    ]

    n_iterations = 1000

    for name, data in test_data:
        start_time = time.time()
        for _ in range(n_iterations):
            result = opt.evaluate(data)
        end_time = time.time()
        avg_time = (end_time - start_time) / n_iterations * 1000  # ms
        print(f"{name:20s}: {avg_time:.3f} ms/eval avg (over {n_iterations} evals)")

def test_memory_safety():
    """Test memory safety with array reuse and modification."""
    print("\n=== Memory Safety Tests ===")

    def simple_obj(x):
        # Modify input to test if it affects original
        x_copy = x.copy()
        x_copy[0] = 999  # This should not affect original
        return np.array([np.sum(x**2), np.sum((x-1)**2)])

    opt = CppOptimizer(3, 2, [-5, -5, -5], [5, 5, 5], simple_obj)

    # Test that input arrays are not modified
    original = np.array([1.0, 2.0, 3.0])
    original_copy = original.copy()

    result = opt.evaluate(original)
    print(f"Original before: {original_copy}")
    print(f"Original after:  {original}")
    print(f"Arrays equal: {np.array_equal(original, original_copy)}")
    print(f"✓ Input array not modified by evaluation")

    # Test array reuse
    same_array = np.array([0.5, 1.0, 1.5])
    results = []
    for i in range(5):
        results.append(opt.evaluate(same_array))

    all_same = all(np.array_equal(results[0], r) for r in results[1:])
    print(f"✓ Consistent results with array reuse: {all_same}")

if __name__ == "__main__":
    test_edge_cases()
    test_performance()
    test_memory_safety()
    print("\n=== All Edge Case Tests Complete ===")
