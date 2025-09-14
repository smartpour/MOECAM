"""
Test the direct WFG interface with manual memory management.
"""

import ctypes
import os

# Load the shared library
script_dir = os.path.dirname(os.path.abspath(__file__))
library_path = os.path.join(script_dir, 'libwfg_simple.so')

lib = ctypes.CDLL(library_path)

# Define function signatures
lib.wfg_simple_init.argtypes = [ctypes.c_int, ctypes.c_int]
lib.wfg_simple_init.restype = ctypes.c_int

lib.wfg_simple_calculate.argtypes = [
    ctypes.POINTER(ctypes.c_double),  # pts
    ctypes.POINTER(ctypes.c_double),  # nadir
    ctypes.c_int,                     # nPoints
    ctypes.c_int                      # dimensions
]
lib.wfg_simple_calculate.restype = ctypes.c_double

lib.wfg_simple_cleanup.argtypes = []
lib.wfg_simple_cleanup.restype = None

print("Testing WFG Direct C Interface")

# Test 1: Simple 2D case
print("\n=== Test 1: Simple 2D Points ===")

# Initialize for 2D
result = lib.wfg_simple_init(2, 100)
print(f"Init result: {result}")

# Points: [(2,8), (4,6), (6,4), (8,2)]
points_data = [2.0, 8.0, 4.0, 6.0, 6.0, 4.0, 8.0, 2.0]
points_c = (ctypes.c_double * len(points_data))(*points_data)

# Reference point: (10, 10)
ref_data = [10.0, 10.0]
ref_c = (ctypes.c_double * len(ref_data))(*ref_data)

hv = lib.wfg_simple_calculate(points_c, ref_c, 4, 2)
print(f"2D Hypervolume: {hv}")

lib.wfg_simple_cleanup()

# Test 2: Different data
print("\n=== Test 2: Kursawe-like Data ===")

result = lib.wfg_simple_init(2, 100)
print(f"Init result: {result}")

# Kursawe-like points
points_data = [-19.0, 2.0, -18.5, 3.0, -17.8, 4.5, -16.2, 7.0, -14.5, 10.5]
points_c = (ctypes.c_double * len(points_data))(*points_data)

ref_data = [0.0, 15.0]
ref_c = (ctypes.c_double * len(ref_data))(*ref_data)

hv = lib.wfg_simple_calculate(points_c, ref_c, 5, 2)
print(f"Kursawe-like Hypervolume: {hv}")

lib.wfg_simple_cleanup()

print("\nAll tests completed!")
