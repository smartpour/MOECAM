#!/usr/bin/env python3
"""
Test Simple Callback Wrapper
============================

Test the callback mechanism with the simplest possible C++ wrapper.
"""

import sys
from pathlib import Path
import numpy as np
import cffi

scripts_dir = Path(__file__).parent.parent / "src" / "tests" / "scripts"
sys.path.insert(0, str(scripts_dir))

# Setup CFFI for simple test
ffi = cffi.FFI()
ffi.cdef("""
    typedef void (*objective_function_t)(double* x, int n_dim, double* f, int n_obj);
    int test_callback(objective_function_t callback, double* test_x, int n_dim, double* result_f, int n_obj);
""")

# Load simple test library
lib_path = scripts_dir / 'libsimpletest.so'
try:
    lib = ffi.dlopen(str(lib_path))
    print("✓ Simple test library loaded")
except Exception as e:
    print(f"❌ Failed to load simple test library: {e}")
    sys.exit(1)

def test_objective(x):
    """Simple test objective function."""
    print(f"Python: test_objective called with x = {list(x)}")
    return [x[0]**2, x[1]**2]

print("Testing simple callback mechanism...")

try:
    # Create callback
    @ffi.callback("void(double*, int, double*, int)")
    def c_callback(x_ptr, n_dim, f_ptr, n_obj):
        print(f"Python: Callback invoked with n_dim={n_dim}, n_obj={n_obj}")

        # Extract input
        x_arr = np.frombuffer(ffi.buffer(x_ptr, n_dim * 8), dtype=np.float64).copy()
        print(f"Python: Input array = {x_arr}")

        # Call Python function
        result = test_objective(x_arr)
        print(f"Python: Result = {result}")

        # Write back
        for i in range(n_obj):
            f_ptr[i] = float(result[i])

        print("Python: Callback completed")

    print("✓ Callback created")

    # Test data
    test_x = np.array([0.5, 0.3], dtype=np.float64)
    result_f = np.zeros(2, dtype=np.float64)

    # Convert to CFFI pointers
    x_ptr = ffi.cast("double *", test_x.ctypes.data)
    f_ptr = ffi.cast("double *", result_f.ctypes.data)

    print("✓ Test data prepared")

    # Call C++ function
    print("Calling C++ test function...")
    status = lib.test_callback(c_callback, x_ptr, 2, f_ptr, 2)

    if status == 0:
        print("✓ C++ test function completed successfully")
        print(f"✓ Result: {result_f}")
    else:
        print(f"❌ C++ test function failed with status {status}")

except Exception as e:
    print(f"❌ Test failed: {e}")
    import traceback
    traceback.print_exc()

print("Simple callback test complete.")
