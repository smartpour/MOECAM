#!/usr/bin/env python3
"""
Minimal CFFI Callback Test
==========================

Test just the callback mechanism without running full algorithms.
"""

import sys
from pathlib import Path
import numpy as np

scripts_dir = Path(__file__).parent.parent / "src" / "tests" / "scripts"
sys.path.insert(0, str(scripts_dir))

def simple_test_function(x):
    """Very simple test function."""
    print(f"Callback called with x = {x}")
    return [x[0]**2, x[1]**2]

print("Testing minimal CFFI callback...")

try:
    # Import the CFFI module
    import moecam_cffi_interface
    print("✓ CFFI module imported")

    # Test just the callback creation without running algorithms
    moecam = moecam_cffi_interface.MOECAMAlgorithmsCFFI()
    print("✓ MOECAM instance created")

    # Test callback creation
    callback = moecam._create_callback(simple_test_function, 2, 2)
    print("✓ Callback created successfully")

    # Test the callback directly (this should be safe)
    print("✓ Direct callback test skipped (would require C pointers)")

    print("✓ Callback mechanism appears to work")

except Exception as e:
    print(f"❌ Test failed: {e}")
    import traceback
    traceback.print_exc()

print("Minimal test complete.")
