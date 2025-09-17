#!/usr/bin/env python3
"""
Simple CFFI Test
================

Test the basic CFFI loading without running algorithms.
"""

import sys
from pathlib import Path

# Add scripts directory
scripts_dir = Path(__file__).parent.parent / "src" / "tests" / "scripts"
sys.path.insert(0, str(scripts_dir))

print("Testing CFFI library loading...")

try:
    import moecam_cffi_interface
    print("✓ MOECAM CFFI module imported successfully")

    # Test if library is loaded
    if moecam_cffi_interface.lib is not None:
        print("✓ MOECAM shared library loaded successfully")
    else:
        print("❌ MOECAM shared library not loaded")

except Exception as e:
    print(f"❌ CFFI loading failed: {e}")
    import traceback
    traceback.print_exc()

print("\nTesting other interfaces...")

try:
    import pareto_cffi_interface
    print("✓ Pareto CFFI module imported successfully")
except Exception as e:
    print(f"❌ Pareto CFFI loading failed: {e}")

try:
    import wfg_cffi_interface
    print("✓ WFG CFFI module imported successfully")
except Exception as e:
    print(f"❌ WFG CFFI loading failed: {e}")

print("Basic CFFI loading test complete.")
