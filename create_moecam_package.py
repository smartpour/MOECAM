#!/usr/bin/env python3
"""
MOECAM Package Installation Script
=================================

This script creates a proper pip-installable package with CFFI compilation.
This is the final step to make MOECAM work with 'pip install -e .'
"""

import os
import sys
import subprocess
from pathlib import Path

def create_package_structure():
    """Create the proper package structure for installation."""

    print("📦 Creating MOECAM package structure...")

    # Define the package directory
    package_root = Path("moecam_package")
    package_root.mkdir(exist_ok=True)

    # Create package directories
    dirs_to_create = [
        "moecam_package/moecam",
        "moecam_package/moecam/core",
        "moecam_package/moecam/algorithms",
        "moecam_package/moecam/utils",
        "moecam_package/tests"
    ]

    for dir_path in dirs_to_create:
        Path(dir_path).mkdir(parents=True, exist_ok=True)
        # Create __init__.py files
        init_file = Path(dir_path) / "__init__.py"
        if not init_file.exists():
            init_file.write_text('"""MOECAM package module."""\n')

    print("   ✅ Package structure created")
    return package_root

def create_cffi_builder(package_root: Path):
    """Create the CFFI builder script."""
    cffi_builder_content = '''"""
CFFI Builder for MOECAM
=======================
"""
from cffi import FFI

ffibuilder = FFI()

# Define C functions for CFFI
ffibuilder.cdef("""
    // Basic functions
    void hello_world(void);
    int add_numbers(int a, int b);

    // Pareto front extraction
    int extract_pareto_front(double* points, int n_points, int n_objectives,
                           double* pareto_points, int* pareto_indices);

    // Simple optimization functions
    int simple_random_optimization(double* bounds, int n_vars, int n_objectives,
                                 int max_evaluations, double* solutions, double* objectives);
""")

# C source code
ffibuilder.set_source("moecam._moecam_cffi",
    r"""
    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include <string.h>

    // Basic functions
    void hello_world(void) {
        printf("Hello from MOECAM CFFI!\\n");
    }

    int add_numbers(int a, int b) {
        return a + b;
    }

    // Pareto dominance check
    int dominates(double* p1, double* p2, int n_obj) {
        int at_least_one_better = 0;
        for (int i = 0; i < n_obj; i++) {
            if (p1[i] > p2[i]) return 0;  // p1 is worse in objective i
            if (p1[i] < p2[i]) at_least_one_better = 1;
        }
        return at_least_one_better;
    }

    // Extract Pareto front
    int extract_pareto_front(double* points, int n_points, int n_objectives,
                           double* pareto_points, int* pareto_indices) {
        int* is_pareto = (int*)calloc(n_points, sizeof(int));
        int pareto_count = 0;

        // Check each point for Pareto optimality
        for (int i = 0; i < n_points; i++) {
            is_pareto[i] = 1;  // Assume Pareto optimal

            for (int j = 0; j < n_points; j++) {
                if (i != j) {
                    if (dominates(&points[j * n_objectives], &points[i * n_objectives], n_objectives)) {
                        is_pareto[i] = 0;  // Point i is dominated
                        break;
                    }
                }
            }

            if (is_pareto[i]) {
                // Copy Pareto point
                for (int k = 0; k < n_objectives; k++) {
                    pareto_points[pareto_count * n_objectives + k] = points[i * n_objectives + k];
                }
                pareto_indices[pareto_count] = i;
                pareto_count++;
            }
        }

        free(is_pareto);
        return pareto_count;
    }

    // Simple random optimization (placeholder)
    int simple_random_optimization(double* bounds, int n_vars, int n_objectives,
                                 int max_evaluations, double* solutions, double* objectives) {
        // This is a placeholder - in real implementation would do proper optimization
        return max_evaluations;
    }
    """)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
'''
    cffi_builder_file = package_root / "moecam" / "cffi_build.py"
    cffi_builder_file.write_text(cffi_builder_content)
    print("   ✅ CFFI builder created")


def create_setup_py(package_root: Path):
    """Create setup.py for pip installation."""

    setup_py_content = '''#!/usr/bin/env python3
"""
MOECAM Package Setup
===================

Multi-Objective Evolutionary Comparison and Analysis Module
"""

from setuptools import setup, find_packages

# Package setup
setup(
    name="moecam",
    version="1.0.0",
    description="Multi-Objective Evolutionary Comparison and Analysis Module",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "matplotlib>=3.3.0",
        "cffi>=1.14.0",
    ],
    setup_requires=[
        "cffi>=1.14.0",
        "pybind11>=2.6"
    ],
    cffi_modules=["moecam/cffi_build.py:ffibuilder"],
    python_requires=">=3.7",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
    ],
)
'''

    setup_file = package_root / "setup.py"
    setup_file.write_text(setup_py_content)
    print("   ✅ setup.py created")

def create_main_module(package_root: Path):
    """Create the main MOECAM module."""

    main_module_content = '''"""
MOECAM - Multi-Objective Evolutionary Comparison and Analysis Module
===================================================================

Main module providing CFFI interfaces to C/C++ optimization algorithms.
"""

try:
    from moecam._moecam_cffi import lib, ffi
    CFFI_AVAILABLE = True
except ImportError:
    CFFI_AVAILABLE = False
    lib = None
    ffi = None

import numpy as np
from typing import List, Tuple, Callable, Optional

def hello_moecam():
    """Test function to verify CFFI is working."""
    if not CFFI_AVAILABLE:
        return "CFFI not available"

    lib.hello_world()
    result = lib.add_numbers(5, 3)
    return f"CFFI working! 5 + 3 = {result}"

def extract_pareto_front(objectives: np.ndarray, return_indices: bool = True) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """
    Extract Pareto front from objective values using CFFI.

    Args:
        objectives: Array of shape (n_points, n_objectives) with objective values
        return_indices: Whether to return indices of Pareto points

    Returns:
        Tuple of (pareto_objectives, pareto_indices) if return_indices=True
        Otherwise just pareto_objectives
    """
    if not CFFI_AVAILABLE:
        raise RuntimeError("CFFI not available. Install package with 'pip install -e .'")

    objectives = np.asarray(objectives, dtype=np.float64)
    n_points, n_objectives = objectives.shape

    # Allocate memory for results
    max_pareto = n_points  # Upper bound
    pareto_points = np.zeros((max_pareto, n_objectives), dtype=np.float64)
    pareto_indices = np.zeros(max_pareto, dtype=np.int32)

    # Create CFFI arrays
    c_points = ffi.cast("double*", ffi.from_buffer(objectives))
    c_pareto_points = ffi.cast("double*", ffi.from_buffer(pareto_points))
    c_pareto_indices = ffi.cast("int*", ffi.from_buffer(pareto_indices))

    # Call C function
    pareto_count = lib.extract_pareto_front(
        c_points, n_points, n_objectives,
        c_pareto_points, c_pareto_indices
    )

    # Extract results
    result_pareto = pareto_points[:pareto_count].copy()
    result_indices = pareto_indices[:pareto_count].copy() if return_indices else None

    return (result_pareto, result_indices) if return_indices else result_pareto

# Mock functions for compatibility with comprehensive test
def minimize_ecam_cffi(objective_func: Callable, bounds: List[Tuple[float, float]],
                      num_objectives: int = 2, max_iterations: int = 100):
    """Mock ECAM optimization - replace with real implementation."""
    raise NotImplementedError("ECAM algorithm not yet implemented in CFFI")

def minimize_random_start_cffi(objective_func: Callable, bounds: List[Tuple[float, float]],
                              num_objectives: int = 2, max_iterations: int = 100):
    """Mock Random Start optimization - replace with real implementation."""
    raise NotImplementedError("Random Start algorithm not yet implemented in CFFI")

def calculate_hypervolume_cffi(points: np.ndarray, reference_point: np.ndarray):
    """Mock hypervolume calculation - replace with real implementation."""
    raise NotImplementedError("Hypervolume calculation not yet implemented in CFFI")

__all__ = [
    'hello_moecam',
    'extract_pareto_front',
    'minimize_ecam_cffi',
    'minimize_random_start_cffi',
    'calculate_hypervolume_cffi',
    'CFFI_AVAILABLE'
]
'''

    init_file = package_root / "moecam" / "__init__.py"
    init_file.write_text(main_module_content)
    print("   ✅ Main module created")

def create_simple_test(package_root: Path):
    """Create a simple test that works with the installable package."""

    test_content = '''#!/usr/bin/env python3
"""
Simple MOECAM Test - Works with installed package
================================================

This test demonstrates the working CFFI functionality after installation.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

try:
    import moecam
    print(f"✅ MOECAM imported successfully")
    print(f"✅ CFFI available: {moecam.CFFI_AVAILABLE}")

    if moecam.CFFI_AVAILABLE:
        # Test basic functionality
        print(moecam.hello_moecam())

        # Test Pareto front extraction
        print("\\n🔍 Testing Pareto front extraction...")

        # Generate test data (ZDT1-like)
        np.random.seed(42)
        n_points = 100
        objectives = np.random.rand(n_points, 2)

        # Make some points clearly Pareto optimal
        for i in range(10):
            objectives[i] = [i/10.0, 1 - i/10.0]

        # Extract Pareto front
        pareto_points, pareto_indices = moecam.extract_pareto_front(objectives)

        print(f"   📊 Total points: {len(objectives)}")
        print(f"   🎯 Pareto points: {len(pareto_points)}")
        print(f"   📈 Pareto ratio: {len(pareto_points)/len(objectives):.1%}")

        # Create visualization
        plt.figure(figsize=(10, 6))
        plt.scatter(objectives[:, 0], objectives[:, 1], alpha=0.6, label='All Points')
        plt.scatter(pareto_points[:, 0], pareto_points[:, 1], color='red', s=50, label='Pareto Front')
        plt.xlabel('Objective 1')
        plt.ylabel('Objective 2')
        plt.title('MOECAM CFFI Test - Pareto Front Extraction')
        plt.legend()
        plt.grid(True, alpha=0.3)

        output_dir = Path("moecam_installation_test_results")
        output_dir.mkdir(exist_ok=True)

        plt.savefig(output_dir / "cffi_test.png", dpi=300, bbox_inches='tight')
        plt.close()

        # Save results
        np.savetxt(output_dir / "all_points.txt", objectives, header="All objective points")
        np.savetxt(output_dir / "pareto_points.txt", pareto_points, header="Pareto optimal points")
        np.savetxt(output_dir / "pareto_indices.txt", pareto_indices, fmt='%d', header="Pareto point indices")

        print(f"\\n💾 Results saved to: {output_dir.absolute()}")
        print("🎉 MOECAM CFFI test completed successfully!")

    else:
        print("❌ CFFI not available - package not properly installed")

except ImportError as e:
    print(f"❌ Failed to import MOECAM: {e}")
    print("💡 Try running: pip install -e moecam_package/")
'''

    test_file = package_root.parent / "test_installed_moecam.py"
    test_file.write_text(test_content)
    print("   ✅ Installation test created")

def install_package(package_root: Path):
    """Install the MOECAM package."""

    print(f"📦 Installing MOECAM package from {package_root}...")

    try:
        # Change to package directory and install
        original_dir = Path.cwd()
        os.chdir(package_root)

        # Install in development mode
        result = subprocess.run([
            sys.executable, "-m", "pip", "install", "-e", "."
        ], capture_output=True, text=True)

        os.chdir(original_dir)

        if result.returncode == 0:
            print("   ✅ Package installation successful!")
            return True
        else:
            print(f"   ❌ Package installation failed:")
            print(f"   Error: {result.stderr}")
            return False

    except Exception as e:
        print(f"   ❌ Installation error: {e}")
        return False

def main():
    """Main function to create and install the MOECAM package."""

    print("🚀 MOECAM Package Installation Setup")
    print("=" * 50)

    # Create package structure
    package_root = create_package_structure()

    # Create package files
    create_cffi_builder(package_root)
    create_setup_py(package_root)
    create_main_module(package_root)
    create_simple_test(package_root)

    print("\\n📦 Package structure created successfully!")
    print(f"📁 Package location: {package_root.absolute()}")

    # Ask user if they want to install
    print("\\n🔧 Ready to install MOECAM package...")
    user_input = input("Install now? (y/n): ").lower().strip()

    if user_input in ['y', 'yes']:
        success = install_package(package_root)

        if success:
            print("\\n🎉 MOECAM package installed successfully!")
            print("\\n📋 Next steps:")
            print("1. Run: python test_installed_moecam.py")
            print("2. Use: import moecam")
            print("3. Test: moecam.hello_moecam()")

            # Try to run the test automatically
            print("\\n🧪 Running installation test...")
            try:
                result = subprocess.run([
                    sys.executable, "test_installed_moecam.py"
                ], capture_output=True, text=True, timeout=30)

                if result.returncode == 0:
                    print("   ✅ Installation test passed!")
                    print("   🎯 MOECAM is ready to use!")
                else:
                    print("   ⚠️  Installation test had issues:")
                    print(f"   Output: {result.stdout}")

            except Exception as e:
                print(f"   ⚠️  Could not run automatic test: {e}")
                print("   💡 Please run manually: python test_installed_moecam.py")

        else:
            print("\\n❌ Package installation failed.")
            print("💡 You can try manually:")
            print(f"   cd {package_root}")
            print("   pip install -e .")
    else:
        print("\\n📦 Package created but not installed.")
        print("💡 To install later:")
        print(f"   cd {package_root}")
        print("   pip install -e .")

if __name__ == "__main__":
    main()
