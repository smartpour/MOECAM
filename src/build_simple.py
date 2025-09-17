#!/usr/bin/env python3
"""
Simple MOECAM CFFI Build Script
==============================

Minimal working CFFI build for testing basic functionality.
"""

import cffi
import os
from pathlib import Path

def build_simple_extension():
    """Build a simple MOECAM extension with basic functionality."""

    ffibuilder = cffi.FFI()

    # Define the basic C interface
    ffibuilder.cdef("""
        // Toy example functions
        double add(double a, double b);
        void print_message(const char* message);

        // Pareto Front Extraction
        int extract_pareto_front(double* all_points, int num_points, int num_objectives,
                                double* pareto_points, int max_pareto_points, int strict_mode);
        int extract_pareto_indices(double* all_points, int num_points, int num_objectives,
                                  int* pareto_indices, int max_pareto_points, int strict_mode);
    """)

    # Find source directory
    this_dir = Path(__file__).parent
    src_dir = this_dir / "src" / "moecam" / "core"

    # Collect C++ source files
    sources = []

    # Add toy example (existing working code)
    toy_cpp = src_dir / "toy_cpp_lib.cpp"
    if toy_cpp.exists():
        sources.append(str(toy_cpp))

    # Add pareto wrapper
    pareto_wrapper = src_dir / "pareto_wrapper.cpp"
    if pareto_wrapper.exists():
        sources.append(str(pareto_wrapper))

    print(f"Building simple MOECAM extension with sources: {sources}")

    # Build the extension
    ffibuilder.set_source(
        "moecam._moecam_cffi",
        """
        // Simple inline implementation
        #include <iostream>

        extern "C" {
            double add(double a, double b) { return a + b; }
            void print_message(const char* message) {
                std::cout << "C++ says: " << message << std::endl;
            }
        }
        """,
        sources=sources,
        extra_compile_args=["-O3", "-std=c++17"],
        extra_link_args=["-lstdc++"]
    )

    return ffibuilder

if __name__ == "__main__":
    ffibuilder = build_simple_extension()
    ffibuilder.compile(verbose=True)
