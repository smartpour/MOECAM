#!/usr/bin/env python3
"""
MOECAM CFFI Build Script
=======================

Builds all C++ extensions using CFFI during package installation.
This follows the standard CFFI installation process and is platform-independent.
"""

import cffi
import os
from pathlib import Path

def build_moecam_extension():
    """Build the main MOECAM extension with all C++ sources included."""

    ffibuilder = cffi.FFI()

    # Define the complete C interface
    ffibuilder.cdef("""
        // Callback type for objective functions
        typedef void (*objective_callback_t)(const double* x, int n_dim, double* f, int n_obj);

        // Toy example functions (legacy)
        double add(double a, double b);
        void print_message(const char* message);

        // Pareto Front Extraction
        int extract_pareto_front(double* all_points, int num_points, int num_objectives,
                                double* pareto_points, int max_pareto_points, int strict_mode);
        int extract_pareto_indices(double* all_points, int num_points, int num_objectives,
                                  int* pareto_indices, int max_pareto_points, int strict_mode);

        // WFG Hypervolume
        int wfg_init(int num_objectives, int max_points);
        double wfg_calculate_hypervolume(double* points, double* reference_point, int num_points, int num_objectives);
        void wfg_cleanup();
    """)

    # Find source directory
    this_dir = Path(__file__).parent
    src_dir = this_dir / "src" / "moecam" / "core"
    csources_dir = this_dir.parent / "csources"
    moecam_sources_dir = csources_dir / "moecam"
    wfg_dir = csources_dir / "wfg" / "WFG_1.15"

    # Collect all C++ source files that will be compiled together
    sources = []

    # Add toy example (existing working code)
    toy_cpp = src_dir / "toy_cpp_lib.cpp"
    if toy_cpp.exists():
        sources.append(str(toy_cpp))

    # Add our wrapper files (skip MOECAM for now due to complex dependencies)
    pareto_wrapper = src_dir / "pareto_wrapper.cpp"
    wfg_wrapper = src_dir / "wfg_wrapper.cpp"

    if pareto_wrapper.exists():
        sources.append(str(pareto_wrapper))
    if wfg_wrapper.exists():
        sources.append(str(wfg_wrapper))

    # Add WFG sources
    wfg_c = wfg_dir / "wfg.c"
    read_c = wfg_dir / "read.c"
    if wfg_c.exists():
        sources.append(str(wfg_c))
    if read_c.exists():
        sources.append(str(read_c))

    # TODO: Add MOECAM sources when dependencies are resolved
    # gansoc_cpp = moecam_sources_dir / "gansoc.cpp"
    # if gansoc_cpp.exists():
    #     sources.append(str(gansoc_cpp))

    print(f"Building MOECAM extension with sources: {sources}")

    # Set compiler arguments and include directories
    include_dirs = []
    if wfg_dir.exists():
        include_dirs.append(str(wfg_dir))
    if moecam_sources_dir.exists():
        include_dirs.append(str(moecam_sources_dir))

    # Prepare source content for inline compilation
    source_content = """
        // Include all necessary headers and implementations

        // Toy example implementation
        extern "C" {
            double add(double a, double b) { return a + b; }
            void print_message(const char* message) {
                std::cout << "C++ says: " << message << std::endl;
            }
        }

        // Include wrapper implementations
        // Note: The actual implementations will be compiled from source files
    """

    # Build the extension - separate C and C++ compilation
    ffibuilder.set_source(
        "moecam._moecam_cffi",
        source_content,
        sources=sources,
        include_dirs=include_dirs,
        extra_compile_args=["-O3", "-DNDEBUG"],
        extra_link_args=["-lstdc++", "-lm"],
        libraries=["m"]  # Math library
    )

    return ffibuilder

if __name__ == "__main__":
    ffibuilder = build_moecam_extension()
    ffibuilder.compile(verbose=True)
