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

        // Optimizer interface
        int create_optimizer(int n_dim, int n_obj, double* lower_bounds, double* upper_bounds, objective_callback_t callback);
        void destroy_optimizer(int handle);
        void generate_random_solution(int handle, double* solution);
        void crossover_solutions(int handle, const double* parent1, const double* parent2, double* child1, double* child2);
        void mutate_solution(int handle, double* individual, double mutation_rate);
        void evaluate_solution(int handle, const double* solution, double* objectives);

        // MOECAM Algorithms
        typedef struct {
            int algorithm_id;
            int n_dim;
            int n_obj;
            int max_iterations;
            double* lower_bounds;
            double* upper_bounds;
            objective_callback_t objective_func;
        } moecam_config_t;

        int moecam_init(moecam_config_t* config);
        int moecam_run_ecam(int config_id, double* best_solution, double* best_objectives,
                           double* all_solutions, double* all_objectives, int* num_evaluations);
        int moecam_run_random_start(int config_id, double* best_solution, double* best_objectives,
                                   double* all_solutions, double* all_objectives, int* num_evaluations);
        int moecam_run_direct(int config_id, double* best_solution, double* best_objective, int* num_evaluations);
        void moecam_cleanup(int config_id);

        // Pareto Front Extraction
        int extract_pareto_front(double* all_points, int num_points, int num_objectives,
                                double* pareto_points, int max_pareto_points);

        // WFG Hypervolume
        int wfg_init(int num_objectives, int max_points);
        double wfg_calculate_hypervolume(double* points, double* reference_point, int num_points, int num_objectives);
        void wfg_cleanup();
    """)

    # Find source directory
    this_dir = Path(__file__).parent
    src_dir = this_dir / "src" / "src" / "moecam" / "core"
    scripts_dir = this_dir / "src" / "tests" / "scripts"
    wfg_dir = this_dir / "csources" / "wfg" / "WFG_1.15"

    # Collect all C++ source files
    sources = []

    # Add toy example (existing working code)
    toy_cpp = src_dir / "toy_cpp_lib.cpp"
    if toy_cpp.exists():
        sources.append(str(toy_cpp))

    # Add MOECAM algorithms wrapper
    moecam_cpp = scripts_dir / "moecam_wrapper.cpp"
    if moecam_cpp.exists():
        sources.append(str(moecam_cpp))

    # Add Pareto extraction wrapper
    pareto_cpp = scripts_dir / "pareto_wrapper.cpp"
    if pareto_cpp.exists():
        sources.append(str(pareto_cpp))

    # Add WFG sources
    wfg_c = wfg_dir / "wfg.c"
    read_c = wfg_dir / "read.c"
    if wfg_c.exists() and read_c.exists():
        sources.extend([str(wfg_c), str(read_c)])

    # Add WFG wrapper
    wfg_cpp = scripts_dir / "wfg_wrapper.cpp"
    if wfg_cpp.exists():
        sources.append(str(wfg_cpp))

    print(f"Building MOECAM extension with sources: {sources}")

    # Set compiler arguments
    include_dirs = [str(wfg_dir)] if wfg_dir.exists() else []

    # Build the extension
    ffibuilder.set_source(
        "moecam._moecam_cffi",
        """
        // Include all necessary headers and implementations
        #include "toy_cpp_lib.cpp"
        #include "moecam_wrapper.cpp"
        #include "pareto_wrapper.cpp"
        #include "wfg_wrapper.cpp"
        """,
        sources=sources,
        include_dirs=include_dirs,
        extra_compile_args=["-std=c++17", "-O3"],
        extra_link_args=["-lstdc++"]
    )

    return ffibuilder

if __name__ == "__main__":
    ffibuilder = build_moecam_extension()
    ffibuilder.compile(verbose=True)
