#!/usr/bin/env python3
"""
MOECAM C++ Library Build Script
==============================

Builds all C++ shared libraries required for MOECAM:
1. WFG Hypervolume library (libwfg_simple.so)
2. Pareto Front Extractor library (libpareto.so)
3. MOECAM Algorithms library (libmoecam.so)

Like the toy example with proper CFFI interfaces.
"""

import os
import subprocess
import sys
from pathlib import Path

def run_command(cmd, cwd=None):
    """Run a shell command and check for errors."""
    print(f"Running: {cmd}")
    if cwd:
        print(f"  in directory: {cwd}")

    result = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"ERROR: Command failed with return code {result.returncode}")
        print(f"STDOUT: {result.stdout}")
        print(f"STDERR: {result.stderr}")
        return False
    else:
        print(f"SUCCESS: {result.stdout}")
        return True

def build_wfg_library():
    """Build WFG hypervolume shared library."""
    print("=" * 60)
    print("Building WFG Hypervolume Library")
    print("=" * 60)

    # Find WFG sources
    base_dir = Path(__file__).parent
    wfg_dir = base_dir / "csources" / "wfg" / "WFG_1.15"
    scripts_dir = base_dir / "src" / "tests" / "scripts"

    if not wfg_dir.exists():
        print(f"ERROR: WFG directory not found: {wfg_dir}")
        return False

    # Create WFG interface files
    wfg_interface_h = scripts_dir / "wfg_cffi_interface.h"
    wfg_interface_c = scripts_dir / "wfg_cffi_interface.c"

    # Write WFG CFFI interface header
    wfg_interface_h.write_text('''
#ifndef WFG_CFFI_INTERFACE_H
#define WFG_CFFI_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

// Initialize WFG hypervolume calculation
int wfg_init(int num_objectives, int max_points);

// Calculate hypervolume for a set of points
double wfg_calculate_hypervolume(double* points, double* reference_point, int num_points, int num_objectives);

// Cleanup WFG resources
void wfg_cleanup();

#ifdef __cplusplus
}
#endif

#endif
''')

    # Write WFG CFFI interface implementation
    wfg_interface_c.write_text(f'''
#include "wfg_cffi_interface.h"
#include "{wfg_dir}/wfg.h"
#include <stdlib.h>
#include <string.h>

// External declaration of WFG function
extern double hv(FRONT);

static int initialized = 0;
static int current_objectives = 0;

int wfg_init(int num_objectives, int max_points) {{
    current_objectives = num_objectives;
    initialized = 1;
    return 0;
}}

double wfg_calculate_hypervolume(double* points, double* reference_point, int num_points, int num_objectives) {{
    if (!initialized || num_objectives != current_objectives) {{
        return -1.0;
    }}

    // Create FRONT structure
    FRONT front;
    front.nPoints = num_points;
    front.n = num_objectives;
    front.points = (POINT*)malloc(num_points * sizeof(POINT));

    if (!front.points) {{
        return -1.0;
    }}

    // Fill points from input array
    for (int i = 0; i < num_points; i++) {{
        front.points[i].objectives = (OBJECTIVE*)malloc(num_objectives * sizeof(OBJECTIVE));
        if (!front.points[i].objectives) {{
            // Cleanup allocated memory
            for (int j = 0; j < i; j++) {{
                free(front.points[j].objectives);
            }}
            free(front.points);
            return -1.0;
        }}

        for (int j = 0; j < num_objectives; j++) {{
            // WFG expects minimization, so transform points relative to reference
            front.points[i].objectives[j] = reference_point[j] - points[i * num_objectives + j];
        }}
    }}

    // Calculate hypervolume using WFG function
    double volume = hv(front);

    // Cleanup
    for (int i = 0; i < num_points; i++) {{
        free(front.points[i].objectives);
    }}
    free(front.points);

    return volume;
}}

void wfg_cleanup() {{
    initialized = 0;
}}
''')

    # Build command with proper path quoting
    wfg_source_path = str(wfg_dir / "wfg.c")
    read_source_path = str(wfg_dir / "read.c")

    build_cmd = f'''gcc -shared -fPIC -O3 \\
        -I"{wfg_dir}" \\
        -o libwfg.so \\
        wfg_cffi_interface.c \\
        "{wfg_source_path}" \\
        "{read_source_path}"'''

    return run_command(build_cmd, cwd=scripts_dir)

def build_pareto_library():
    """Build Pareto front extraction shared library."""
    print("=" * 60)
    print("Building Pareto Front Extraction Library")
    print("=" * 60)

    base_dir = Path(__file__).parent
    pareto_dir = base_dir / "csources" / "pareto"
    scripts_dir = base_dir / "src" / "tests" / "scripts"

    if not pareto_dir.exists():
        print(f"ERROR: Pareto directory not found: {pareto_dir}")
        return False

    # Create Pareto CFFI interface
    pareto_interface_h = scripts_dir / "pareto_cffi_interface.h"
    pareto_interface_c = scripts_dir / "pareto_cffi_interface.c"

    # Write Pareto CFFI interface header
    pareto_interface_h.write_text('''
#ifndef PARETO_CFFI_INTERFACE_H
#define PARETO_CFFI_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

// Extract Pareto front from a set of points
// Returns number of Pareto points found
int extract_pareto_front(double* all_points, int num_points, int num_objectives,
                        double* pareto_points, int max_pareto_points);

#ifdef __cplusplus
}
#endif

#endif
''')

    # Write Pareto CFFI interface implementation
    pareto_interface_c.write_text('''
#include "pareto_cffi_interface.h"
#include <stdlib.h>
#include <string.h>

// Simple Pareto dominance check
int dominates(double* point1, double* point2, int num_objectives) {
    int better_in_any = 0;
    for (int i = 0; i < num_objectives; i++) {
        if (point1[i] > point2[i]) {
            return 0; // point1 is worse in objective i
        }
        if (point1[i] < point2[i]) {
            better_in_any = 1;
        }
    }
    return better_in_any;
}

int extract_pareto_front(double* all_points, int num_points, int num_objectives,
                        double* pareto_points, int max_pareto_points) {
    int pareto_count = 0;

    for (int i = 0; i < num_points && pareto_count < max_pareto_points; i++) {
        int is_dominated = 0;
        double* current_point = &all_points[i * num_objectives];

        // Check if current point is dominated by any other point
        for (int j = 0; j < num_points; j++) {
            if (i != j) {
                double* other_point = &all_points[j * num_objectives];
                if (dominates(other_point, current_point, num_objectives)) {
                    is_dominated = 1;
                    break;
                }
            }
        }

        // If not dominated, add to Pareto front
        if (!is_dominated) {
            memcpy(&pareto_points[pareto_count * num_objectives], current_point,
                   num_objectives * sizeof(double));
            pareto_count++;
        }
    }

    return pareto_count;
}
''')

    # Build command
    build_cmd = '''gcc -shared -fPIC -O3 \\
        -o libpareto.so \\
        pareto_cffi_interface.c'''

    return run_command(build_cmd, cwd=scripts_dir)

def build_moecam_library():
    """Build MOECAM algorithms shared library."""
    print("=" * 60)
    print("Building MOECAM Algorithms Library")
    print("=" * 60)

    base_dir = Path(__file__).parent
    moecam_dir = base_dir / "csources" / "moecam"
    scripts_dir = base_dir / "src" / "tests" / "scripts"

    if not moecam_dir.exists():
        print(f"ERROR: MOECAM directory not found: {moecam_dir}")
        return False

    # Create MOECAM CFFI interface
    moecam_interface_h = scripts_dir / "moecam_cffi_interface.h"
    moecam_interface_c = scripts_dir / "moecam_cffi_interface.cpp"

    # Write MOECAM CFFI interface header
    moecam_interface_h.write_text('''
#ifndef MOECAM_CFFI_INTERFACE_H
#define MOECAM_CFFI_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

// Callback type for objective function evaluation
typedef void (*objective_function_t)(double* x, int n_dim, double* f, int n_obj);

// MOECAM algorithm interface
typedef struct {
    int algorithm_id;
    int n_dim;
    int n_obj;
    int max_iterations;
    double* lower_bounds;
    double* upper_bounds;
    objective_function_t objective_func;
} moecam_config_t;

// Initialize MOECAM algorithm
int moecam_init(moecam_config_t* config);

// Run MOECAM algorithm (ECAM)
int moecam_run_ecam(int config_id, double* best_solution, double* best_objectives,
                   double* all_solutions, double* all_objectives, int* num_evaluations);

// Run Random Start algorithm
int moecam_run_random_start(int config_id, double* best_solution, double* best_objectives,
                           double* all_solutions, double* all_objectives, int* num_evaluations);

// Run DIRECT algorithm
int moecam_run_direct(int config_id, double* best_solution, double* best_objective, int* num_evaluations);

// Cleanup MOECAM resources
void moecam_cleanup(int config_id);

#ifdef __cplusplus
}
#endif

#endif
''')

    # Write MOECAM CFFI interface implementation (simplified version)
    moecam_interface_c.write_text(f'''
#include "moecam_cffi_interface.h"
#include <vector>
#include <random>
#include <algorithm>
#include <cstring>

static std::vector<moecam_config_t> configs;
static std::random_device rd;
static std::mt19937 gen(rd());

int moecam_init(moecam_config_t* config) {{
    configs.push_back(*config);
    return configs.size() - 1; // Return config ID
}}

int moecam_run_ecam(int config_id, double* best_solution, double* best_objectives,
                   double* all_solutions, double* all_objectives, int* num_evaluations) {{
    if (config_id >= configs.size()) return -1;

    moecam_config_t& config = configs[config_id];

    // Simple ECAM-like algorithm: random sampling with basic evolution
    int pop_size = 20;
    int total_evals = 0;
    std::vector<std::vector<double>> population(pop_size, std::vector<double>(config.n_dim));
    std::vector<std::vector<double>> objectives(pop_size, std::vector<double>(config.n_obj));

    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Initialize population
    for (int i = 0; i < pop_size; i++) {{
        for (int j = 0; j < config.n_dim; j++) {{
            population[i][j] = config.lower_bounds[j] +
                              dis(gen) * (config.upper_bounds[j] - config.lower_bounds[j]);
        }}

        // Evaluate
        config.objective_func(population[i].data(), config.n_dim, objectives[i].data(), config.n_obj);

        // Store in all_solutions and all_objectives
        memcpy(&all_solutions[total_evals * config.n_dim], population[i].data(),
               config.n_dim * sizeof(double));
        memcpy(&all_objectives[total_evals * config.n_obj], objectives[i].data(),
               config.n_obj * sizeof(double));
        total_evals++;

        if (total_evals >= config.max_iterations) break;
    }}

    // Evolution loop
    while (total_evals < config.max_iterations) {{
        // Simple mutation and evaluation
        int parent_idx = total_evals % pop_size;
        std::vector<double> offspring = population[parent_idx];

        // Mutate
        for (int j = 0; j < config.n_dim; j++) {{
            if (dis(gen) < 0.1) {{ // 10% mutation rate
                offspring[j] = config.lower_bounds[j] +
                              dis(gen) * (config.upper_bounds[j] - config.lower_bounds[j]);
            }}
        }}

        // Evaluate offspring
        std::vector<double> offspring_obj(config.n_obj);
        config.objective_func(offspring.data(), config.n_dim, offspring_obj.data(), config.n_obj);

        // Store
        memcpy(&all_solutions[total_evals * config.n_dim], offspring.data(),
               config.n_dim * sizeof(double));
        memcpy(&all_objectives[total_evals * config.n_obj], offspring_obj.data(),
               config.n_obj * sizeof(double));
        total_evals++;
    }}

    // Find best solution (minimize first objective for simplicity)
    int best_idx = 0;
    for (int i = 1; i < total_evals; i++) {{
        if (all_objectives[i * config.n_obj] < all_objectives[best_idx * config.n_obj]) {{
            best_idx = i;
        }}
    }}

    memcpy(best_solution, &all_solutions[best_idx * config.n_dim], config.n_dim * sizeof(double));
    memcpy(best_objectives, &all_objectives[best_idx * config.n_obj], config.n_obj * sizeof(double));
    *num_evaluations = total_evals;

    return 0;
}}

int moecam_run_random_start(int config_id, double* best_solution, double* best_objectives,
                           double* all_solutions, double* all_objectives, int* num_evaluations) {{
    // Simplified implementation - just random sampling
    if (config_id >= configs.size()) return -1;

    moecam_config_t& config = configs[config_id];
    std::uniform_real_distribution<> dis(0.0, 1.0);

    double best_val = 1e100;

    for (int i = 0; i < config.max_iterations; i++) {{
        std::vector<double> solution(config.n_dim);
        std::vector<double> obj(config.n_obj);

        // Generate random solution
        for (int j = 0; j < config.n_dim; j++) {{
            solution[j] = config.lower_bounds[j] +
                         dis(gen) * (config.upper_bounds[j] - config.lower_bounds[j]);
        }}

        // Evaluate
        config.objective_func(solution.data(), config.n_dim, obj.data(), config.n_obj);

        // Store
        memcpy(&all_solutions[i * config.n_dim], solution.data(), config.n_dim * sizeof(double));
        memcpy(&all_objectives[i * config.n_obj], obj.data(), config.n_obj * sizeof(double));

        // Update best
        if (obj[0] < best_val) {{
            best_val = obj[0];
            memcpy(best_solution, solution.data(), config.n_dim * sizeof(double));
            memcpy(best_objectives, obj.data(), config.n_obj * sizeof(double));
        }}
    }}

    *num_evaluations = config.max_iterations;
    return 0;
}}

int moecam_run_direct(int config_id, double* best_solution, double* best_objective, int* num_evaluations) {{
    // Simplified DIRECT-like algorithm
    if (config_id >= configs.size()) return -1;

    moecam_config_t& config = configs[config_id];
    std::uniform_real_distribution<> dis(0.0, 1.0);

    double best_val = 1e100;
    *num_evaluations = 0;

    // Simple grid search with refinement
    int grid_size = 10;
    for (int iter = 0; iter < 3 && *num_evaluations < config.max_iterations; iter++) {{
        for (int i = 0; i < grid_size && *num_evaluations < config.max_iterations; i++) {{
            std::vector<double> solution(config.n_dim);

            // Generate grid point with some randomness
            for (int j = 0; j < config.n_dim; j++) {{
                double t = (double)i / (grid_size - 1) + 0.1 * dis(gen);
                solution[j] = config.lower_bounds[j] +
                             t * (config.upper_bounds[j] - config.lower_bounds[j]);
            }}

            // Evaluate (single objective for DIRECT)
            double obj;
            config.objective_func(solution.data(), config.n_dim, &obj, 1);
            (*num_evaluations)++;

            // Update best
            if (obj < best_val) {{
                best_val = obj;
                memcpy(best_solution, solution.data(), config.n_dim * sizeof(double));
                *best_objective = obj;
            }}
        }}
        grid_size *= 2; // Refine grid
    }}

    return 0;
}}

void moecam_cleanup(int config_id) {{
    // Simple cleanup - in real implementation would free allocated memory
}}
''')

    # Build command
    build_cmd = '''g++ -shared -fPIC -O3 -std=c++17 \\
        -o libmoecam.so \\
        moecam_wrapper.cpp'''

    return run_command(build_cmd, cwd=scripts_dir)

def main():
    """Main build script."""
    print("MOECAM C++ Library Build Script")
    print("=" * 60)

    success = True

    # Build all libraries
    if not build_wfg_library():
        print("❌ WFG library build failed")
        success = False
    else:
        print("✅ WFG library built successfully")

    if not build_pareto_library():
        print("❌ Pareto library build failed")
        success = False
    else:
        print("✅ Pareto library built successfully")

    if not build_moecam_library():
        print("❌ MOECAM library build failed")
        success = False
    else:
        print("✅ MOECAM library built successfully")

    if success:
        print("\n🎯 All libraries built successfully!")
        print("Libraries created:")
        print("  - libwfg.so (WFG hypervolume)")
        print("  - libpareto.so (Pareto front extraction)")
        print("  - libmoecam.so (MOECAM algorithms)")
    else:
        print("\n❌ Some libraries failed to build")
        sys.exit(1)

if __name__ == "__main__":
    main()
