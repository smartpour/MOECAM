/**
 * MOECAM Algorithms CFFI Wrapper
 * ==============================
 *
 * C++ wrapper for MOECAM optimization algorithms that can be compiled
 * as a shared library using CFFI during package installation.
 *
 * Based on example3.cpp but converted to library functions without main().
 */

#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <set>
#include <cstdlib>
#include <string>

using namespace std;

// Include MOECAM headers - need to adapt paths during compilation
#include "gansoc.h"

extern "C" {

// Callback type for objective function evaluation from Python
typedef void (*objective_callback_t)(const double* x, int n_dim, double* f, int n_obj);

// Global state for MOECAM algorithms
static objective_callback_t g_python_callback = nullptr;
static int g_current_objective = 0;
static int g_total_objectives = 1;
static int g_evaluation_count = 0;
static vector<vector<double>> g_all_solutions;
static vector<vector<double>> g_all_objectives;

// Wrapper function to bridge MOECAM's function signature to our callback
void moecam_objective_wrapper(int* n, double* x, double* f, int* objective_index) {
    if (!g_python_callback) {
        *f = 1e6; // Large value on error
        return;
    }

    g_evaluation_count++;

    // Store solution for later retrieval
    vector<double> solution(*n);
    for (int i = 0; i < *n; i++) {
        solution[i] = x[i];
    }

    // If this is the first objective, allocate space for all objectives
    if (*objective_index == 0) {
        g_all_solutions.push_back(solution);
        g_all_objectives.push_back(vector<double>(g_total_objectives));
    }

    // Calculate all objectives at once on first call
    if (*objective_index == 0) {
        double* all_f = new double[g_total_objectives];
        g_python_callback(x, *n, all_f, g_total_objectives);

        // Store all objective values
        for (int i = 0; i < g_total_objectives; i++) {
            g_all_objectives.back()[i] = all_f[i];
        }

        *f = all_f[0]; // Return first objective
        delete[] all_f;
    } else if (*objective_index < g_total_objectives) {
        // Return cached objective value
        *f = g_all_objectives.back()[*objective_index];
    } else {
        *f = 1e6; // Invalid objective index
    }
}

/**
 * Initialize MOECAM algorithm configuration
 * @param n_dim Number of decision variables
 * @param n_obj Number of objectives
 * @param lower_bounds Array of lower bounds
 * @param upper_bounds Array of upper bounds
 * @param callback Python objective function callback
 * @return Configuration ID (0) on success, -1 on error
 */
int moecam_init(int n_dim, int n_obj, double* lower_bounds, double* upper_bounds,
                objective_callback_t callback) {
    if (n_dim <= 0 || n_obj <= 0 || !lower_bounds || !upper_bounds || !callback) {
        return -1;
    }

    // Store global state
    g_python_callback = callback;
    g_total_objectives = n_obj;
    g_evaluation_count = 0;
    g_all_solutions.clear();
    g_all_objectives.clear();

    return 0; // Simple single-configuration system
}

/**
 * Run ECAM algorithm
 * @param config_id Configuration ID (ignored in simple implementation)
 * @param n_dim Number of decision variables
 * @param lower_bounds Array of lower bounds
 * @param upper_bounds Array of upper bounds
 * @param max_iterations Maximum number of iterations
 * @param best_solution Output: best solution found
 * @param best_objectives Output: best objective values
 * @param all_solutions_out Output: all evaluated solutions (flattened)
 * @param all_objectives_out Output: all objective evaluations (flattened)
 * @param max_solutions Maximum number of solutions to return
 * @param num_evaluations Output: actual number of function evaluations
 * @return 0 on success, -1 on error
 */
int moecam_run_ecam(int config_id, int n_dim, double* lower_bounds, double* upper_bounds,
                   int max_iterations, double* best_solution, double* best_objectives,
                   double* all_solutions_out, double* all_objectives_out,
                   int max_solutions, int* num_evaluations) {

    if (!g_python_callback || n_dim <= 0) {
        return -1;
    }

    // Initialize Ganso optimizer
    Ganso optimizer;
    optimizer.m_objectives = g_total_objectives;

    // Set up initial solution
    double* x0 = new double[n_dim];
    for (int i = 0; i < n_dim; i++) {
        x0[i] = (lower_bounds[i] + upper_bounds[i]) / 2.0; // Start at middle
    }

    double val;
    int retcode;

    // Run ECAM algorithm
    retcode = optimizer.MinimizeECAM(n_dim, x0, &val, moecam_objective_wrapper,
                                   0, 0, 70, NULL, NULL, NULL, NULL,
                                   lower_bounds, upper_bounds, NULL, -max_iterations);

    // Copy best solution
    for (int i = 0; i < n_dim; i++) {
        best_solution[i] = x0[i];
    }

    // For single-objective case, copy the best objective
    if (g_total_objectives == 1) {
        best_objectives[0] = val;
    } else {
        // For multi-objective, find the last evaluation's objectives
        if (!g_all_objectives.empty()) {
            for (int i = 0; i < g_total_objectives; i++) {
                best_objectives[i] = g_all_objectives.back()[i];
            }
        }
    }

    // Copy all solutions and objectives to output arrays
    int solutions_to_copy = min((int)g_all_solutions.size(), max_solutions);
    for (int i = 0; i < solutions_to_copy; i++) {
        // Copy solution
        for (int j = 0; j < n_dim; j++) {
            all_solutions_out[i * n_dim + j] = g_all_solutions[i][j];
        }

        // Copy objectives
        for (int j = 0; j < g_total_objectives; j++) {
            all_objectives_out[i * g_total_objectives + j] = g_all_objectives[i][j];
        }
    }

    *num_evaluations = g_evaluation_count;

    delete[] x0;
    return retcode == 0 ? 0 : -1;
}

/**
 * Run Random Start algorithm
 * @param config_id Configuration ID (ignored in simple implementation)
 * @param n_dim Number of decision variables
 * @param lower_bounds Array of lower bounds
 * @param upper_bounds Array of upper bounds
 * @param max_iterations Maximum number of iterations
 * @param best_solution Output: best solution found
 * @param best_objectives Output: best objective values
 * @param all_solutions_out Output: all evaluated solutions (flattened)
 * @param all_objectives_out Output: all objective evaluations (flattened)
 * @param max_solutions Maximum number of solutions to return
 * @param num_evaluations Output: actual number of function evaluations
 * @return 0 on success, -1 on error
 */
int moecam_run_random_start(int config_id, int n_dim, double* lower_bounds, double* upper_bounds,
                           int max_iterations, double* best_solution, double* best_objectives,
                           double* all_solutions_out, double* all_objectives_out,
                           int max_solutions, int* num_evaluations) {

    if (!g_python_callback || n_dim <= 0) {
        return -1;
    }

    // Initialize Ganso optimizer
    Ganso optimizer;
    optimizer.m_objectives = g_total_objectives;

    // Set up initial solution
    double* x0 = new double[n_dim];
    for (int i = 0; i < n_dim; i++) {
        x0[i] = (lower_bounds[i] + upper_bounds[i]) / 2.0; // Start at middle
    }

    double val;
    int retcode;

    // Run Random Start algorithm
    retcode = optimizer.MinimizeRandomStart(n_dim, x0, &val, moecam_objective_wrapper,
                                          0, 0, NULL, NULL, NULL, NULL,
                                          lower_bounds, upper_bounds, NULL, max_iterations);

    // Copy best solution
    for (int i = 0; i < n_dim; i++) {
        best_solution[i] = x0[i];
    }

    // For single-objective case, copy the best objective
    if (g_total_objectives == 1) {
        best_objectives[0] = val;
    } else {
        // For multi-objective, find the last evaluation's objectives
        if (!g_all_objectives.empty()) {
            for (int i = 0; i < g_total_objectives; i++) {
                best_objectives[i] = g_all_objectives.back()[i];
            }
        }
    }

    // Copy all solutions and objectives to output arrays
    int solutions_to_copy = min((int)g_all_solutions.size(), max_solutions);
    for (int i = 0; i < solutions_to_copy; i++) {
        // Copy solution
        for (int j = 0; j < n_dim; j++) {
            all_solutions_out[i * n_dim + j] = g_all_solutions[i][j];
        }

        // Copy objectives
        for (int j = 0; j < g_total_objectives; j++) {
            all_objectives_out[i * g_total_objectives + j] = g_all_objectives[i][j];
        }
    }

    *num_evaluations = g_evaluation_count;

    delete[] x0;
    return retcode == 0 ? 0 : -1;
}

/**
 * Clean up MOECAM resources
 * @param config_id Configuration ID (ignored in simple implementation)
 */
void moecam_cleanup(int config_id) {
    g_python_callback = nullptr;
    g_evaluation_count = 0;
    g_all_solutions.clear();
    g_all_objectives.clear();
}

} // extern "C"
