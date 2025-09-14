#include "moecam_python_interface.h"
#include "../../../csources/moecam/gansoc.h"
#include "../../../csources/moecam/DIRECT-MASTER/src/direct.h"
#include <cstring>
#include <vector>

// Global state for callback handling
static python_callback_func g_python_callback = nullptr;
static int g_num_objectives = 2;
static std::vector<double> g_objective_values;

// Convert Python callback to MOECAM format
void moecam_callback_wrapper(int *n, double *x, double *f, int *t) {
    if (!g_python_callback) {
        *f = 1e20; // Error value
        return;
    }

    // If this is the first call for this point, compute all objectives
    if (*t == 0) {
        g_objective_values.resize(g_num_objectives);
        g_python_callback(*n, x, g_objective_values.data(), g_num_objectives);
    }

    // Return the requested objective
    if (*t < g_num_objectives) {
        *f = g_objective_values[*t];
    } else {
        // This is called for Pareto front handling - return dummy value
        *f = 0.0;
    }
}

// Convert Python callback to DIRECT format (single objective)
double direct_callback_wrapper(int n, double *x, int *undefined_flag, void *data) {
    if (!g_python_callback) {
        *undefined_flag = 1;
        return 1e20;
    }

    std::vector<double> objectives(g_num_objectives);
    g_python_callback(n, x, objectives.data(), g_num_objectives);

    *undefined_flag = 0;
    return objectives[0]; // Use first objective for single-objective DIRECT
}

int moecam_init() {
    g_objective_values.clear();
    return 0;
}

void moecam_cleanup() {
    g_python_callback = nullptr;
    g_objective_values.clear();
}

int moecam_minimize_ecam(
    int dim,
    double *x0,
    double *result_value,
    python_callback_func callback,
    int num_objectives,
    const double *bounds_lower,
    const double *bounds_upper,
    int max_iterations
) {
    // Set up global state
    g_python_callback = callback;
    g_num_objectives = num_objectives;

    try {
        Ganso optimizer;

        // Copy bounds
        std::vector<double> lower_bounds(bounds_lower, bounds_lower + dim);
        std::vector<double> upper_bounds(bounds_upper, bounds_upper + dim);

        // Call MinimizeECAM
        // Parameters: dim, x0, val, function, lineq, linineq, LC, AE, AI, RHSE, RHSI, Xl, Xu, basic, maxiter, iterLocal
        int retcode = optimizer.MinimizeECAM(
            dim,                              // dimension
            x0,                              // initial point (in/out)
            result_value,                    // result objective value
            moecam_callback_wrapper,         // objective function
            0,                               // linear equality constraints
            0,                               // linear inequality constraints
            70.0,                            // LC parameter (default from example)
            nullptr,                         // equality constraint matrix
            nullptr,                         // inequality constraint matrix
            nullptr,                         // equality RHS
            nullptr,                         // inequality RHS
            lower_bounds.data(),             // lower bounds
            upper_bounds.data(),             // upper bounds
            nullptr,                         // basic variables
            max_iterations < 0 ? -max_iterations : max_iterations,  // max iterations
            1                                // iterLocal parameter
        );

        return retcode;

    } catch (...) {
        return -1;
    }
}

int moecam_minimize_random_start(
    int dim,
    double *x0,
    double *result_value,
    python_callback_func callback,
    int num_objectives,
    const double *bounds_lower,
    const double *bounds_upper,
    int max_iterations
) {
    // Set up global state
    g_python_callback = callback;
    g_num_objectives = num_objectives;

    try {
        Ganso optimizer;

        // Copy bounds
        std::vector<double> lower_bounds(bounds_lower, bounds_lower + dim);
        std::vector<double> upper_bounds(bounds_upper, bounds_upper + dim);

        // Call MinimizeRandomStart
        // Parameters: dim, x0, val, function, lineq, linineq, AE, AI, RHSE, RHSI, Xl, Xu, basic, maxiter
        int retcode = optimizer.MinimizeRandomStart(
            dim,                              // dimension
            x0,                              // initial point (in/out)
            result_value,                    // result objective value
            moecam_callback_wrapper,         // objective function
            0,                               // linear equality constraints
            0,                               // linear inequality constraints
            nullptr,                         // equality constraint matrix
            nullptr,                         // inequality constraint matrix
            nullptr,                         // equality RHS
            nullptr,                         // inequality RHS
            lower_bounds.data(),             // lower bounds
            upper_bounds.data(),             // upper bounds
            nullptr,                         // basic variables
            max_iterations                   // max iterations
        );

        return retcode;

    } catch (...) {
        return -1;
    }
}

int moecam_direct_optimize(
    int dim,
    double *x0,
    double *result_value,
    python_callback_func callback,
    int num_objectives,
    const double *bounds_lower,
    const double *bounds_upper,
    int max_function_evaluations
) {
    // Set up global state
    g_python_callback = callback;
    g_num_objectives = num_objectives;

    try {
        int force_stop = 0;

        // Call direct_optimize
        direct_return_code info = direct_optimize(
            direct_callback_wrapper,         // objective function
            nullptr,                         // function data
            dim,                             // dimension
            bounds_lower,                    // lower bounds
            bounds_upper,                    // upper bounds
            x0,                              // initial point (in/out)
            result_value,                    // result objective value
            max_function_evaluations,        // max function evaluations
            500,                             // max iterations (from example)
            0,                               // start parameter
            0,                               // max time
            0,                               // magic_eps
            0,                               // magic_eps_abs
            0.0,                             // vol_tol
            -1.0,                            // sigma_tol
            &force_stop,                     // force stop flag
            DIRECT_UNKNOWN_FGLOBAL,          // f_global
            0,                               // f_global_reltol
            nullptr,                         // logfile
            DIRECT_ORIGINAL                  // algorithm variant
        );

        return (int)info;

    } catch (...) {
        return -1;
    }
}
