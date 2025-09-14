#ifndef MOECAM_PYTHON_INTERFACE_H
#define MOECAM_PYTHON_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

// Python callback function type (standard Python MOO callback)
// Returns array of objective values for given input vector
typedef void (*python_callback_func)(int n, const double *x, double *objectives, int num_objectives);

// Initialize MOECAM system
int moecam_init();

// Cleanup MOECAM system
void moecam_cleanup();

// MinimizeECAM wrapper
// dim: problem dimension
// x0: initial point (input/output - final solution)
// bounds_lower, bounds_upper: variable bounds
// callback: Python objective function
// num_objectives: number of objectives
// max_iterations: maximum iterations (negative for internal default)
// Returns: 0 on success, error code otherwise
int moecam_minimize_ecam(
    int dim,
    double *x0,
    double *result_value,
    python_callback_func callback,
    int num_objectives,
    const double *bounds_lower,
    const double *bounds_upper,
    int max_iterations
);

// MinimizeRandomStart wrapper
int moecam_minimize_random_start(
    int dim,
    double *x0,
    double *result_value,
    python_callback_func callback,
    int num_objectives,
    const double *bounds_lower,
    const double *bounds_upper,
    int max_iterations
);

// DIRECT optimization wrapper
// Uses single-objective optimization (first objective)
int moecam_direct_optimize(
    int dim,
    double *x0,
    double *result_value,
    python_callback_func callback,
    int num_objectives,
    const double *bounds_lower,
    const double *bounds_upper,
    int max_function_evaluations
);

#ifdef __cplusplus
}
#endif

#endif // MOECAM_PYTHON_INTERFACE_H
