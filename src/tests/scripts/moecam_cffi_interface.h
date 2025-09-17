
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
