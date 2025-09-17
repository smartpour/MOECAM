"""
MOECAM Algorithms CFFI Interface
===============================

Proper Python-to-C++ interface for MOECAM optimization algorithms
using CFFI like the toy example.
"""

import cffi
import os
import numpy as np
from pathlib import Path

# Setup CFFI
ffi = cffi.FFI()

# Find the compiled library
_this_dir = Path(__file__).parent
_lib_path = _this_dir / 'libmoecam.so'

# Define the C interface
ffi.cdef("""
// Callback type for objective function evaluation
typedef void (*objective_function_t)(double* x, int n_dim, double* f, int n_obj);

// MOECAM algorithm configuration
typedef struct {
    int algorithm_id;
    int n_dim;
    int n_obj;
    int max_iterations;
    double* lower_bounds;
    double* upper_bounds;
    objective_function_t objective_func;
} moecam_config_t;

// MOECAM algorithm interface
int moecam_init(moecam_config_t* config);
int moecam_run_ecam(int config_id, double* best_solution, double* best_objectives,
                   double* all_solutions, double* all_objectives, int* num_evaluations);
int moecam_run_random_start(int config_id, double* best_solution, double* best_objectives,
                           double* all_solutions, double* all_objectives, int* num_evaluations);
int moecam_run_direct(int config_id, double* best_solution, double* best_objective, int* num_evaluations);
void moecam_cleanup(int config_id);
""")

# Load the library (will be None if not built yet)
try:
    lib = ffi.dlopen(str(_lib_path))
except OSError:
    lib = None

class MOECAMAlgorithmsCFFI:
    """CFFI-based MOECAM algorithms interface."""

    def __init__(self):
        if lib is None:
            raise RuntimeError(
                f"MOECAM library not found at {_lib_path}. "
                f"Run 'python build_libraries.py' to build it."
            )
        self.configs = {}
        self.next_id = 0

    def _create_callback(self, python_func, n_dim, n_obj):
        """Create CFFI callback from Python function."""
        @ffi.callback("void(double*, int, double*, int)")
        def c_callback(x_ptr, n_dim_c, f_ptr, n_obj_c):
            # Use the parameters passed from C, not the outer scope
            actual_n_dim = int(n_dim_c)
            actual_n_obj = int(n_obj_c)

            # Extract input array
            x_arr = np.frombuffer(ffi.buffer(x_ptr, actual_n_dim * 8), dtype=np.float64).copy()

            try:
                # Call Python function
                result = python_func(x_arr)
                result_arr = np.asarray(result, dtype=np.float64)

                if result_arr.shape[0] != actual_n_obj:
                    raise ValueError(f"Expected {actual_n_obj} objectives, got {result_arr.shape[0]}")

                # Write results back
                for i in range(actual_n_obj):
                    f_ptr[i] = float(result_arr[i])

            except Exception as e:
                print(f"Objective function error: {e}")
                # Fill with large values on error
                for i in range(actual_n_obj):
                    f_ptr[i] = 1e6

        return c_callback

    def minimize_ecam(self, objective_func, bounds, num_objectives=2, max_iterations=100):
        """Run ECAM algorithm using CFFI interface.

        Args:
            objective_func: Python function f(x) -> [f1, f2, ...]
            bounds: List of (lower, upper) bounds for each dimension
            num_objectives: Number of objectives
            max_iterations: Maximum function evaluations

        Returns:
            dict: Results with success, x, objectives, all_solutions, all_objectives, etc.
        """
        bounds = np.asarray(bounds, dtype=np.float64)
        n_dim = len(bounds)

        # Prepare bounds arrays
        lower_bounds = np.ascontiguousarray(bounds[:, 0], dtype=np.float64)
        upper_bounds = np.ascontiguousarray(bounds[:, 1], dtype=np.float64)

        # Create callback
        callback = self._create_callback(objective_func, n_dim, num_objectives)

        # Create configuration
        config = ffi.new("moecam_config_t *")
        config.algorithm_id = 1  # ECAM
        config.n_dim = n_dim
        config.n_obj = num_objectives
        config.max_iterations = max_iterations
        config.lower_bounds = ffi.cast("double *", lower_bounds.ctypes.data)
        config.upper_bounds = ffi.cast("double *", upper_bounds.ctypes.data)
        config.objective_func = callback

        # Initialize
        config_id = lib.moecam_init(config)
        if config_id < 0:
            raise RuntimeError("Failed to initialize MOECAM")

        # Store callback to prevent garbage collection
        self.configs[config_id] = callback

        try:
            # Prepare output arrays - make sure we have enough space
            # ECAM starts with population of 20, so we need at least that much space
            max_evals = max(max_iterations, 50)  # Ensure minimum space for population
            best_solution = np.zeros(n_dim, dtype=np.float64)
            best_objectives = np.zeros(num_objectives, dtype=np.float64)
            all_solutions = np.zeros(max_evals * n_dim, dtype=np.float64)
            all_objectives = np.zeros(max_evals * num_objectives, dtype=np.float64)
            num_evaluations = ffi.new("int *")

            # Convert to CFFI pointers
            best_sol_ptr = ffi.cast("double *", best_solution.ctypes.data)
            best_obj_ptr = ffi.cast("double *", best_objectives.ctypes.data)
            all_sol_ptr = ffi.cast("double *", all_solutions.ctypes.data)
            all_obj_ptr = ffi.cast("double *", all_objectives.ctypes.data)

            # Run algorithm
            result = lib.moecam_run_ecam(config_id, best_sol_ptr, best_obj_ptr,
                                        all_sol_ptr, all_obj_ptr, num_evaluations)

            if result != 0:
                raise RuntimeError("ECAM algorithm failed")

            # Extract results
            actual_evals = num_evaluations[0]

            return {
                'success': True,
                'x': best_solution,
                'f': best_objectives,
                'function_evaluations': actual_evals,
                'all_solutions': all_solutions[:actual_evals * n_dim].reshape(actual_evals, n_dim),
                'all_objectives': all_objectives[:actual_evals * num_objectives].reshape(actual_evals, num_objectives),
                'message': 'ECAM optimization completed successfully'
            }

        finally:
            # Cleanup
            lib.moecam_cleanup(config_id)
            if config_id in self.configs:
                del self.configs[config_id]

    def minimize_random_start(self, objective_func, bounds, num_objectives=2, max_iterations=100):
        """Run Random Start algorithm using CFFI interface."""
        bounds = np.asarray(bounds, dtype=np.float64)
        n_dim = len(bounds)

        # Similar implementation to minimize_ecam but calls moecam_run_random_start
        lower_bounds = np.ascontiguousarray(bounds[:, 0], dtype=np.float64)
        upper_bounds = np.ascontiguousarray(bounds[:, 1], dtype=np.float64)

        callback = self._create_callback(objective_func, n_dim, num_objectives)

        config = ffi.new("moecam_config_t *")
        config.algorithm_id = 2  # Random Start
        config.n_dim = n_dim
        config.n_obj = num_objectives
        config.max_iterations = max_iterations
        config.lower_bounds = ffi.cast("double *", lower_bounds.ctypes.data)
        config.upper_bounds = ffi.cast("double *", upper_bounds.ctypes.data)
        config.objective_func = callback

        config_id = lib.moecam_init(config)
        if config_id < 0:
            raise RuntimeError("Failed to initialize Random Start")

        self.configs[config_id] = callback

        try:
            # Prepare output arrays - Random Start doesn't need as much space as ECAM
            max_evals = max_iterations  # Random Start evaluates one at a time
            best_solution = np.zeros(n_dim, dtype=np.float64)
            best_objectives = np.zeros(num_objectives, dtype=np.float64)
            all_solutions = np.zeros(max_evals * n_dim, dtype=np.float64)
            all_objectives = np.zeros(max_evals * num_objectives, dtype=np.float64)
            num_evaluations = ffi.new("int *")

            best_sol_ptr = ffi.cast("double *", best_solution.ctypes.data)
            best_obj_ptr = ffi.cast("double *", best_objectives.ctypes.data)
            all_sol_ptr = ffi.cast("double *", all_solutions.ctypes.data)
            all_obj_ptr = ffi.cast("double *", all_objectives.ctypes.data)

            result = lib.moecam_run_random_start(config_id, best_sol_ptr, best_obj_ptr,
                                               all_sol_ptr, all_obj_ptr, num_evaluations)

            if result != 0:
                raise RuntimeError("Random Start algorithm failed")

            actual_evals = num_evaluations[0]

            return {
                'success': True,
                'x': best_solution,
                'f': best_objectives,
                'function_evaluations': actual_evals,
                'all_solutions': all_solutions[:actual_evals * n_dim].reshape(actual_evals, n_dim),
                'all_objectives': all_objectives[:actual_evals * num_objectives].reshape(actual_evals, num_objectives),
                'message': 'Random Start optimization completed successfully'
            }

        finally:
            lib.moecam_cleanup(config_id)
            if config_id in self.configs:
                del self.configs[config_id]

    def direct_optimize(self, objective_func, bounds, max_iterations=100):
        """Run DIRECT algorithm using CFFI interface."""
        bounds = np.asarray(bounds, dtype=np.float64)
        n_dim = len(bounds)

        lower_bounds = np.ascontiguousarray(bounds[:, 0], dtype=np.float64)
        upper_bounds = np.ascontiguousarray(bounds[:, 1], dtype=np.float64)

        callback = self._create_callback(objective_func, n_dim, 1)  # Single objective

        config = ffi.new("moecam_config_t *")
        config.algorithm_id = 3  # DIRECT
        config.n_dim = n_dim
        config.n_obj = 1
        config.max_iterations = max_iterations
        config.lower_bounds = ffi.cast("double *", lower_bounds.ctypes.data)
        config.upper_bounds = ffi.cast("double *", upper_bounds.ctypes.data)
        config.objective_func = callback

        config_id = lib.moecam_init(config)
        if config_id < 0:
            raise RuntimeError("Failed to initialize DIRECT")

        self.configs[config_id] = callback

        try:
            best_solution = np.zeros(n_dim, dtype=np.float64)
            best_objective = ffi.new("double *")
            num_evaluations = ffi.new("int *")

            best_sol_ptr = ffi.cast("double *", best_solution.ctypes.data)

            result = lib.moecam_run_direct(config_id, best_sol_ptr, best_objective, num_evaluations)

            if result != 0:
                raise RuntimeError("DIRECT algorithm failed")

            return {
                'success': True,
                'x': best_solution,
                'fun': best_objective[0],
                'function_evaluations': num_evaluations[0],
                'message': 'DIRECT optimization completed successfully'
            }

        finally:
            lib.moecam_cleanup(config_id)
            if config_id in self.configs:
                del self.configs[config_id]

# Convenience functions
def minimize_ecam_cffi(objective_func, bounds, num_objectives=2, max_iterations=100):
    """Run ECAM algorithm using CFFI interface."""
    algorithms = MOECAMAlgorithmsCFFI()
    return algorithms.minimize_ecam(objective_func, bounds, num_objectives, max_iterations)

def minimize_random_start_cffi(objective_func, bounds, num_objectives=2, max_iterations=100):
    """Run Random Start algorithm using CFFI interface."""
    algorithms = MOECAMAlgorithmsCFFI()
    return algorithms.minimize_random_start(objective_func, bounds, num_objectives, max_iterations)

def direct_optimize_cffi(objective_func, bounds, max_iterations=100):
    """Run DIRECT algorithm using CFFI interface."""
    algorithms = MOECAMAlgorithmsCFFI()
    return algorithms.direct_optimize(objective_func, bounds, max_iterations)
