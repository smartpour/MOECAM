"""
MOECAM Unified CFFI Interface
============================

Proper Python-to-C++ interface for all MOECAM functionality
using CFFI compiled at installation time.
"""

import numpy as np
from typing import Callable, List, Tuple, Dict, Any, Optional
import warnings

try:
    from moecam._moecam_cffi import ffi, lib
except ImportError:
    # Fallback for development
    try:
        import cffi
        ffi = cffi.FFI()
        lib = None
        warnings.warn("MOECAM CFFI library not found. Run 'pip install -e .' to build it.")
    except ImportError:
        ffi = None
        lib = None
        warnings.warn("CFFI not available. Install cffi to use MOECAM.")


class MOECAMInterface:
    """Unified interface for all MOECAM functionality."""

    def __init__(self):
        if lib is None:
            raise RuntimeError("MOECAM CFFI library not available. Run 'pip install -e .' to build it.")
        self._callbacks = {}  # Keep callbacks alive

    def _create_callback(self, python_func: Callable, n_dim: int, n_obj: int):
        """Create CFFI callback from Python function."""
        @ffi.callback("void(const double*, int, double*, int)")
        def c_callback(x_ptr, n_dim_c, f_ptr, n_obj_c):
            try:
                # Extract input array
                x_arr = np.frombuffer(ffi.buffer(x_ptr, n_dim_c * 8), dtype=np.float64).copy()

                # Call Python function
                result = python_func(x_arr)
                result_arr = np.asarray(result, dtype=np.float64)

                if len(result_arr) != n_obj_c:
                    raise ValueError(f"Expected {n_obj_c} objectives, got {len(result_arr)}")

                # Write results back
                for i in range(n_obj_c):
                    f_ptr[i] = float(result_arr[i])

            except Exception as e:
                print(f"Objective function error: {e}")
                # Fill with large values on error
                for i in range(n_obj_c):
                    f_ptr[i] = 1e6

        return c_callback

    def optimize_ecam(self, objective_func: Callable, bounds: List[Tuple[float, float]],
                     num_objectives: int = 1, max_iterations: int = 100,
                     max_solutions: int = 1000) -> Dict[str, Any]:
        """
        Run ECAM optimization algorithm.

        Args:
            objective_func: Function f(x) -> [f1, f2, ...] returning objective values
            bounds: List of (min, max) bounds for each variable
            num_objectives: Number of objectives
            max_iterations: Maximum number of iterations
            max_solutions: Maximum solutions to store

        Returns:
            Dictionary with optimization results
        """
        n_dim = len(bounds)

        # Prepare bounds arrays
        lower_bounds = ffi.new("double[]", [b[0] for b in bounds])
        upper_bounds = ffi.new("double[]", [b[1] for b in bounds])

        # Create callback
        callback = self._create_callback(objective_func, n_dim, num_objectives)
        callback_id = id(callback)
        self._callbacks[callback_id] = callback  # Keep alive

        try:
            # Initialize MOECAM
            config_id = lib.moecam_init(n_dim, num_objectives, lower_bounds, upper_bounds, callback)
            if config_id < 0:
                raise RuntimeError("Failed to initialize MOECAM")

            # Prepare output arrays
            best_solution = ffi.new("double[]", n_dim)
            best_objectives = ffi.new("double[]", num_objectives)
            all_solutions = ffi.new("double[]", max_solutions * n_dim)
            all_objectives = ffi.new("double[]", max_solutions * num_objectives)
            num_evaluations = ffi.new("int*")

            # Run ECAM
            result = lib.moecam_run_ecam(config_id, n_dim, lower_bounds, upper_bounds,
                                       max_iterations, best_solution, best_objectives,
                                       all_solutions, all_objectives, max_solutions, num_evaluations)

            if result < 0:
                raise RuntimeError("ECAM optimization failed")

            # Extract results
            best_x = np.array([best_solution[i] for i in range(n_dim)])
            best_f = np.array([best_objectives[i] for i in range(num_objectives)])
            n_evals = num_evaluations[0]

            # Extract all solutions
            all_x = np.zeros((n_evals, n_dim))
            all_f = np.zeros((n_evals, num_objectives))

            for i in range(min(n_evals, max_solutions)):
                for j in range(n_dim):
                    all_x[i, j] = all_solutions[i * n_dim + j]
                for j in range(num_objectives):
                    all_f[i, j] = all_objectives[i * num_objectives + j]

            # Clean up
            lib.moecam_cleanup(config_id)

            return {
                'success': True,
                'best_solution': best_x,
                'best_objectives': best_f,
                'all_solutions': all_x,
                'all_objectives': all_f,
                'function_evaluations': n_evals,
                'algorithm': 'ECAM'
            }

        finally:
            # Clean up callback
            if callback_id in self._callbacks:
                del self._callbacks[callback_id]

    def optimize_random_start(self, objective_func: Callable, bounds: List[Tuple[float, float]],
                             num_objectives: int = 1, max_iterations: int = 100,
                             max_solutions: int = 1000) -> Dict[str, Any]:
        """
        Run Random Start optimization algorithm.

        Args:
            objective_func: Function f(x) -> [f1, f2, ...] returning objective values
            bounds: List of (min, max) bounds for each variable
            num_objectives: Number of objectives
            max_iterations: Maximum number of iterations
            max_solutions: Maximum solutions to store

        Returns:
            Dictionary with optimization results
        """
        n_dim = len(bounds)

        # Prepare bounds arrays
        lower_bounds = ffi.new("double[]", [b[0] for b in bounds])
        upper_bounds = ffi.new("double[]", [b[1] for b in bounds])

        # Create callback
        callback = self._create_callback(objective_func, n_dim, num_objectives)
        callback_id = id(callback)
        self._callbacks[callback_id] = callback  # Keep alive

        try:
            # Initialize MOECAM
            config_id = lib.moecam_init(n_dim, num_objectives, lower_bounds, upper_bounds, callback)
            if config_id < 0:
                raise RuntimeError("Failed to initialize MOECAM")

            # Prepare output arrays
            best_solution = ffi.new("double[]", n_dim)
            best_objectives = ffi.new("double[]", num_objectives)
            all_solutions = ffi.new("double[]", max_solutions * n_dim)
            all_objectives = ffi.new("double[]", max_solutions * num_objectives)
            num_evaluations = ffi.new("int*")

            # Run Random Start
            result = lib.moecam_run_random_start(config_id, n_dim, lower_bounds, upper_bounds,
                                               max_iterations, best_solution, best_objectives,
                                               all_solutions, all_objectives, max_solutions, num_evaluations)

            if result < 0:
                raise RuntimeError("Random Start optimization failed")

            # Extract results
            best_x = np.array([best_solution[i] for i in range(n_dim)])
            best_f = np.array([best_objectives[i] for i in range(num_objectives)])
            n_evals = num_evaluations[0]

            # Extract all solutions
            all_x = np.zeros((n_evals, n_dim))
            all_f = np.zeros((n_evals, num_objectives))

            for i in range(min(n_evals, max_solutions)):
                for j in range(n_dim):
                    all_x[i, j] = all_solutions[i * n_dim + j]
                for j in range(num_objectives):
                    all_f[i, j] = all_objectives[i * num_objectives + j]

            # Clean up
            lib.moecam_cleanup(config_id)

            return {
                'success': True,
                'best_solution': best_x,
                'best_objectives': best_f,
                'all_solutions': all_x,
                'all_objectives': all_f,
                'function_evaluations': n_evals,
                'algorithm': 'Random Start'
            }

        finally:
            # Clean up callback
            if callback_id in self._callbacks:
                del self._callbacks[callback_id]

    def extract_pareto_front(self, points: np.ndarray, strict_mode: bool = False) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extract Pareto front from a set of objective vectors.

        Args:
            points: Array of shape (n_points, n_objectives)
            strict_mode: Use strict Pareto dominance if True, weak if False

        Returns:
            Tuple of (pareto_points, pareto_indices)
        """
        points = np.asarray(points, dtype=np.float64)
        if points.ndim != 2:
            raise ValueError("Points must be a 2D array")

        n_points, n_objectives = points.shape
        max_pareto = n_points  # At most all points can be Pareto optimal

        # Flatten points for C interface
        points_flat = ffi.new("double[]", points.flatten())
        pareto_points_flat = ffi.new("double[]", max_pareto * n_objectives)
        pareto_indices = ffi.new("int[]", max_pareto)

        # Extract Pareto front
        n_pareto = lib.extract_pareto_front(points_flat, n_points, n_objectives,
                                          pareto_points_flat, max_pareto, int(strict_mode))

        if n_pareto < 0:
            raise RuntimeError("Pareto front extraction failed")

        # Extract results
        pareto_points = np.zeros((n_pareto, n_objectives))
        for i in range(n_pareto):
            for j in range(n_objectives):
                pareto_points[i, j] = pareto_points_flat[i * n_objectives + j]

        # Get indices too
        n_indices = lib.extract_pareto_indices(points_flat, n_points, n_objectives,
                                             pareto_indices, max_pareto, int(strict_mode))

        indices = np.array([pareto_indices[i] for i in range(n_indices)])

        return pareto_points, indices

    def calculate_hypervolume(self, points: np.ndarray, reference_point: np.ndarray) -> float:
        """
        Calculate hypervolume of a set of points.

        Args:
            points: Array of shape (n_points, n_objectives)
            reference_point: Reference point for hypervolume calculation

        Returns:
            Hypervolume value
        """
        points = np.asarray(points, dtype=np.float64)
        reference_point = np.asarray(reference_point, dtype=np.float64)

        if points.ndim != 2:
            raise ValueError("Points must be a 2D array")

        n_points, n_objectives = points.shape

        if len(reference_point) != n_objectives:
            raise ValueError("Reference point must have same dimension as objectives")

        # Initialize WFG
        if lib.wfg_init(n_objectives, n_points) < 0:
            raise RuntimeError("Failed to initialize WFG hypervolume")

        try:
            # Prepare arrays
            points_flat = ffi.new("double[]", points.flatten())
            ref_point = ffi.new("double[]", reference_point)

            # Calculate hypervolume
            volume = lib.wfg_calculate_hypervolume(points_flat, ref_point, n_points, n_objectives)

            if volume < 0:
                raise RuntimeError("Hypervolume calculation failed")

            return float(volume)

        finally:
            lib.wfg_cleanup()


# Create singleton instance
_moecam_interface = None

def get_interface() -> MOECAMInterface:
    """Get the MOECAM interface singleton."""
    global _moecam_interface
    if _moecam_interface is None:
        _moecam_interface = MOECAMInterface()
    return _moecam_interface


# Convenience functions
def minimize_ecam_cffi(objective_func: Callable, bounds: List[Tuple[float, float]],
                      num_objectives: int = 2, max_iterations: int = 100) -> Dict[str, Any]:
    """Convenience function for ECAM optimization."""
    return get_interface().optimize_ecam(objective_func, bounds, num_objectives, max_iterations)


def minimize_random_start_cffi(objective_func: Callable, bounds: List[Tuple[float, float]],
                              num_objectives: int = 2, max_iterations: int = 100) -> Dict[str, Any]:
    """Convenience function for Random Start optimization."""
    return get_interface().optimize_random_start(objective_func, bounds, num_objectives, max_iterations)


def extract_pareto_front_cffi(points: np.ndarray, strict_mode: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """Convenience function for Pareto front extraction."""
    return get_interface().extract_pareto_front(points, strict_mode)


def calculate_hypervolume_cffi(points: np.ndarray, reference_point: np.ndarray) -> float:
    """Convenience function for hypervolume calculation."""
    return get_interface().calculate_hypervolume(points, reference_point)
