"""
MOECAM Python Interface

Provides Python bindings for the three main MOECAM optimization algorithms:
1. MinimizeECAM - Evolutionary Constrained Augmented Minimization
2. MinimizeRandomStart - Multi-start optimization
3. DirectOptimize - DIRECT global optimization algorithm

Usage:
    from moecam_interface import MOECAMOptimizer

    def objective_function(x):
        # Return list of objective values
        f1 = sum(xi**2 for xi in x)
        f2 = sum((xi-1)**2 for xi in x)
        return [f1, f2]

    optimizer = MOECAMOptimizer()

    # Example: ECAM optimization
    result = optimizer.minimize_ecam(
        objective_function,
        bounds=[(-5, 5), (-5, 5)],  # bounds for each variable
        num_objectives=2,
        max_iterations=1000
    )

    print(f"Best solution: {result['x']}")
    print(f"Best objective: {result['f']}")
"""

import ctypes
import numpy as np
import os
from typing import List, Tuple, Callable, Dict, Optional

class MOECAMOptimizer:
    """Python interface to MOECAM optimization algorithms."""

    def __init__(self, library_path: Optional[str] = None):
        """Initialize MOECAM optimizer.

        Args:
            library_path: Path to MOECAM shared library. If None, auto-detect.
        """
        if library_path is None:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            library_path = os.path.join(script_dir, 'libmoecam.so')

        if not os.path.exists(library_path):
            raise FileNotFoundError(
                f"MOECAM library not found at {library_path}. "
                f"Build the library first with the provided Makefile."
            )

        self.lib = ctypes.CDLL(library_path)
        self._setup_function_signatures()
        self._initialize()

    def _setup_function_signatures(self):
        """Set up ctypes function signatures."""

        # Callback function type
        self.CALLBACK_FUNC = ctypes.CFUNCTYPE(
            None,                           # return type
            ctypes.c_int,                   # n (dimension)
            ctypes.POINTER(ctypes.c_double), # x (input vector)
            ctypes.POINTER(ctypes.c_double), # objectives (output)
            ctypes.c_int                    # num_objectives
        )

        # moecam_init
        self.lib.moecam_init.argtypes = []
        self.lib.moecam_init.restype = ctypes.c_int

        # moecam_cleanup
        self.lib.moecam_cleanup.argtypes = []
        self.lib.moecam_cleanup.restype = None

        # moecam_minimize_ecam
        self.lib.moecam_minimize_ecam.argtypes = [
            ctypes.c_int,                    # dim
            ctypes.POINTER(ctypes.c_double), # x0
            ctypes.POINTER(ctypes.c_double), # result_value
            self.CALLBACK_FUNC,              # callback
            ctypes.c_int,                    # num_objectives
            ctypes.POINTER(ctypes.c_double), # bounds_lower
            ctypes.POINTER(ctypes.c_double), # bounds_upper
            ctypes.c_int                     # max_iterations
        ]
        self.lib.moecam_minimize_ecam.restype = ctypes.c_int

        # moecam_minimize_random_start
        self.lib.moecam_minimize_random_start.argtypes = [
            ctypes.c_int,                    # dim
            ctypes.POINTER(ctypes.c_double), # x0
            ctypes.POINTER(ctypes.c_double), # result_value
            self.CALLBACK_FUNC,              # callback
            ctypes.c_int,                    # num_objectives
            ctypes.POINTER(ctypes.c_double), # bounds_lower
            ctypes.POINTER(ctypes.c_double), # bounds_upper
            ctypes.c_int                     # max_iterations
        ]
        self.lib.moecam_minimize_random_start.restype = ctypes.c_int

        # moecam_direct_optimize
        self.lib.moecam_direct_optimize.argtypes = [
            ctypes.c_int,                    # dim
            ctypes.POINTER(ctypes.c_double), # x0
            ctypes.POINTER(ctypes.c_double), # result_value
            self.CALLBACK_FUNC,              # callback
            ctypes.c_int,                    # num_objectives
            ctypes.POINTER(ctypes.c_double), # bounds_lower
            ctypes.POINTER(ctypes.c_double), # bounds_upper
            ctypes.c_int                     # max_function_evaluations
        ]
        self.lib.moecam_direct_optimize.restype = ctypes.c_int

    def _initialize(self):
        """Initialize MOECAM system."""
        result = self.lib.moecam_init()
        if result != 0:
            raise RuntimeError(f"Failed to initialize MOECAM: error {result}")

    def _create_callback_wrapper(self, python_func: Callable[[List[float]], List[float]]):
        """Create a C callback wrapper for Python objective function."""

        def callback_wrapper(n, x_ptr, objectives_ptr, num_objectives):
            try:
                # Convert C array to Python list
                x = [x_ptr[i] for i in range(n)]

                # Call Python objective function
                objective_values = python_func(x)

                # Convert result back to C array
                for i in range(min(len(objective_values), num_objectives)):
                    objectives_ptr[i] = objective_values[i]

            except Exception as e:
                print(f"Error in objective function: {e}")
                # Set error values
                for i in range(num_objectives):
                    objectives_ptr[i] = 1e20

        return self.CALLBACK_FUNC(callback_wrapper)

    def _prepare_bounds(self, bounds: List[Tuple[float, float]]) -> Tuple[np.ndarray, np.ndarray]:
        """Convert bounds to lower and upper bound arrays."""
        lower_bounds = np.array([b[0] for b in bounds], dtype=np.float64)
        upper_bounds = np.array([b[1] for b in bounds], dtype=np.float64)
        return lower_bounds, upper_bounds

    def minimize_ecam(self,
                     objective_func: Callable[[List[float]], List[float]],
                     bounds: List[Tuple[float, float]],
                     num_objectives: int,
                     initial_point: Optional[List[float]] = None,
                     max_iterations: int = 1000) -> Dict:
        """Minimize using ECAM algorithm.

        Args:
            objective_func: Function that takes a list of variables and returns list of objectives
            bounds: List of (lower, upper) bounds for each variable
            num_objectives: Number of objectives
            initial_point: Initial guess (if None, uses center of bounds)
            max_iterations: Maximum iterations

        Returns:
            dict: {'x': solution, 'f': objective_value, 'success': bool, 'message': str}
        """
        dim = len(bounds)

        # Prepare bounds
        lower_bounds, upper_bounds = self._prepare_bounds(bounds)

        # Prepare initial point
        if initial_point is None:
            x0 = 0.5 * (lower_bounds + upper_bounds)
        else:
            x0 = np.array(initial_point, dtype=np.float64)

        # Create callback
        callback = self._create_callback_wrapper(objective_func)

        # Prepare output
        result_value = ctypes.c_double(0.0)

        # Convert to ctypes
        x0_c = (ctypes.c_double * dim)(*x0)
        lower_c = (ctypes.c_double * dim)(*lower_bounds)
        upper_c = (ctypes.c_double * dim)(*upper_bounds)

        # Call optimization
        retcode = self.lib.moecam_minimize_ecam(
            dim,
            x0_c,
            ctypes.byref(result_value),
            callback,
            num_objectives,
            lower_c,
            upper_c,
            max_iterations
        )

        # Extract results
        solution = [x0_c[i] for i in range(dim)]
        success = (retcode == 0)

        return {
            'x': solution,
            'f': result_value.value,
            'success': success,
            'message': f'ECAM optimization completed with return code {retcode}',
            'return_code': retcode
        }

    def minimize_random_start(self,
                             objective_func: Callable[[List[float]], List[float]],
                             bounds: List[Tuple[float, float]],
                             num_objectives: int,
                             initial_point: Optional[List[float]] = None,
                             max_iterations: int = 1000) -> Dict:
        """Minimize using random start algorithm.

        Args:
            objective_func: Function that takes a list of variables and returns list of objectives
            bounds: List of (lower, upper) bounds for each variable
            num_objectives: Number of objectives
            initial_point: Initial guess (if None, uses center of bounds)
            max_iterations: Maximum iterations

        Returns:
            dict: {'x': solution, 'f': objective_value, 'success': bool, 'message': str}
        """
        dim = len(bounds)

        # Prepare bounds
        lower_bounds, upper_bounds = self._prepare_bounds(bounds)

        # Prepare initial point
        if initial_point is None:
            x0 = 0.5 * (lower_bounds + upper_bounds)
        else:
            x0 = np.array(initial_point, dtype=np.float64)

        # Create callback
        callback = self._create_callback_wrapper(objective_func)

        # Prepare output
        result_value = ctypes.c_double(0.0)

        # Convert to ctypes
        x0_c = (ctypes.c_double * dim)(*x0)
        lower_c = (ctypes.c_double * dim)(*lower_bounds)
        upper_c = (ctypes.c_double * dim)(*upper_bounds)

        # Call optimization
        retcode = self.lib.moecam_minimize_random_start(
            dim,
            x0_c,
            ctypes.byref(result_value),
            callback,
            num_objectives,
            lower_c,
            upper_c,
            max_iterations
        )

        # Extract results
        solution = [x0_c[i] for i in range(dim)]
        success = (retcode == 0)

        return {
            'x': solution,
            'f': result_value.value,
            'success': success,
            'message': f'Random start optimization completed with return code {retcode}',
            'return_code': retcode
        }

    def direct_optimize(self,
                       objective_func: Callable[[List[float]], List[float]],
                       bounds: List[Tuple[float, float]],
                       num_objectives: int = 1,
                       initial_point: Optional[List[float]] = None,
                       max_function_evaluations: int = 1000) -> Dict:
        """Optimize using DIRECT algorithm (single-objective, uses first objective).

        Args:
            objective_func: Function that takes a list of variables and returns list of objectives
            bounds: List of (lower, upper) bounds for each variable
            num_objectives: Number of objectives (DIRECT uses only the first one)
            initial_point: Initial guess (if None, uses center of bounds)
            max_function_evaluations: Maximum function evaluations

        Returns:
            dict: {'x': solution, 'f': objective_value, 'success': bool, 'message': str}
        """
        dim = len(bounds)

        # Prepare bounds
        lower_bounds, upper_bounds = self._prepare_bounds(bounds)

        # Prepare initial point
        if initial_point is None:
            x0 = 0.5 * (lower_bounds + upper_bounds)
        else:
            x0 = np.array(initial_point, dtype=np.float64)

        # Create callback
        callback = self._create_callback_wrapper(objective_func)

        # Prepare output
        result_value = ctypes.c_double(0.0)

        # Convert to ctypes
        x0_c = (ctypes.c_double * dim)(*x0)
        lower_c = (ctypes.c_double * dim)(*lower_bounds)
        upper_c = (ctypes.c_double * dim)(*upper_bounds)

        # Call optimization
        retcode = self.lib.moecam_direct_optimize(
            dim,
            x0_c,
            ctypes.byref(result_value),
            callback,
            max(1, num_objectives),  # Ensure at least 1 objective
            lower_c,
            upper_c,
            max_function_evaluations
        )

        # Extract results
        solution = [x0_c[i] for i in range(dim)]

        # DIRECT return codes: positive = success, negative = error
        success = (retcode > 0)

        return {
            'x': solution,
            'f': result_value.value,
            'success': success,
            'message': f'DIRECT optimization completed with return code {retcode}',
            'return_code': retcode
        }

    def cleanup(self):
        """Clean up MOECAM resources."""
        self.lib.moecam_cleanup()

    def __del__(self):
        """Cleanup on destruction."""
        if hasattr(self, 'lib'):
            self.cleanup()

# Convenience functions for quick usage
def minimize_ecam(objective_func: Callable[[List[float]], List[float]],
                 bounds: List[Tuple[float, float]],
                 num_objectives: int = 2,
                 **kwargs) -> Dict:
    """Quick ECAM optimization."""
    optimizer = MOECAMOptimizer()
    try:
        return optimizer.minimize_ecam(objective_func, bounds, num_objectives, **kwargs)
    finally:
        optimizer.cleanup()

def minimize_random_start(objective_func: Callable[[List[float]], List[float]],
                         bounds: List[Tuple[float, float]],
                         num_objectives: int = 2,
                         **kwargs) -> Dict:
    """Quick random start optimization."""
    optimizer = MOECAMOptimizer()
    try:
        return optimizer.minimize_random_start(objective_func, bounds, num_objectives, **kwargs)
    finally:
        optimizer.cleanup()

def direct_optimize(objective_func: Callable[[List[float]], List[float]],
                   bounds: List[Tuple[float, float]],
                   **kwargs) -> Dict:
    """Quick DIRECT optimization."""
    optimizer = MOECAMOptimizer()
    try:
        return optimizer.direct_optimize(objective_func, bounds, **kwargs)
    finally:
        optimizer.cleanup()

if __name__ == "__main__":
    # Test with simple multi-objective functions
    print("Testing MOECAM Python Interface")
    print("=" * 40)

    # Define test functions
    def simple_bi_objective(x):
        """Simple bi-objective function."""
        f1 = sum(xi**2 for xi in x)
        f2 = sum((xi - 1)**2 for xi in x)
        return [f1, f2]

    def kursawe(x):
        """Kursawe function."""
        import math
        f1 = sum(-10 * math.exp(-0.2 * math.sqrt(x[i]**2 + x[i+1]**2)) for i in range(len(x)-1))
        f2 = sum(abs(xi)**0.8 + 5 * math.sin(xi**3) for xi in x)
        return [f1, f2]

    bounds = [(-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0)]

    print("\nTest 1: ECAM with simple bi-objective")
    try:
        result = minimize_ecam(simple_bi_objective, bounds, num_objectives=2, max_iterations=100)
        print(f"  Success: {result['success']}")
        print(f"  Solution: {[f'{x:.4f}' for x in result['x']]}")
        print(f"  Objective: {result['f']:.6f}")
    except Exception as e:
        print(f"  Error: {e}")

    print("\nTest 2: Random Start with Kursawe")
    try:
        result = minimize_random_start(kursawe, bounds, num_objectives=2, max_iterations=100)
        print(f"  Success: {result['success']}")
        print(f"  Solution: {[f'{x:.4f}' for x in result['x']]}")
        print(f"  Objective: {result['f']:.6f}")
    except Exception as e:
        print(f"  Error: {e}")

    print("\nTest 3: DIRECT with single objective")
    try:
        result = direct_optimize(lambda x: [sum(xi**2 for xi in x)], bounds, max_function_evaluations=100)
        print(f"  Success: {result['success']}")
        print(f"  Solution: {[f'{x:.4f}' for x in result['x']]}")
        print(f"  Objective: {result['f']:.6f}")
    except Exception as e:
        print(f"  Error: {e}")

    print("\n✅ MOECAM Python interface testing complete!")
