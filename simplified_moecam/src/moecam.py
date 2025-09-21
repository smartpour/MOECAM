"""
MOECAM - Multi-Objective Evolutionary Comparison and Analysis Module
====================================================================

A simplified Python interface to MOECAM optimization algorithms.
"""

import numpy as np
from typing import Tuple, Optional, Union

try:
    from _moecam import ffi, lib
    CFFI_AVAILABLE = True
except ImportError:
    CFFI_AVAILABLE = False
    ffi = None
    lib = None

__version__ = "1.0.0"
__all__ = ['minimize_ecam', 'minimize_dfbm', 'minimize_random_start',
           'extract_pareto_front', 'extract_pareto_indices',
           'evaluate_zdt1', 'evaluate_zdt2', 'evaluate_zdt3',
           'wfg_hypervolume', 'hello_moecam', 'CFFI_AVAILABLE']


def hello_moecam() -> str:
    """Return a hello message from MOECAM."""
    if not CFFI_AVAILABLE:
        return "❌ MOECAM CFFI library not available"

    lib.print_message(b"Hello from MOECAM!")
    return "✅ MOECAM CFFI library loaded successfully"


def minimize_ecam(x0: np.ndarray, bounds: np.ndarray,
                 lipschitz_constant: float = 1.0,
                 max_iterations: int = 1000,
                 local_iterations: int = 100) -> Tuple[np.ndarray, float, int]:
    """
    Minimize a function using the ECAM algorithm.

    Parameters:
    -----------
    x0 : np.ndarray
        Initial solution vector
    bounds : np.ndarray
        Array of shape (n_vars, 2) with lower and upper bounds
    lipschitz_constant : float
        Lipschitz constant for the function
    max_iterations : int
        Maximum number of iterations
    local_iterations : int
        Number of local iterations

    Returns:
    --------
    Tuple[np.ndarray, float, int]
        Optimized solution, objective value, status code
    """
    if not CFFI_AVAILABLE:
        raise RuntimeError("MOECAM CFFI library not available")

    x0 = np.asarray(x0, dtype=np.float64)
    bounds = np.asarray(bounds, dtype=np.float64)

    n_vars = len(x0)
    xl = bounds[:, 0].copy()
    xu = bounds[:, 1].copy()

    # Create C arrays
    x0_c = ffi.new("double[]", x0.tolist())
    xl_c = ffi.new("double[]", xl.tolist())
    xu_c = ffi.new("double[]", xu.tolist())
    val_c = ffi.new("double *")

    # Call CFFI function
    status = lib.minimize_ecam(n_vars, x0_c, val_c, lipschitz_constant,
                              xl_c, xu_c, max_iterations, local_iterations)

    return np.array([x0_c[i] for i in range(n_vars)]), val_c[0], status


def minimize_dfbm(x0: np.ndarray, bounds: np.ndarray,
                 max_iterations: int = 1000) -> Tuple[np.ndarray, float, int]:
    """
    Minimize a function using the DFBM algorithm.

    Parameters:
    -----------
    x0 : np.ndarray
        Initial solution vector
    bounds : np.ndarray
        Array of shape (n_vars, 2) with lower and upper bounds
    max_iterations : int
        Maximum number of iterations

    Returns:
    --------
    Tuple[np.ndarray, float, int]
        Optimized solution, objective value, status code
    """
    if not CFFI_AVAILABLE:
        raise RuntimeError("MOECAM CFFI library not available")

    x0 = np.asarray(x0, dtype=np.float64)
    bounds = np.asarray(bounds, dtype=np.float64)

    n_vars = len(x0)
    xl = bounds[:, 0].copy()
    xu = bounds[:, 1].copy()

    # Create C arrays
    x0_c = ffi.new("double[]", x0.tolist())
    xl_c = ffi.new("double[]", xl.tolist())
    xu_c = ffi.new("double[]", xu.tolist())
    val_c = ffi.new("double *")

    # Call CFFI function
    status = lib.minimize_dfbm(n_vars, x0_c, val_c, xl_c, xu_c, max_iterations)

    return np.array([x0_c[i] for i in range(n_vars)]), val_c[0], status


def minimize_random_start(x0: np.ndarray, bounds: np.ndarray,
                         max_iterations: int = 1000) -> Tuple[np.ndarray, float, int]:
    """
    Minimize a function using random starts.

    Parameters:
    -----------
    x0 : np.ndarray
        Initial solution vector
    bounds : np.ndarray
        Array of shape (n_vars, 2) with lower and upper bounds
    max_iterations : int
        Maximum number of iterations

    Returns:
    --------
    Tuple[np.ndarray, float, int]
        Optimized solution, objective value, status code
    """
    if not CFFI_AVAILABLE:
        raise RuntimeError("MOECAM CFFI library not available")

    x0 = np.asarray(x0, dtype=np.float64)
    bounds = np.asarray(bounds, dtype=np.float64)

    n_vars = len(x0)
    xl = bounds[:, 0].copy()
    xu = bounds[:, 1].copy()

    # Create C arrays
    x0_c = ffi.new("double[]", x0.tolist())
    xl_c = ffi.new("double[]", xl.tolist())
    xu_c = ffi.new("double[]", xu.tolist())
    val_c = ffi.new("double *")

    # Call CFFI function
    status = lib.minimize_random_start(n_vars, x0_c, val_c, xl_c, xu_c, max_iterations)

    return np.array([x0_c[i] for i in range(n_vars)]), val_c[0], status


def extract_pareto_front(points: np.ndarray, strict_mode: bool = False) -> np.ndarray:
    """
    Extract Pareto front from a set of points.

    Parameters:
    -----------
    points : np.ndarray
        Array of shape (n_points, n_objectives) containing objective values
    strict_mode : bool
        If True, use strict dominance (exclude equal points)

    Returns:
    --------
    np.ndarray
        Pareto optimal points
    """
    if not CFFI_AVAILABLE:
        raise RuntimeError("MOECAM CFFI library not available")

    points = np.asarray(points, dtype=np.float64)
    n_points, n_objectives = points.shape

    # Flatten points for C interface
    points_flat = points.flatten()
    points_c = ffi.new("double[]", points_flat.tolist())

    # Allocate output array (worst case: all points are Pareto optimal)
    pareto_points_c = ffi.new("double[]", n_points * n_objectives)

    # Call CFFI function
    n_pareto = lib.extract_pareto_front(points_c, n_points, n_objectives,
                                       pareto_points_c, n_points, int(strict_mode))

    # Convert back to numpy array
    if n_pareto > 0:
        pareto_flat = [pareto_points_c[i] for i in range(n_pareto * n_objectives)]
        return np.array(pareto_flat).reshape(n_pareto, n_objectives)
    else:
        return np.empty((0, n_objectives))


def extract_pareto_indices(points: np.ndarray, strict_mode: bool = False) -> np.ndarray:
    """
    Extract indices of Pareto optimal points.

    Parameters:
    -----------
    points : np.ndarray
        Array of shape (n_points, n_objectives) containing objective values
    strict_mode : bool
        If True, use strict dominance (exclude equal points)

    Returns:
    --------
    np.ndarray
        Indices of Pareto optimal points
    """
    if not CFFI_AVAILABLE:
        raise RuntimeError("MOECAM CFFI library not available")

    points = np.asarray(points, dtype=np.float64)
    n_points, n_objectives = points.shape

    # Flatten points for C interface
    points_flat = points.flatten()
    points_c = ffi.new("double[]", points_flat.tolist())

    # Allocate output array
    indices_c = ffi.new("int[]", n_points)

    # Call CFFI function
    n_pareto = lib.extract_pareto_indices(points_c, n_points, n_objectives,
                                         indices_c, n_points, int(strict_mode))

    # Convert back to numpy array
    if n_pareto > 0:
        return np.array([indices_c[i] for i in range(n_pareto)])
    else:
        return np.empty(0, dtype=int)


def evaluate_zdt1(x: np.ndarray) -> np.ndarray:
    """Evaluate ZDT1 test function."""
    if not CFFI_AVAILABLE:
        raise RuntimeError("MOECAM CFFI library not available")

    x = np.asarray(x, dtype=np.float64)
    x_c = ffi.new("double[]", x.tolist())
    objectives_c = ffi.new("double[]", 2)

    lib.evaluate_zdt1(x_c, len(x), objectives_c)

    return np.array([objectives_c[0], objectives_c[1]])


def evaluate_zdt2(x: np.ndarray) -> np.ndarray:
    """Evaluate ZDT2 test function."""
    if not CFFI_AVAILABLE:
        raise RuntimeError("MOECAM CFFI library not available")

    x = np.asarray(x, dtype=np.float64)
    x_c = ffi.new("double[]", x.tolist())
    objectives_c = ffi.new("double[]", 2)

    lib.evaluate_zdt2(x_c, len(x), objectives_c)

    return np.array([objectives_c[0], objectives_c[1]])


def evaluate_zdt3(x: np.ndarray) -> np.ndarray:
    """Evaluate ZDT3 test function."""
    if not CFFI_AVAILABLE:
        raise RuntimeError("MOECAM CFFI library not available")

    x = np.asarray(x, dtype=np.float64)
    x_c = ffi.new("double[]", x.tolist())
    objectives_c = ffi.new("double[]", 2)

    lib.evaluate_zdt3(x_c, len(x), objectives_c)

    return np.array([objectives_c[0], objectives_c[1]])


def wfg_hypervolume(points: np.ndarray, reference: np.ndarray) -> float:
    """
    Calculate hypervolume using WFG algorithm.

    Parameters:
    -----------
    points : np.ndarray
        Pareto front points of shape (n_points, n_objectives)
    reference : np.ndarray
        Reference point for hypervolume calculation

    Returns:
    --------
    float
        Hypervolume value
    """
    if not CFFI_AVAILABLE:
        raise RuntimeError("MOECAM CFFI library not available")

    points = np.asarray(points, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)

    n_points, n_objectives = points.shape

    # Flatten points for C interface
    points_flat = points.flatten()
    points_c = ffi.new("double[]", points_flat.tolist())
    reference_c = ffi.new("double[]", reference.tolist())

    # Call CFFI function
    hv = lib.wfg_hypervolume(points_c, reference_c, n_points, n_objectives)

    return float(hv)


# Module initialization
if CFFI_AVAILABLE:
    lib.wfg_prepare(10, 1000)  # Initialize WFG with reasonable defaults
