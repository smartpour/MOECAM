import numpy as np
import os
import sys

# Ensure package path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
from moecam.core.cffi_interface import CppOptimizer  # noqa: E402

"""Demonstration of registering a Python objective callback for the C++ optimizer
and visualizing sampled objective values.

The objective callback is now REQUIRED - the optimizer cannot be created without one.

Run inside the repo (after building toy_cpp_lib.so):
    python examples/callback_objective_demo.py
"""

def objectives(x: np.ndarray) -> np.ndarray:
    # Example bi-objective: sphere and shifted sphere
    return np.array([
        np.sum(x**2),
        np.sum((x-1.0)**2)
    ], dtype=np.float64)

if __name__ == "__main__":
    n_dim = 5
    lower = -5 * np.ones(n_dim)
    upper = 5 * np.ones(n_dim)

    # Callback is now required - cannot pass None
    opt = CppOptimizer(n_dim, 2, lower, upper, objective_callback=objectives)
    print("Optimizer info:", opt.get_info())
    print("Evaluate at zeros:", opt.evaluate(np.zeros(n_dim)))

    # Plot samples
    ax = opt.plot_objectives(n_samples=200)
    import matplotlib.pyplot as plt
    plt.show()
