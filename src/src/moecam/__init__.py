"""
Multi-Objective Evolutionary Comparison and Analysis Module (MOECAM)
"""

__version__ = "0.1.0"
__author__ = "smartpour"
__email__ = "s224247523@deakin.edu.au"

# Import main modules
try:
    from . import algorithms
    from . import problems
    from . import metrics
    from . import core
    from . import utils
    from . import visualization
except ImportError:
    # Handle missing modules gracefully during development
    pass

# Import main MOECAM interfaces for easy access
try:
    from .algorithms import (
        minimize_ecam,
        minimize_random_start,
        direct_optimize,
        MOECAMDirectInterface
    )
    from .core import (
        extract_pareto_front,
        calculate_hypervolume,
        WFGHypervolume
    )
    __all__ = [
        'minimize_ecam',
        'minimize_random_start',
        'direct_optimize',
        'MOECAMDirectInterface',
        'extract_pareto_front',
        'calculate_hypervolume',
        'WFGHypervolume'
    ]
except ImportError:
    __all__ = []

__all__ = [
    "algorithms",
    "problems",
    "metrics",
    "core",
    "utils",
    "visualization"
]