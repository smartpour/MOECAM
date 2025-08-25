"""
Multi-Objective Evolutionary Comparison and Analysis Module (MOECAM)
"""

__version__ = "0.1.0"
__author__ = "Kartik Sangwan"
__email__ = "s224247523@deakin.edu.in"

# Import main modules
try:
    from . import algorithms
    from . import problems
    from . import metrics
    from . import core
    from . import utils
except ImportError:
    # Handle missing modules gracefully during development
    pass

__all__ = [
    "algorithms",
    "problems",
    "metrics",
    "core",
    "utils"
]