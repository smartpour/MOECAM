"""
Core functionality for MOECAM
"""

# Import core modules if they exist
try:
    from .cffi_interface import *
except ImportError:
    pass

try:
    from .cpp_bindings import *
except ImportError:
    pass

try:
    from .moecam_direct import (
        extract_pareto_front,
        calculate_hypervolume,
        WFGHypervolume
    )
except ImportError:
    pass

__all__ = [
    "cffi_interface",
    "cpp_bindings",
    "extract_pareto_front",
    "calculate_hypervolume",
    "WFGHypervolume"
]