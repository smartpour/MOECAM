The MOECAM project has been successfully implemented and tested. The current architecture and implementation details are outlined below. Please refer to the `PROJECT_SUMMARY.md` for a detailed report on the project's progress and future plans.

# MOECAM Project Architecture

## 1. Project Structure

```
MOECAM/
├── src/
│   ├── moecam/
│   │   ├── __init__.py
│   │   ├── core/
│   │   │   ├── __init__.py
│   │   │   ├── cffi_interface.py  # CFFI definitions and loading C++ library
│   │   │   └── cpp_bindings.py    # Python classes wrapping C++ objects
│   │   ├── algorithms/
│   │   │   ├── __init__.py
│   │   │   └── moea_algorithms.py  # Python implementations of MOEAs
│   │   ├── problems/
│   │   │   ├── __init__.py
│     │   ├── test_functions.py  # Multi-objective test functions
│   │   │   ├── scalable_test_problems.py # Scalable test problem generators
│   │   │   ├── wfg_test_suite.py  # WFG test suite
│   │   │   └── constraints.py     # Constraint handling
│   │   ├── metrics/
│   │   │   ├── __init__.py
│   │   │   └── performance_metrics.py # Pareto front, hypervolume, etc.
│   │   └── utils/
│   │       ├── __init__.py
│   │       └── visualization.py   # Plotting Pareto fronts, etc.
├── tests/
│   ├── __init__.py
│   ├── test_cffi_interface.py
│   ├── test_algorithms.py
│   ├── test_problems.py
│   └── test_metrics.py
├── docs/
│   ├── conf.py
│   ├── index.rst
│   └── ...
├── examples/
│   ├── basic_usage.py
│   └── custom_problem.py
├── setup.py
├── README.md
└── requirements.txt
```

## 2. C++ to Python Interface (CFFI)

- **`cffi_interface.py`**: This module handles the CFFI definitions. It defines the C structures and function signatures that correspond to the C++ classes and methods. It is also responsible for loading the compiled C++ library.
- **`cpp_bindings.py`**: This module contains Python classes that act as wrappers around the C++ objects. Each Python class manages a pointer to its corresponding C++ instance and provides methods that call the C functions defined in `cffi_interface.py`. This abstracts away the CFFI details from the user.

### Key Considerations for CFFI:
- **Data Conversion**: Handles conversion between Python types (especially NumPy arrays for efficiency) and C++ types (e.g., `std::vector` or `double*`).
- **Memory Management**: Implements proper memory deallocation in C++ and ensures Python wrappers handle object lifetimes correctly to prevent memory leaks.
- **Callback Functions**: Provides a mechanism for C++ code to call back into Python functions (e.g., for objective function evaluations).

## 3. Core Components

- **`algorithms/moea_algorithms.py`**: This module contains Python implementations of various multi-objective evolutionary algorithms (MOEAs). These algorithms interact with the test functions defined in `problems/test_functions.py`.
- **`problems/test_functions.py`**: This module houses the implementations of standard multi-objective test functions (e.g., ZDT, DTLZ, WFG).
- **`metrics/performance_metrics.py`**: This module implements performance indicators such as Pareto front calculation, hypervolume, generational distance, etc. It also includes the C-code wrappers for hypervolume calculation as specified.

## 4. Testing Framework

- The `tests/` directory contains unit tests for each module. These tests ensure the correctness of the CFFI interface, algorithm implementations, and metric calculations. Pytest is used for testing.

## 5. Documentation

- The `docs/` directory contains reStructuredText files for Sphinx documentation. This includes API documentation, user guides, and examples.


```json
{
  

