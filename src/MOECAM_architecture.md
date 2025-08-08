# MOECAM Python Package Architecture

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
│   │   │   └── test_functions.py  # Multi-objective test functions
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

- **`cffi_interface.py`**: This module will handle the CFFI definitions. It will define the C structures and function signatures that correspond to the C++ classes and methods. It will also be responsible for loading the compiled C++ library.
- **`cpp_bindings.py`**: This module will contain Python classes that act as wrappers around the C++ objects. Each Python class will manage a pointer to its corresponding C++ instance and provide methods that call the C functions defined in `cffi_interface.py`. This will abstract away the CFFI details from the user.

### Key Considerations for CFFI:
- **Data Conversion**: Handle conversion between Python types (especially NumPy arrays for efficiency) and C++ types (e.g., `std::vector` or `double*`).
- **Memory Management**: Implement proper memory deallocation in C++ and ensure Python wrappers handle object lifetimes correctly to prevent memory leaks.
- **Callback Functions**: Design a mechanism for C++ code to call back into Python functions (e.g., for objective function evaluations).

## 3. Core Components

- **`algorithms/moea_algorithms.py`**: This will contain Python implementations of various multi-objective evolutionary algorithms (MOEAs). These algorithms will interact with the test functions defined in `problems/test_functions.py`.
- **`problems/test_functions.py`**: This module will house the implementations of standard multi-objective test functions (e.g., ZDT, DTLZ, WFG series) as well as the simpler ones mentioned in the requirements (e.g., from Wikipedia).
- **`metrics/performance_metrics.py`**: This module will implement performance indicators such as Pareto front calculation, hypervolume, generational distance, etc. It will also include the C-code wrappers for hypervolume calculation as specified.
- **`utils/visualization.py`**: This module will provide utilities for visualizing Pareto fronts, objective space plots, and other relevant data.

## 4. Testing Framework

- The `tests/` directory will contain unit tests for each module, ensuring the correctness of the CFFI interface, algorithm implementations, and metric calculations. `pytest` will be used for testing.

## 5. Documentation

- The `docs/` directory will contain reStructuredText files for Sphinx documentation generation. This will include API documentation generated from docstrings, a user manual, and examples.

## 6. Examples

- The `examples/` directory will provide clear and concise examples demonstrating how to use the MOECAM package for various tasks, including running algorithms on test problems and visualizing results.

## 7. Setup and Distribution

- `setup.py`: Standard Python setup file for package installation and distribution on PyPI.
- `requirements.txt`: Lists all Python dependencies.

## 8. C++ Code Integration

- The C++ numerical algorithms will be compiled into a shared library (e.g., `.so` on Linux, `.dylib` on macOS, `.dll` on Windows) that `cffi_interface.py` will load. The C++ code will expose a C-compatible API for CFFI to interact with. This will involve writing C wrappers around the C++ classes and methods, as CFFI directly interfaces with C, not C++ objects. These C wrappers will manage C++ object instances (e.g., by returning opaque pointers to Python) and handle method calls.



## 9. Class Hierarchies and Implementation Roadmap

### 9.1. Class Hierarchies

- **`moecam.core.cpp_bindings.CppOptimizer`**: Base class for wrapping C++ optimizer instances.
  - `__init__(self, cpp_optimizer_ptr)`: Initializes with a pointer to the C++ object.
  - `run_optimization(self, problem, ...)`: Calls the C++ optimization routine.
  - `get_results(self)`: Retrieves results from the C++ object.
  - `cleanup(self)`: Deallocates C++ memory.

- **`moecam.problems.test_functions.MultiObjectiveProblem`**: Abstract base class for multi-objective problems.
  - `evaluate(self, x)`: Abstract method to evaluate objective functions.
  - `get_bounds(self)`: Returns variable bounds.
  - `get_num_objectives(self)`: Returns number of objectives.

- **`moecam.algorithms.moea_algorithms.MOEA`**: Abstract base class for MOEAs.
  - `optimize(self, problem)`: Abstract method to run the optimization.

- **`moecam.metrics.performance_metrics.Metric`**: Abstract base class for performance metrics.
  - `calculate(self, pareto_front_approximation, true_pareto_front)`: Abstract method to calculate metric.

### 9.2. Implementation Roadmap

1.  **Set up basic project structure**: Create directories and initial `__init__.py` files.
2.  **CFFI Toy Example**: Create a simple C++ library with a few functions and a corresponding CFFI Python wrapper to demonstrate basic data passing and function calls. This will validate the CFFI setup.
3.  **Implement Basic Test Functions**: Code a few simple multi-objective test functions in `moecam/problems/test_functions.py` (e.g., ZDT1, ZDT2 from Wikipedia).
4.  **Develop C++ Wrapper for a Toy Class**: Create a simple C++ class with member variables and methods, and then write the C-compatible wrapper functions and the Python `CppOptimizer` class using CFFI. Include memory management and a basic callback mechanism.
5.  **Integrate Test Functions with C++ Wrapper**: Modify the C++ wrapper to accept a Python callback for objective function evaluation, allowing C++ algorithms to use Python-defined test functions.
6.  **Implement a Simple MOEA in Python**: Start with a basic MOEA (e.g., a simple NSGA-II or MOEA/D) in `moecam/algorithms/moea_algorithms.py` that uses the Python-wrapped C++ objective function evaluation.
7.  **Implement Basic Performance Metrics**: Add Pareto front calculation and a simple diversity metric in `moecam/metrics/performance_metrics.py`.
8.  **Visualization**: Implement basic plotting for Pareto fronts in `moecam/utils/visualization.py`.
9.  **Testing**: Write unit tests for each implemented component.
10. **Documentation**: Start drafting the documentation for the implemented modules.
11. **Iterative Refinement**: Continuously add more complex algorithms, test functions, and metrics, refining the C++ interface and Python wrappers as needed.


