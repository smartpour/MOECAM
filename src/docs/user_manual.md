# MOECAM User Manual

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Test Problems](#test-problems)
5. [Algorithms](#algorithms)
6. [Performance Metrics](#performance-metrics)
7. [C++ Integration](#cpp-integration)
8. [Examples](#examples)
9. [API Reference](#api-reference)

## Introduction

MOECAM (Multi-Objective Evolutionary Comparison and Analysis Module) is a comprehensive Python package designed for multi-objective optimization research and applications. It provides:

- Multiple state-of-the-art multi-objective evolutionary algorithms
- Standard test problem suites (ZDT, DTLZ, WFG)
- Performance evaluation metrics
- Visualization tools
- C++ integration capabilities via CFFI

## Installation

### From PyPI (Recommended)

```bash
pip install moecam
```

### From Source

```bash
git clone https://github.com/smartpour/moecam.git
cd moecam
pip install -e .
```

### Development Installation

```bash
git clone https://github.com/smartpour/moecam.git
cd moecam
pip install -e .[dev]
```

## Quick Start

Here's a simple example to get you started:

```python
from moecam.problems.test_functions import ZDT1
from moecam.algorithms.moea_algorithms import NSGAII
from moecam.metrics.performance_metrics import hypervolume

# Define problem
problem = ZDT1(n_dim=30)

# Initialize algorithm
algorithm = NSGAII(problem, pop_size=100, num_generations=250)

# Run optimization
pareto_front = algorithm.optimize()

# Evaluate performance
reference_point = [1.1, 1.1]
hv = hypervolume(pareto_front, reference_point)
print(f"Hypervolume: {hv}")
```

## Test Problems

MOECAM includes several standard test problem suites:

### ZDT Test Suite

The ZDT test suite consists of six test functions (ZDT1-ZDT6) with known Pareto fronts:

```python
from moecam.problems.test_functions import ZDT1, ZDT2, ZDT3

# ZDT1: Convex Pareto front
problem1 = ZDT1(n_dim=30)

# ZDT2: Non-convex Pareto front
problem2 = ZDT2(n_dim=30)

# ZDT3: Disconnected Pareto front
problem3 = ZDT3(n_dim=30)
```

### DTLZ Test Suite

The DTLZ test suite is designed for many-objective optimization:

```python
from moecam.problems.test_functions import DTLZ1, DTLZ2

# DTLZ1: Linear Pareto front
problem1 = DTLZ1(n_dim=7, n_obj=3)

# DTLZ2: Spherical Pareto front
problem2 = DTLZ2(n_dim=12, n_obj=3)
```

### Custom Problems

You can define custom problems by inheriting from the base problem class:

```python
import numpy as np

class CustomProblem:
    def __init__(self, n_dim=10):
        self.n_dim = n_dim
        self.n_obj = 2
        self.xl = np.zeros(n_dim)
        self.xu = np.ones(n_dim)

    def evaluate(self, x):
        f1 = np.sum(x**2)
        f2 = np.sum((x - 1)**2)
        return np.array([f1, f2])
```

## Algorithms

### NSGA-II

The Non-dominated Sorting Genetic Algorithm II is one of the most popular MOEAs:

```python
from moecam.algorithms.moea_algorithms import NSGAII

algorithm = NSGAII(
    problem=problem,
    pop_size=100,
    num_generations=250,
    crossover_rate=0.9,
    mutation_rate=0.01
)

pareto_front = algorithm.optimize()
```

### MOEA/D

Multi-Objective Evolutionary Algorithm based on Decomposition:

```python
from moecam.algorithms.moea_algorithms import MOEAD

algorithm = MOEAD(
    problem=problem,
    pop_size=100,
    num_generations=250
)

pareto_front = algorithm.optimize()
```

## Performance Metrics

MOECAM provides several performance metrics for evaluating algorithm performance:

### Hypervolume

The hypervolume indicator measures the volume of objective space dominated by a solution set:

```python
from moecam.metrics.performance_metrics import hypervolume

reference_point = [1.1, 1.1]  # Must dominate all Pareto optimal solutions
hv = hypervolume(pareto_front, reference_point)
```

### Pareto Front Extraction

Extract the non-dominated solutions from a set of objective vectors:

```python
from moecam.metrics.performance_metrics import pareto_front

# Assuming 'objectives' is a numpy array of objective vectors
pf = pareto_front(objectives)
```

### Evaluation Counting and Timing

Track the number of function evaluations and execution time:

```python
from moecam.metrics.performance_metrics import EvaluationCounter

counter = EvaluationCounter()
counter.start()

# Run your algorithm here
# counter.increment() should be called for each function evaluation

counter.stop()
print(f"Evaluations: {counter.get_count()}")
print(f"Time: {counter.get_time():.2f}s")
```

## C++ Integration

MOECAM supports integration with C++ algorithms via CFFI. This allows you to:

1. Wrap existing C++ optimization libraries
2. Implement performance-critical algorithms in C++
3. Use callback functions to evaluate Python-defined objective functions from C++

### Basic C++ Integration Example

```python
from moecam.core.cpp_bindings import ToyCppLib

# Example of calling C++ functions
cpp_lib = ToyCppLib()
result = cpp_lib.add(3.14, 2.86)
cpp_lib.print_message("Hello from Python!")
```

## Examples

The `examples/` directory contains several complete examples:

- `basic_usage.py`: Basic usage of MOECAM
- `algorithm_comparison.py`: Comparing different algorithms
- `custom_problem.py`: Defining and solving custom problems

## API Reference

### Problems Module (`moecam.problems`)

#### `test_functions.py`

- `ZDT1(n_dim=30)`: ZDT1 test function
- `ZDT2(n_dim=30)`: ZDT2 test function
- `ZDT3(n_dim=30)`: ZDT3 test function
- `DTLZ1(n_dim=7, n_obj=3)`: DTLZ1 test function
- `DTLZ2(n_dim=12, n_obj=3)`: DTLZ2 test function

### Algorithms Module (`moecam.algorithms`)

#### `moea_algorithms.py`

- `NSGAII(problem, pop_size=100, num_generations=250, crossover_rate=0.9, mutation_rate=0.01)`
- `MOEAD(problem, pop_size=100, num_generations=250)`

### Metrics Module (`moecam.metrics`)

#### `performance_metrics.py`

- `pareto_front(points)`: Extract Pareto front from objective vectors
- `hypervolume(pareto_front_points, reference_point)`: Calculate hypervolume
- `EvaluationCounter()`: Track evaluations and timing

### Core Module (`moecam.core`)

#### `cffi_interface.py`

- CFFI definitions and library loading

#### `cpp_bindings.py`

- Python wrappers for C++ functionality

## Troubleshooting

### Common Issues

1. **Import Errors**: Make sure MOECAM is properly installed and in your Python path
2. **C++ Library Issues**: Ensure the C++ shared library is compiled and accessible
3. **Memory Issues**: For large problems, consider reducing population size or number of generations

### Getting Help

- Check the examples in the `examples/` directory
- Review the API documentation
- Submit issues on GitHub: https://github.com/smartpour/moecam/issues

## Contributing

We welcome contributions! Please see `CONTRIBUTING.md` for guidelines on:

- Reporting bugs
- Suggesting enhancements
- Contributing code
- Writing documentation

## License

MOECAM is released under the MIT License. See `LICENSE` for details.
