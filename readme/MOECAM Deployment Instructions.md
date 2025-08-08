# MOECAM Deployment Instructions

## Overview

This document provides instructions for deploying and using the MOECAM (Multi-Objective Evolutionary Comparison and Analysis Module) package.

## Package Contents

The MOECAM project includes:

- `src/moecam/` - Main package source code
- `examples/` - Usage examples and demonstrations
- `tests/` - Unit tests and validation scripts
- `docs/` - Documentation and user manual
- `README.md` - Project overview and quick start guide
- `setup.py` - Package installation script
- `MOECAM_architecture.md` - Detailed architecture documentation

## Installation Options

### Option 1: Development Installation

For development and testing:

```bash
# Extract the project archive
cd MOECAM_project/

# Install in development mode
pip install -e .

# Install development dependencies
pip install -e .[dev]
```

### Option 2: Package Installation

For production use:

```bash
# Build the package
python setup.py sdist bdist_wheel

# Install the built package
pip install dist/moecam-0.1.0-py3-none-any.whl
```

### Option 3: PyPI Installation (Future)

Once published to PyPI:

```bash
pip install moecam
```

## Dependencies

The package requires:

- Python >= 3.7
- NumPy >= 1.19.0
- Matplotlib >= 3.3.0
- CFFI >= 1.14.0

## Quick Start

After installation, you can run the basic example:

```bash
cd examples/
python basic_usage.py
```

This will:

1. Run NSGA-II on the ZDT1 problem
2. Calculate performance metrics
3. Generate a Pareto front visualization

## Running Tests

To validate the installation:

```bash
cd tests/
python test_basic_functionality.py
```

All tests should pass, indicating the package is working correctly.

## C++ Integration

The package includes a toy C++ library demonstrating CFFI integration:

- `src/moecam/core/toy_cpp_lib.cpp` - C++ source code
- `src/moecam/core/toy_cpp_lib.so` - Compiled shared library
- `src/moecam/core/cffi_interface.py` - CFFI interface definitions
- `src/moecam/core/cpp_bindings.py` - Python wrapper classes

To recompile the C++ library:

```bash
cd src/moecam/core/
g++ -shared -fPIC -o toy_cpp_lib.so toy_cpp_lib.cpp
```

## Usage Examples

### Basic Optimization

```python
from moecam.problems.test_functions import ZDT1
from moecam.algorithms.moea_algorithms import NSGAII

problem = ZDT1(n_dim=30)
algorithm = NSGAII(problem, pop_size=100, num_generations=250)
pareto_front = algorithm.optimize()
```

### Algorithm Comparison

```python
from moecam.algorithms.moea_algorithms import NSGAII, MOEAD
from moecam.metrics.performance_metrics import hypervolume

# Compare algorithms
algorithms = {
    'NSGA-II': NSGAII(problem),
    'MOEA/D': MOEAD(problem)
}

for name, alg in algorithms.items():
    pf = alg.optimize()
    hv = hypervolume(pf, [1.1, 1.1])
    print(f"{name}: HV = {hv:.4f}")
```

## Extending the Package

### Adding New Test Problems

Create a new class in `src/moecam/problems/test_functions.py`:

```python
class NewProblem:
    def __init__(self, n_dim=10):
        self.n_dim = n_dim
        self.n_obj = 2
        self.xl = np.zeros(n_dim)
        self.xu = np.ones(n_dim)

    def evaluate(self, x):
        # Define your objective functions here
        f1 = np.sum(x**2)
        f2 = np.sum((x - 1)**2)
        return np.array([f1, f2])
```

### Adding New Algorithms

Create a new class in `src/moecam/algorithms/moea_algorithms.py`:

```python
class NewAlgorithm:
    def __init__(self, problem, **kwargs):
        self.problem = problem
        # Initialize algorithm parameters

    def optimize(self):
        # Implement your algorithm here
        # Return Pareto front as numpy array
        pass
```

## Troubleshooting

### Common Issues

1. **Import Errors**: Ensure the package is properly installed and in your Python path
2. **C++ Library Issues**: Recompile the shared library if needed
3. **Memory Issues**: Reduce population size or generations for large problems

### Getting Help

- Check the documentation in `docs/user_manual.md`
- Review examples in the `examples/` directory
- Run tests to verify installation

## Performance Considerations

- For large problems, consider using smaller population sizes initially
- The current implementation is simplified for demonstration purposes
- For production use, consider implementing more sophisticated algorithms

## Future Enhancements

Potential improvements include:

1. More sophisticated algorithm implementations
2. Additional test problem suites
3. Advanced performance metrics
4. Parallel processing support
5. GUI interface
6. Integration with optimization libraries

## License

This project is released under the MIT License. See the LICENSE file for details.

## Support

For questions, issues, or contributions, please contact the development team or submit issues through the project repository.
