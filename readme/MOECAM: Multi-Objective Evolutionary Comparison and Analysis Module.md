# MOECAM: Multi-Objective Evolutionary Comparison and Analysis Module

MOECAM is a comprehensive Python package for multi-objective optimization that provides a unified interface to various multi-objective evolutionary algorithms (MOEAs) and test problems. The package is designed to facilitate research and comparison of different multi-objective optimization approaches.

## Features

- **Multiple Algorithms**: Implementations of popular MOEAs including NSGA-II and MOEA/D
- **Test Problems**: Comprehensive suite of test functions including ZDT, DTLZ, and WFG series
- **Performance Metrics**: Built-in metrics for evaluating algorithm performance (Pareto front, hypervolume, etc.)
- **C++ Integration**: CFFI-based interface for integrating high-performance C++ algorithms
- **Visualization**: Tools for plotting Pareto fronts and analyzing results
- **Scalable**: Support for problems with varying dimensions and objectives

## Installation

```bash
pip install moecam
```

## Quick Start

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

# Calculate performance metrics
reference_point = [1.1, 1.1]
hv = hypervolume(pareto_front, reference_point)
print(f"Hypervolume: {hv}")
```

## Documentation

For detailed documentation, examples, and API reference, please visit our [documentation site](https://moecam.readthedocs.io).

## Contributing

We welcome contributions! Please see our [contributing guidelines](CONTRIBUTING.md) for details.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

```bibtex
@software{moecam2024,
  title={MOECAM: Multi-Objective Evolutionary Comparison and Analysis Module},
  author={Kartik Sangwan},
  year={2024},
  url={https://github.com/smartpour/moecam}
}
```
