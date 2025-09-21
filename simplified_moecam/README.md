# MOECAM - Multi-Objective Evolutionary Comparison and Analysis Module

A simplified implementation of MOECAM algorithms with CFFI compilation.

## Features

- Multi-objective optimization algorithms (ECAM, DFBM, Random Start)
- Pareto front extraction
- WFG hypervolume calculation
- Simple installation with CFFI

## Installation

```bash
pip install -e .
```

## Usage

```python
import moecam

# Use MOECAM algorithms
result = moecam.minimize_ecam(...)

# Extract Pareto front
pareto_points = moecam.extract_pareto_front(points)
```
