# Test Scripts Directory

This directory contains demonstration and utility scripts for the MOECAM project.

## Files

### Pareto Analysis Scripts

- **`quick_pareto_demo.py`** - Comprehensive Pareto front demonstration with multiple test problems
- **`console_pareto_demo.py`** - Simple console-only Pareto analysis with detailed output
- **`objective_recorder.py`** - Utility class for recording and analyzing objective evaluations

## Usage

Run scripts from this directory:

```bash
# From src/tests/scripts/
python quick_pareto_demo.py      # Visual analysis with plots
python console_pareto_demo.py    # Console-only analysis
python objective_recorder.py     # Demo of recorder utility
```

All outputs (plots, data files) are saved to `../test_results/`

## Features

- ✅ Random sampling via C++/Python callback interface
- ✅ Pareto front extraction using non-dominated sorting
- ✅ Multiple test problems (sphere, ZDT1, Kursawe)
- ✅ Automatic file saving with timestamps
- ✅ Console output with statistics and sample points
- ✅ Visualization (2D/3D scatter plots, parallel coordinates)

## Test Problems

1. **Sphere objectives**: `f1 = ||x||², f2 = ||x-1||²` (convex Pareto front)
2. **ZDT1**: Classic test problem with known analytical Pareto front
3. **Kursawe**: Multi-modal problem with disconnected Pareto front
