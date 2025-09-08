# MOECAM Tests Directory Structure

Organized test suite for the Multi-Objective Evolutionary Comparison and Analysis Module.

## Directory Structure

```
src/tests/
├── scripts/                    # Demo and utility scripts
│   ├── quick_pareto_demo.py   # Comprehensive Pareto analysis with plots
│   ├── console_pareto_demo.py # Console-only Pareto demonstration
│   ├── objective_recorder.py  # Utility class for objective recording
│   └── README.md              # Script documentation
├── array_conversion/          # Array conversion testing
│   ├── test_array_conversion.py  # Comprehensive conversion tests
│   ├── test_edge_cases.py        # Edge cases and performance
│   └── README.md                 # Test documentation
├── test_results/              # Generated output files
│   ├── *.txt                  # Objective data files
│   ├── *.png                  # Plot files
│   └── README.md              # File format documentation
└── test_basic_functionality.py  # Original basic tests
```

## Quick Start

```bash
# Run Pareto analysis demo
cd src/tests/scripts/
python console_pareto_demo.py

# Test array conversion utilities
cd ../array_conversion/
python test_array_conversion.py

# Check generated results
ls ../test_results/
```

## What's Included

### ✅ Pareto Front Analysis

- Random sampling with objective recording
- Non-dominated sorting (Pareto front extraction)
- Multiple test problems (sphere, ZDT1, Kursawe)
- Console output and file saving
- Visualization capabilities

### ✅ Array Conversion Testing

- Robust Python↔C++ array handling
- Input type validation (lists, NumPy arrays, various dtypes)
- Memory safety and performance testing
- Edge case handling (large arrays, empty arrays, None values)

### ✅ C++/Python Integration

- CFFI callback interface validation
- Objective function evaluation pipeline
- Memory management testing
- Performance benchmarking

## Expected Outputs

All scripts generate timestamped files in `test_results/`:

- `*_objectives.txt` - Objective function evaluations
- `*_pareto_*.png` - Pareto front visualizations
- `recorded_objectives.txt` - Detailed evaluation logs

## Integration with Main Package

These tests validate the core functionality implemented in:

- `src/moecam/core/cffi_interface.py` - C++/Python interface
- `src/moecam/core/toy_cpp_lib.cpp` - C++ optimizer backend
