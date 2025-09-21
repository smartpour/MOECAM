# Simplified MOECAM Installation/Compilation

## Overview

Following the suggestion to simplify the installation process based on the `examplelibrary` template, this simplified version of MOECAM demonstrates a much cleaner approach to CFFI compilation.

## Key Improvements

### 1. Simple Structure

The simplified structure follows the `examplelibrary` pattern:

```
simplified_moecam/
├── setup.py                 # Simple setup script
├── README.md                # Documentation
├── src/
│   ├── build_moecam.py      # CFFI build script (like build_liblip.py)
│   └── moecam.py           # Python wrapper module
└── test_simplified_moecam.py # Test script
```

### 2. Minimal setup.py

The `setup.py` file is now extremely simple, just like the `examplelibrary`:

- Lists package info
- Calls the install module with `cffi_modules=['./src/build_moecam.py:ffibuilder']`
- Nothing else complex

### 3. Simple CFFI Build Script

The `build_moecam.py` follows the `build_liblip.py` pattern:

#### C Function Definitions

```python
ffibuilder.cdef(r"""
    // Basic functions
    double add(double a, double b);
    void print_message(const char* message);

    // Pareto front extraction
    int extract_pareto_front(double* all_points, int num_points, int num_objectives,
                            double* pareto_points, int max_pareto_points, int strict_mode);

    // MOECAM optimization functions
    int minimize_ecam(int dim, double* x0, double *val, double lc, double* xl, double* xu, int maxiter, int iterlocal);
    // ... more functions
""")
```

#### Source File List

```python
moecam_src = [
    '../csources/wfg/WFG_1.15/wfg.c',
    '../csources/wfg/WFG_1.15/read.c',
    # Add more source files as needed
]
```

#### Simple Compilation

```python
ffibuilder.set_source("_moecam",
    # Inline C code with function implementations
    sources=moecam_src,
    include_dirs=include_dirs,
    extra_compile_args=['-O3'],
    extra_link_args=['-lm']
)
```

### 4. Clean Python Wrapper

The Python module (`moecam.py`) provides a clean interface:

- Imports the compiled CFFI library
- Provides high-level Python functions
- Handles numpy array conversions
- Graceful error handling when CFFI is not available

## Benefits of This Approach

### 1. Simplicity

- Much less complex wrapper code
- Easy to understand and maintain
- Follows standard CFFI patterns

### 2. Maintainability

- Clear separation of concerns
- Minimal technical wrapper code
- Easy to update in 6 months

### 3. Extensibility

- Easy to add new C/C++ source files to `moecam_src` list
- Simple to add new function definitions to `ffibuilder.cdef()`
- Offloads complexity to CFFI package itself

### 4. Error Reduction

- Less custom code = fewer opportunities for errors
- Standard CFFI patterns are well-tested
- Simple structure is easier to debug

## Current Implementation Status

### ✅ Working Features

- Basic utility functions (add, print_message)
- Pareto front extraction (pure C implementation)
- ZDT test functions (ZDT1, ZDT2, ZDT3)
- Placeholder optimization algorithms
- Placeholder hypervolume calculation
- Clean Python interface with numpy integration

### 🔄 Next Steps (when needed)

1. **Add real MOECAM algorithms**: Gradually include the actual C++ source files
2. **Add WFG hypervolume**: Include the real WFG implementation
3. **Add more test functions**: Include additional benchmark functions
4. **Add real optimization**: Replace placeholders with actual algorithms

## Installation

```bash
cd simplified_moecam
pip install -e .
```

## Usage

```python
import moecam
import numpy as np

# Test basic functionality
print(moecam.hello_moecam())

# Extract Pareto front
points = np.random.rand(50, 2)
pareto_front = moecam.extract_pareto_front(points)

# Evaluate test functions
x = np.array([0.5, 0.3, 0.2])
objectives = moecam.evaluate_zdt1(x)

# Use optimization algorithms (placeholders for now)
x0 = np.array([0.5, 0.5])
bounds = np.array([[0.0, 1.0], [0.0, 1.0]])
result = moecam.minimize_ecam(x0, bounds)
```

## Philosophy

This approach follows the principle that **the wrapper should be as simple as possible**. All the complexity is offloaded to the CFFI package itself, which:

1. Handles cross-platform compilation
2. Manages memory allocation
3. Provides robust C-Python interfacing
4. Is well-tested and maintained

The result is a much cleaner, more maintainable codebase that follows standard patterns and is easy to understand and modify.

## Comparison with Previous Approach

| Aspect              | Previous Complex Approach       | Simplified Approach           |
| ------------------- | ------------------------------- | ----------------------------- |
| setup.py            | Complex with custom build logic | Simple, standard CFFI pattern |
| Build script        | Multiple complex files          | Single, clear build script    |
| C/C++ handling      | Custom compilation logic        | CFFI handles everything       |
| Maintainability     | Difficult to modify             | Easy to understand and update |
| Error prone         | Many custom components          | Minimal custom code           |
| Following standards | Custom patterns                 | Standard CFFI patterns        |

This simplified approach demonstrates that complex functionality can be achieved with much simpler infrastructure, making the codebase more maintainable and less error-prone.
