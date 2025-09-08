# Array Conversion Tests

This directory contains tests for the robust array conversion utilities that ensure safe Python↔C++ interoperability.

## Files

- **`test_array_conversion.py`** - Comprehensive tests for array conversion functions
- **`test_edge_cases.py`** - Edge case validation and performance testing

## What's Tested

### Array Conversion Utilities

- `convert_py_float_to_cffi()` - Python arrays → contiguous float64 + CFFI pointers
- `convert_py_int_to_cffi()` - Python arrays → contiguous int32 + CFFI pointers

### Input Types Tested

- Python lists
- NumPy arrays (various dtypes: float32, float64, int32, int64)
- Non-contiguous arrays (e.g., `arr[::2]`)
- Empty arrays
- None values
- Large arrays (100+ dimensions)

### CppOptimizer Integration

- Bounds specification with various array types
- Evaluation with different input formats
- Crossover and mutation operations
- Memory safety (no modification of input arrays)

## Usage

```bash
# From src/tests/array_conversion/
python test_array_conversion.py  # Full array conversion test suite
python test_edge_cases.py        # Edge cases and performance tests
```

## Expected Results

All tests should pass, showing:

- ✅ Successful conversion for all input types
- ✅ Guaranteed contiguous memory layout
- ✅ Correct data types (float64, int32)
- ✅ No memory leaks or corruption
- ✅ Consistent performance (~5-6 μs per evaluation)
