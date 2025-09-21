#!/usr/bin/env python3
"""
Example of how to extend the simplified MOECAM with additional functionality.

This shows how easy it is to add new functions following the template pattern.
"""

def add_new_function_example():
    """
    Example of how to add a new function to the simplified MOECAM.

    Steps:
    1. Add function signature to ffibuilder.cdef() in build_moecam.py
    2. Add implementation to the C source code in build_moecam.py
    3. Add Python wrapper function to moecam.py
    4. Test the new function
    """

    # Step 1: Add to ffibuilder.cdef()
    cdef_addition = '''
    // New function signature
    double my_new_function(double x, double y);
    '''

    # Step 2: Add to C source code
    c_implementation = '''
    double my_new_function(double x, double y) {
        return x * y + 1.0;
    }
    '''

    # Step 3: Add Python wrapper
    python_wrapper = '''
    def my_new_function(x: float, y: float) -> float:
        """My new function that multiplies x and y and adds 1."""
        if not CFFI_AVAILABLE:
            raise RuntimeError("MOECAM CFFI library not available")

        return float(lib.my_new_function(x, y))
    '''

    print("Example of how to extend simplified MOECAM:")
    print("1. Add to C function definitions:", cdef_addition)
    print("2. Add C implementation:", c_implementation)
    print("3. Add Python wrapper:", python_wrapper)
    print("4. Rebuild with: pip install -e .")

if __name__ == "__main__":
    add_new_function_example()
