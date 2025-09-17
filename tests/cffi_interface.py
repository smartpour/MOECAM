
import cffi
import os

_this_dir = os.path.dirname(__file__)
_lib_path = os.path.join(_this_dir, 'toy_cpp_lib.so')

ffi = cffi.FFI()
ffi.cdef("""
double add(double a, double b);
void print_message(const char* message);
""")

lib = ffi.dlopen(_lib_path)


