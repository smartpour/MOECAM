
import cffi
import os
import numpy as np

_this_dir = os.path.dirname(__file__)
_lib_path = os.path.join(_this_dir, 'toy_cpp_lib.so')

ffi = cffi.FFI()

# Array conversion utilities for safe C++ interop
def convert_py_float_to_cffi(x):
    """Convert Python array to contiguous float64 array with CFFI pointer.

    Args:
        x: Input array-like (list, numpy array, etc.) or None

    Returns:
        tuple: (numpy_array, cffi_pointer)
               numpy_array is guaranteed contiguous float64
               cffi_pointer can be passed to C++ as double*
    """
    if x is not None:
        if isinstance(x, np.ndarray) and x.flags.c_contiguous and x.dtype == np.float64:
            px = x
        else:
            px = np.ascontiguousarray(x, dtype='float64')
        pxcffi = ffi.cast("double *", px.ctypes.data)
    else:
        # Handle None case
        px = np.array([0.0], dtype='float64')
        pxcffi = ffi.cast("double *", 0)
    return px, pxcffi

def convert_py_int_to_cffi(x):
    """Convert Python array to contiguous int32 array with CFFI pointer.

    Args:
        x: Input array-like (list, numpy array, etc.) or None

    Returns:
        tuple: (numpy_array, cffi_pointer)
               numpy_array is guaranteed contiguous int32
               cffi_pointer can be passed to C++ as int*
    """
    if x is not None:
        if isinstance(x, np.ndarray) and x.flags.c_contiguous and x.dtype == np.int32:
            px = x
        else:
            px = np.ascontiguousarray(x, dtype='int32')
        pxcffi = ffi.cast("int *", px.ctypes.data)
    else:
        # Handle None case
        px = np.array([0], dtype='int32')
        pxcffi = ffi.cast("int *", 0)
    return px, pxcffi
ffi.cdef("""
// Legacy functions
double add(double a, double b);
void print_message(const char* message);

// Callback type
typedef void (*objective_callback_t)(const double* x, int n_dim, double* f, int n_obj);

// New optimizer interface
int create_optimizer(int n_dim, int n_obj, double* lower_bounds, double* upper_bounds, objective_callback_t callback);
void destroy_optimizer(int handle);
void generate_random_solution(int handle, double* solution);
void crossover_solutions(int handle, const double* parent1, const double* parent2,
                       double* child1, double* child2);
void mutate_solution(int handle, double* individual, double mutation_rate);
int get_optimizer_dimensions(int handle);
int get_optimizer_objectives(int handle);
void set_objective_callback(int handle, objective_callback_t cb);
void evaluate_solution(int handle, const double* solution, double* objectives);
void evaluate_random_solution_with_objectives(int handle, double* solution, double* objectives);
""")

lib = ffi.dlopen(_lib_path)

class CppOptimizer:
    """Python wrapper for C++ MOOptimizer class"""

    def __init__(self, n_dim, n_obj, lower_bounds, upper_bounds, objective_callback):
        if objective_callback is None:
            raise ValueError("objective_callback is required and cannot be None")

        self.n_dim = n_dim
        self.n_obj = n_obj

        # Convert bounds using safe conversion utilities
        self.lower_bounds, self.lb_c = convert_py_float_to_cffi(lower_bounds)
        self.upper_bounds, self.ub_c = convert_py_float_to_cffi(upper_bounds)

        # Validate bounds dimensions
        if self.lower_bounds.shape[0] != n_dim or self.upper_bounds.shape[0] != n_dim:
            raise ValueError("Bounds must match n_dim")

        # Set up callback first
        self._py_objective = objective_callback
        self._setup_callback()

        # Create C++ optimizer instance with callback
        self.handle = lib.create_optimizer(n_dim, n_obj, self.lb_c, self.ub_c, self._c_callback)
        if self.handle == -1:
            raise RuntimeError("Failed to create optimizer - invalid callback")

    def __del__(self):
        """Cleanup C++ resources"""
        if hasattr(self, 'handle'):
            lib.destroy_optimizer(self.handle)

    def generate_random_solution(self):
        """Generate a random solution within bounds"""
        solution_c = ffi.new("double[]", self.n_dim)
        lib.generate_random_solution(self.handle, solution_c)
        return np.array([solution_c[i] for i in range(self.n_dim)])

    def crossover(self, parent1, parent2):
        """Perform crossover between two parents"""
        # Convert to C arrays using safe conversion
        p1_np, p1_c = convert_py_float_to_cffi(parent1)
        p2_np, p2_c = convert_py_float_to_cffi(parent2)

        # Validate dimensions
        if p1_np.shape[0] != self.n_dim or p2_np.shape[0] != self.n_dim:
            raise ValueError("Parent dimensions must match optimizer n_dim")

        child1_c = ffi.new("double[]", self.n_dim)
        child2_c = ffi.new("double[]", self.n_dim)

        # Call C++ crossover
        lib.crossover_solutions(self.handle, p1_c, p2_c, child1_c, child2_c)

        # Convert back to numpy arrays
        child1 = np.array([child1_c[i] for i in range(self.n_dim)])
        child2 = np.array([child2_c[i] for i in range(self.n_dim)])

        return child1, child2

    def mutate(self, individual, mutation_rate=0.01):
        """Mutate an individual solution"""
        # Convert to C array using safe conversion
        ind_np, ind_c = convert_py_float_to_cffi(individual)

        # Validate dimensions
        if ind_np.shape[0] != self.n_dim:
            raise ValueError("Individual dimensions must match optimizer n_dim")

        # Call C++ mutation (modifies array in place)
        lib.mutate_solution(self.handle, ind_c, mutation_rate)

        # Convert back to numpy array
        return np.array([ind_c[i] for i in range(self.n_dim)])

    def get_info(self):
        """Get optimizer information"""
        return {
            'dimensions': lib.get_optimizer_dimensions(self.handle),
            'objectives': lib.get_optimizer_objectives(self.handle)
        }

    # --- Objective evaluation interface ---
    def _setup_callback(self):
        """Internal method to set up the C callback from Python function."""
        if not callable(self._py_objective):
            raise TypeError("objective_callback must be callable")

        @ffi.callback("void(const double*, int, double*, int)")
        def _c_obj(x_ptr, n_dim, f_ptr, n_obj):  # noqa: N803
            # Copy decision vector
            x = np.frombuffer(ffi.buffer(x_ptr, n_dim * 8), dtype=np.float64).copy()
            try:
                f = self._py_objective(x)
                f_arr = np.asarray(f, dtype=np.float64)
                if f_arr.shape[0] != n_obj:
                    raise ValueError("Objective callback returned wrong length")
            except Exception as e:  # Fallback: zeros
                f_arr = np.zeros(n_obj, dtype=np.float64)
                print(f"Objective callback error: {e}")
            # Write back to f_ptr
            for i in range(n_obj):
                f_ptr[i] = float(f_arr[i])

        self._c_callback = _c_obj  # Keep reference to prevent GC

    def set_objective_callback(self, py_func):
        """Update the Python objective callback.

        py_func: callable(x: np.ndarray) -> np.ndarray (length n_obj)
        """
        if not callable(py_func):
            raise TypeError("objective_callback must be callable")

        self._py_objective = py_func
        self._setup_callback()
        lib.set_objective_callback(self.handle, self._c_callback)

    def evaluate(self, x):
        """Evaluate solution using the registered callback"""
        x_np, x_c = convert_py_float_to_cffi(x)
        if x_np.shape[0] != self.n_dim:
            raise ValueError("x length mismatch")
        f_c = ffi.new("double[]", self.n_obj)
        lib.evaluate_solution(self.handle, x_c, f_c)
        return np.array([f_c[i] for i in range(self.n_obj)])

    def sample_random(self, n_samples):
        """Generate n random solutions and evaluate objectives via callback.
        Returns (X, F) where X shape (n_samples, n_dim), F shape (n_samples, n_obj)
        """
        X = np.zeros((n_samples, self.n_dim))
        F = np.zeros((n_samples, self.n_obj))
        for i in range(n_samples):
            x_c = ffi.new("double[]", self.n_dim)
            f_c = ffi.new("double[]", self.n_obj)
            lib.evaluate_random_solution_with_objectives(self.handle, x_c, f_c)
            for d in range(self.n_dim):
                X[i, d] = x_c[d]
            for j in range(self.n_obj):
                F[i, j] = f_c[j]
        return X, F

    def plot_objectives(self, n_samples=100, ax=None):
        import matplotlib.pyplot as plt
        if self.n_obj < 2:
            raise ValueError("Need at least 2 objectives for plotting")
        X, F = self.sample_random(n_samples)
        if ax is None:
            fig = plt.figure()
            if self.n_obj == 2:
                ax = fig.add_subplot(111)
            elif self.n_obj == 3:
                from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
                ax = fig.add_subplot(111, projection='3d')
        if self.n_obj == 2:
            ax.scatter(F[:,0], F[:,1], s=20, alpha=0.7)
            ax.set_xlabel('f1')
            ax.set_ylabel('f2')
        elif self.n_obj == 3:
            ax.scatter(F[:,0], F[:,1], F[:,2], s=20, alpha=0.7)
            ax.set_xlabel('f1')
            ax.set_ylabel('f2')
            ax.set_zlabel('f3')
        else:
            # Parallel coordinates fallback
            from pandas import DataFrame
            from pandas.plotting import parallel_coordinates
            df = DataFrame(F, columns=[f'f{i+1}' for i in range(self.n_obj)])
            df['class'] = 'pts'
            parallel_coordinates(df, 'class', color=['#1f77b4'])
        ax.set_title(f'Objective Samples (n={n_samples})')
        return ax

# Legacy interface functions
def cpp_add(a, b):
    """Legacy add function"""
    return lib.add(a, b)

def cpp_print(message):
    """Legacy print function"""
    lib.print_message(message.encode('utf-8'))


