
import cffi
import os
import numpy as np

_this_dir = os.path.dirname(__file__)
_lib_path = os.path.join(_this_dir, 'toy_cpp_lib.so')

ffi = cffi.FFI()
ffi.cdef("""
// Legacy functions
double add(double a, double b);
void print_message(const char* message);

// New optimizer interface
int create_optimizer(int n_dim, int n_obj, double* lower_bounds, double* upper_bounds);
void destroy_optimizer(int handle);
void generate_random_solution(int handle, double* solution);
void crossover_solutions(int handle, const double* parent1, const double* parent2,
                       double* child1, double* child2);
void mutate_solution(int handle, double* individual, double mutation_rate);
int get_optimizer_dimensions(int handle);
int get_optimizer_objectives(int handle);
""")

lib = ffi.dlopen(_lib_path)

class CppOptimizer:
    """Python wrapper for C++ MOOptimizer class"""

    def __init__(self, n_dim, n_obj, lower_bounds, upper_bounds):
        self.n_dim = n_dim
        self.n_obj = n_obj

        # Convert bounds to C arrays
        self.lb_c = ffi.new("double[]", lower_bounds)
        self.ub_c = ffi.new("double[]", upper_bounds)

        # Create C++ optimizer instance
        self.handle = lib.create_optimizer(n_dim, n_obj, self.lb_c, self.ub_c)

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
        # Convert to C arrays
        p1_c = ffi.new("double[]", parent1.tolist())
        p2_c = ffi.new("double[]", parent2.tolist())
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
        # Convert to C array
        ind_c = ffi.new("double[]", individual.tolist())

        # Call C++ mutation
        lib.mutate_solution(self.handle, ind_c, mutation_rate)

        # Convert back to numpy array
        return np.array([ind_c[i] for i in range(self.n_dim)])

    def get_info(self):
        """Get optimizer information"""
        return {
            'dimensions': lib.get_optimizer_dimensions(self.handle),
            'objectives': lib.get_optimizer_objectives(self.handle)
        }

# Legacy interface functions
def cpp_add(a, b):
    """Legacy add function"""
    return lib.add(a, b)

def cpp_print(message):
    """Legacy print function"""
    lib.print_message(message.encode('utf-8'))


