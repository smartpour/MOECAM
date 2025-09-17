"""
CFFI Builder for MOECAM
=======================
"""
from cffi import FFI

ffibuilder = FFI()

# Define C functions for CFFI
ffibuilder.cdef("""
    // Basic functions
    void hello_world(void);
    int add_numbers(int a, int b);

    // Pareto front extraction
    int extract_pareto_front(double* points, int n_points, int n_objectives,
                           double* pareto_points, int* pareto_indices);

    // Simple optimization functions
    int simple_random_optimization(double* bounds, int n_vars, int n_objectives,
                                 int max_evaluations, double* solutions, double* objectives);
""")

# C source code
ffibuilder.set_source("moecam._moecam_cffi",
    r"""
    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include <string.h>

    // Basic functions
    void hello_world(void) {
        printf("Hello from MOECAM CFFI!\n");
    }

    int add_numbers(int a, int b) {
        return a + b;
    }

    // Pareto dominance check
    int dominates(double* p1, double* p2, int n_obj) {
        int at_least_one_better = 0;
        for (int i = 0; i < n_obj; i++) {
            if (p1[i] > p2[i]) return 0;  // p1 is worse in objective i
            if (p1[i] < p2[i]) at_least_one_better = 1;
        }
        return at_least_one_better;
    }

    // Extract Pareto front
    int extract_pareto_front(double* points, int n_points, int n_objectives,
                           double* pareto_points, int* pareto_indices) {
        int* is_pareto = (int*)calloc(n_points, sizeof(int));
        int pareto_count = 0;

        // Check each point for Pareto optimality
        for (int i = 0; i < n_points; i++) {
            is_pareto[i] = 1;  // Assume Pareto optimal

            for (int j = 0; j < n_points; j++) {
                if (i != j) {
                    if (dominates(&points[j * n_objectives], &points[i * n_objectives], n_objectives)) {
                        is_pareto[i] = 0;  // Point i is dominated
                        break;
                    }
                }
            }

            if (is_pareto[i]) {
                // Copy Pareto point
                for (int k = 0; k < n_objectives; k++) {
                    pareto_points[pareto_count * n_objectives + k] = points[i * n_objectives + k];
                }
                pareto_indices[pareto_count] = i;
                pareto_count++;
            }
        }

        free(is_pareto);
        return pareto_count;
    }

    // Simple random optimization (placeholder)
    int simple_random_optimization(double* bounds, int n_vars, int n_objectives,
                                 int max_evaluations, double* solutions, double* objectives) {
        // This is a placeholder - in real implementation would do proper optimization
        return max_evaluations;
    }
    """)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
