from cffi import FFI
import os

ffibuilder = FFI()
PATH = os.path.dirname(__file__)

ffibuilder.cdef(r"""
    // Basic utility functions
    double add(double a, double b);
    void print_message(const char* message);

    // Pareto Front Extraction Functions
    int extract_pareto_front(double* all_points, int num_points, int num_objectives,
                            double* pareto_points, int max_pareto_points, int strict_mode);
    int extract_pareto_indices(double* all_points, int num_points, int num_objectives,
                              int* pareto_indices, int max_pareto_points, int strict_mode);

    // MOECAM Optimization Functions (simplified C interface)
    int minimize_dfbm(int dim, double* x0, double *val, double* xl, double* xu, int maxiter);
    int minimize_ecam(int dim, double* x0, double *val, double lc, double* xl, double* xu, int maxiter, int iterlocal);
    int minimize_random_start(int dim, double* x0, double *val, double* xl, double* xu, int maxiter);

    // WFG Hypervolume Calculation
    double wfg_hypervolume(double* points, double* reference, int num_points, int num_objectives);
    void wfg_prepare(int max_objectives, int max_points);
    void wfg_cleanup(void);

    // Test function evaluation (ZDT1, ZDT2, etc.)
    void evaluate_zdt1(double* x, int n_vars, double* objectives);
    void evaluate_zdt2(double* x, int n_vars, double* objectives);
    void evaluate_zdt3(double* x, int n_vars, double* objectives);

""", override=True)

# Start with minimal source files (no complex C++ dependencies for now)
moecam_src = [
    # Start with just the basic WFG C files (they are simpler)
    '../csources/wfg/WFG_1.15/wfg.c',
    '../csources/wfg/WFG_1.15/read.c',
]

# Include directories - minimal set
include_dirs = [
    PATH,
    os.path.join(PATH, '../csources/wfg/WFG_1.15'),
]

ffibuilder.set_source("_moecam",r"""
    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include <string.h>
    #include <stdbool.h>

    // Simple C implementations (no C++ syntax)

    // Basic utility functions
    double add(double a, double b) {
        return a + b;
    }

    void print_message(const char* message) {
        printf("MOECAM says: %s\n", message);
    }

    // Simple Pareto front extraction (basic implementation)
    int extract_pareto_front(double* all_points, int num_points, int num_objectives,
                            double* pareto_points, int max_pareto_points, int strict_mode) {
        if (!all_points || !pareto_points || num_points <= 0 || num_objectives <= 0) {
            return 0;
        }

        int pareto_count = 0;

        for (int i = 0; i < num_points && pareto_count < max_pareto_points; i++) {
            bool is_dominated = false;

            // Check if point i is dominated by any other point
            for (int j = 0; j < num_points && !is_dominated; j++) {
                if (i == j) continue;

                bool dominates = true;
                bool strictly_better = false;

                for (int k = 0; k < num_objectives; k++) {
                    double val_i = all_points[i * num_objectives + k];
                    double val_j = all_points[j * num_objectives + k];

                    if (val_j > val_i) {  // Assuming minimization
                        dominates = false;
                        break;
                    }
                    if (val_j < val_i) {
                        strictly_better = true;
                    }
                }

                if (dominates && (strictly_better || !strict_mode)) {
                    is_dominated = true;
                }
            }

            if (!is_dominated) {
                // Copy point to pareto front
                for (int k = 0; k < num_objectives; k++) {
                    pareto_points[pareto_count * num_objectives + k] =
                        all_points[i * num_objectives + k];
                }
                pareto_count++;
            }
        }

        return pareto_count;
    }

    int extract_pareto_indices(double* all_points, int num_points, int num_objectives,
                              int* pareto_indices, int max_pareto_points, int strict_mode) {
        if (!all_points || !pareto_indices || num_points <= 0 || num_objectives <= 0) {
            return 0;
        }

        int pareto_count = 0;

        for (int i = 0; i < num_points && pareto_count < max_pareto_points; i++) {
            bool is_dominated = false;

            // Check if point i is dominated by any other point
            for (int j = 0; j < num_points && !is_dominated; j++) {
                if (i == j) continue;

                bool dominates = true;
                bool strictly_better = false;

                for (int k = 0; k < num_objectives; k++) {
                    double val_i = all_points[i * num_objectives + k];
                    double val_j = all_points[j * num_objectives + k];

                    if (val_j > val_i) {  // Assuming minimization
                        dominates = false;
                        break;
                    }
                    if (val_j < val_i) {
                        strictly_better = true;
                    }
                }

                if (dominates && (strictly_better || !strict_mode)) {
                    is_dominated = true;
                }
            }

            if (!is_dominated) {
                pareto_indices[pareto_count] = i;
                pareto_count++;
            }
        }

        return pareto_count;
    }

    // Simplified optimization function wrappers (placeholders for now)
    int minimize_dfbm(int dim, double* x0, double *val, double* xl, double* xu, int maxiter) {
        // Placeholder - simplified random optimization for now
        *val = 0.0; // Dummy value
        return 0; // Success
    }

    int minimize_ecam(int dim, double* x0, double *val, double lc, double* xl, double* xu, int maxiter, int iterlocal) {
        // Placeholder - simplified random optimization for now
        *val = 0.0; // Dummy value
        return 0; // Success
    }

    int minimize_random_start(int dim, double* x0, double *val, double* xl, double* xu, int maxiter) {
        // Placeholder - simplified random optimization for now
        *val = 0.0; // Dummy value
        return 0; // Success
    }

    // WFG hypervolume calculation (placeholder for now)
    double wfg_hypervolume(double* points, double* reference, int num_points, int num_objectives) {
        // Placeholder - simple calculation for now
        return 1.0; // Dummy value
    }

    void wfg_prepare(int max_objectives, int max_points) {
        // Placeholder for WFG initialization
    }

    void wfg_cleanup(void) {
        // Placeholder for WFG cleanup
    }

    // Test function implementations
    void evaluate_zdt1(double* x, int n_vars, double* objectives) {
        if (n_vars < 1) return;

        double f1 = x[0];

        double g = 0.0;
        for (int i = 1; i < n_vars; i++) {
            g += x[i];
        }
        g = 1.0 + 9.0 * g / (n_vars - 1);

        double f2 = g * (1.0 - sqrt(f1 / g));

        objectives[0] = f1;
        objectives[1] = f2;
    }

    void evaluate_zdt2(double* x, int n_vars, double* objectives) {
        if (n_vars < 1) return;

        double f1 = x[0];

        double g = 0.0;
        for (int i = 1; i < n_vars; i++) {
            g += x[i];
        }
        g = 1.0 + 9.0 * g / (n_vars - 1);

        double f2 = g * (1.0 - (f1 / g) * (f1 / g));

        objectives[0] = f1;
        objectives[1] = f2;
    }

    void evaluate_zdt3(double* x, int n_vars, double* objectives) {
        if (n_vars < 1) return;

        double f1 = x[0];

        double g = 0.0;
        for (int i = 1; i < n_vars; i++) {
            g += x[i];
        }
        g = 1.0 + 9.0 * g / (n_vars - 1);

        double h = 1.0 - sqrt(f1 / g) - (f1 / g) * sin(10.0 * M_PI * f1);
        double f2 = g * h;

        objectives[0] = f1;
        objectives[1] = f2;
    }
""",
    sources=moecam_src,
    include_dirs=include_dirs,
    extra_compile_args=['-O3'],
    extra_link_args=['-lm'],
    )

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
