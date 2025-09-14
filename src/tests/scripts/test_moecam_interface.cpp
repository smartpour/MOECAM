#include "moecam_python_interface.h"
#include <iostream>
#include <vector>
#include <cmath>

// Test objective function: simple bi-objective
void test_bi_objective(int n, const double *x, double *objectives, int num_objectives) {
    if (num_objectives >= 1) {
        // f1 = sum(x^2)
        double f1 = 0.0;
        for (int i = 0; i < n; i++) {
            f1 += x[i] * x[i];
        }
        objectives[0] = f1;
    }

    if (num_objectives >= 2) {
        // f2 = sum((x-1)^2)
        double f2 = 0.0;
        for (int i = 0; i < n; i++) {
            f2 += (x[i] - 1.0) * (x[i] - 1.0);
        }
        objectives[1] = f2;
    }
}

// Test objective function: Kursawe
void test_kursawe(int n, const double *x, double *objectives, int num_objectives) {
    if (num_objectives >= 1) {
        double f1 = 0.0;
        for (int i = 0; i < n - 1; i++) {
            f1 += -10.0 * std::exp(-0.2 * std::sqrt(x[i]*x[i] + x[i+1]*x[i+1]));
        }
        objectives[0] = f1;
    }

    if (num_objectives >= 2) {
        double f2 = 0.0;
        for (int i = 0; i < n; i++) {
            f2 += std::pow(std::abs(x[i]), 0.8) + 5.0 * std::sin(x[i]*x[i]*x[i]);
        }
        objectives[1] = f2;
    }
}

void test_algorithm(const char* name,
                   int (*algorithm_func)(int, double*, double*, python_callback_func, int, const double*, const double*, int),
                   python_callback_func callback,
                   int dim,
                   const std::vector<double>& bounds_lower,
                   const std::vector<double>& bounds_upper,
                   int max_iter) {

    std::cout << "\n=== Testing " << name << " ===" << std::endl;

    // Initial point (center of bounds)
    std::vector<double> x0(dim);
    for (int i = 0; i < dim; i++) {
        x0[i] = 0.5 * (bounds_lower[i] + bounds_upper[i]);
    }

    double result_value = 0.0;

    int retcode = algorithm_func(
        dim,
        x0.data(),
        &result_value,
        callback,
        2,  // num_objectives
        bounds_lower.data(),
        bounds_upper.data(),
        max_iter
    );

    std::cout << "Return code: " << retcode << std::endl;
    std::cout << "Final solution: [";
    for (int i = 0; i < dim; i++) {
        std::cout << x0[i];
        if (i < dim - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    std::cout << "Objective value: " << result_value << std::endl;

    // Verify by calling objective function directly
    std::vector<double> objectives(2);
    callback(dim, x0.data(), objectives.data(), 2);
    std::cout << "Verification - Objectives: [" << objectives[0] << ", " << objectives[1] << "]" << std::endl;
}

int main() {
    std::cout << "MOECAM C++ Interface Test" << std::endl;
    std::cout << "=========================" << std::endl;

    // Initialize MOECAM
    if (moecam_init() != 0) {
        std::cerr << "Failed to initialize MOECAM" << std::endl;
        return 1;
    }

    // Problem setup
    int dim = 3;
    std::vector<double> bounds_lower(dim, -5.0);
    std::vector<double> bounds_upper(dim, 5.0);

    // Test 1: ECAM with simple bi-objective
    std::cout << "\nTest 1: Simple bi-objective function" << std::endl;
    test_algorithm("ECAM", moecam_minimize_ecam, test_bi_objective,
                  dim, bounds_lower, bounds_upper, 100);

    test_algorithm("Random Start", moecam_minimize_random_start, test_bi_objective,
                  dim, bounds_lower, bounds_upper, 100);

    test_algorithm("DIRECT", moecam_direct_optimize, test_bi_objective,
                  dim, bounds_lower, bounds_upper, 100);

    // Test 2: Kursawe function
    std::cout << "\n\nTest 2: Kursawe function" << std::endl;
    test_algorithm("ECAM", moecam_minimize_ecam, test_kursawe,
                  dim, bounds_lower, bounds_upper, 100);

    test_algorithm("Random Start", moecam_minimize_random_start, test_kursawe,
                  dim, bounds_lower, bounds_upper, 100);

    // Cleanup
    moecam_cleanup();

    std::cout << "\n✅ C++ interface testing complete!" << std::endl;
    return 0;
}
