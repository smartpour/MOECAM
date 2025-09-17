/*
 * Simple MOECAM Test Wrapper
 * =========================
 *
 * Minimal C++ wrapper to test callback mechanism
 */

#include <vector>
#include <iostream>

extern "C" {

// Callback type
typedef void (*objective_function_t)(double* x, int n_dim, double* f, int n_obj);

// Simple test function that just calls the callback once
int test_callback(objective_function_t callback, double* test_x, int n_dim, double* result_f, int n_obj) {
    try {
        std::cout << "C++: About to call Python callback" << std::endl;

        // Call the Python callback
        callback(test_x, n_dim, result_f, n_obj);

        std::cout << "C++: Callback completed successfully" << std::endl;
        return 0;

    } catch (...) {
        std::cout << "C++: Exception in callback" << std::endl;
        return -1;
    }
}

} // extern "C"
