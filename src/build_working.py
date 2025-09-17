#!/usr/bin/env python3
"""
MOECAM Working CFFI Build Script
===============================

Minimal working version that compiles successfully.
"""

import cffi
import os
from pathlib import Path

def build_working_extension():
    """Build a working MOECAM extension with basic functionality."""

    ffibuilder = cffi.FFI()

    # Define the C interface
    ffibuilder.cdef("""
        // Toy example functions
        double add(double a, double b);
        void print_message(const char* message);

        // Pareto Front Extraction
        int extract_pareto_front(double* all_points, int num_points, int num_objectives,
                                double* pareto_points, int max_pareto_points, int strict_mode);
        int extract_pareto_indices(double* all_points, int num_points, int num_objectives,
                                  int* pareto_indices, int max_pareto_points, int strict_mode);
    """)

    # Find source directory
    this_dir = Path(__file__).parent
    src_dir = this_dir / "src" / "moecam" / "core"

    # Include the implementation directly in the source
    source_code = """
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Toy example implementation
double add(double a, double b) {
    return a + b;
}

void print_message(const char* message) {
    printf("C says: %s\\n", message);
}

// Simple Pareto Front Extraction Implementation in C
int extract_pareto_front(double* all_points, int num_points, int num_objectives,
                        double* pareto_points, int max_pareto_points, int strict_mode) {

    if (!all_points || !pareto_points || num_points <= 0 || num_objectives <= 0) {
        return -1;
    }

    // Simple implementation: mark dominated points
    int* dominated = (int*)calloc(num_points, sizeof(int));
    if (!dominated) return -1;

    // Check each point against all others
    for (int i = 0; i < num_points; i++) {
        if (dominated[i]) continue;

        for (int j = 0; j < num_points; j++) {
            if (i == j || dominated[j]) continue;

            // Check if point j dominates point i
            int j_dominates_i = 1;
            int j_better_in_some = 0;

            for (int k = 0; k < num_objectives; k++) {
                double i_val = all_points[i * num_objectives + k];
                double j_val = all_points[j * num_objectives + k];

                if (strict_mode) {
                    if (j_val >= i_val) {
                        j_dominates_i = 0;
                        break;
                    }
                    if (j_val < i_val) {
                        j_better_in_some = 1;
                    }
                } else {
                    if (j_val > i_val) {
                        j_dominates_i = 0;
                        break;
                    }
                    if (j_val < i_val) {
                        j_better_in_some = 1;
                    }
                }
            }

            if (j_dominates_i && j_better_in_some) {
                dominated[i] = 1;
                break;
            }
        }
    }

    // Count non-dominated points
    int pareto_count = 0;
    for (int i = 0; i < num_points; i++) {
        if (!dominated[i]) {
            pareto_count++;
        }
    }

    // Check if we have enough space
    if (pareto_count > max_pareto_points) {
        free(dominated);
        return -1;
    }

    // Copy non-dominated points to output
    int output_idx = 0;
    for (int i = 0; i < num_points; i++) {
        if (!dominated[i]) {
            for (int j = 0; j < num_objectives; j++) {
                pareto_points[output_idx * num_objectives + j] = all_points[i * num_objectives + j];
            }
            output_idx++;
        }
    }

    free(dominated);
    return pareto_count;
}

int extract_pareto_indices(double* all_points, int num_points, int num_objectives,
                          int* pareto_indices, int max_pareto_points, int strict_mode) {

    if (!all_points || !pareto_indices || num_points <= 0 || num_objectives <= 0) {
        return -1;
    }

    // Simple implementation: mark dominated points
    int* dominated = (int*)calloc(num_points, sizeof(int));
    if (!dominated) return -1;

    // Check each point against all others (same logic as above)
    for (int i = 0; i < num_points; i++) {
        if (dominated[i]) continue;

        for (int j = 0; j < num_points; j++) {
            if (i == j || dominated[j]) continue;

            int j_dominates_i = 1;
            int j_better_in_some = 0;

            for (int k = 0; k < num_objectives; k++) {
                double i_val = all_points[i * num_objectives + k];
                double j_val = all_points[j * num_objectives + k];

                if (strict_mode) {
                    if (j_val >= i_val) {
                        j_dominates_i = 0;
                        break;
                    }
                    if (j_val < i_val) {
                        j_better_in_some = 1;
                    }
                } else {
                    if (j_val > i_val) {
                        j_dominates_i = 0;
                        break;
                    }
                    if (j_val < i_val) {
                        j_better_in_some = 1;
                    }
                }
            }

            if (j_dominates_i && j_better_in_some) {
                dominated[i] = 1;
                break;
            }
        }
    }

    // Count and copy indices
    int pareto_count = 0;
    for (int i = 0; i < num_points; i++) {
        if (!dominated[i]) {
            if (pareto_count >= max_pareto_points) {
                free(dominated);
                return -1;
            }
            pareto_indices[pareto_count] = i;
            pareto_count++;
        }
    }

    free(dominated);
    return pareto_count;
}
"""

    print("Building working MOECAM extension...")

    # Build the extension
    ffibuilder.set_source(
        "moecam._moecam_cffi",
        source_code,
        sources=[],  # No external sources for now
        extra_compile_args=["-O3"],
        libraries=["m"]
    )

    return ffibuilder

if __name__ == "__main__":
    ffibuilder = build_working_extension()
    ffibuilder.compile(verbose=True)
