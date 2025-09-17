/**
 * Pareto Front Extraction CFFI Wrapper
 * ====================================
 *
 * C++ wrapper for Pareto front extraction that can be compiled
 * as a shared library using CFFI during package installation.
 *
 * Based on paretofront.cpp but converted to library functions.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <random>
#include <chrono>
#include <string>
#include <list>
#include <set>
#include <vector>

using namespace std;

extern "C" {

/**
 * Extract Pareto front from a set of points
 * @param all_points Array of points (flattened: [x1_obj1, x1_obj2, ..., x2_obj1, x2_obj2, ...])
 * @param num_points Number of points
 * @param num_objectives Number of objectives per point
 * @param pareto_points Output array for Pareto front points (must be pre-allocated)
 * @param max_pareto_points Maximum number of Pareto points that can be stored
 * @param strict_mode 0 for weak Pareto dominance, 1 for strict Pareto dominance
 * @return Number of Pareto front points found, or -1 on error
 */
int extract_pareto_front(double* all_points, int num_points, int num_objectives,
                        double* pareto_points, int max_pareto_points, int strict_mode) {

    if (!all_points || !pareto_points || num_points <= 0 || num_objectives <= 0) {
        return -1;
    }

    // Use a set to store indices of Pareto optimal points
    set<int> pareto_indices;
    vector<int> erase_list;

    // Main Pareto dominance check loop
    for (int i = 0; i < num_points; i++) {
        bool dominated = false;
        erase_list.clear();

        // Check against all current Pareto points
        for (auto it = pareto_indices.begin(); it != pareto_indices.end(); it++) {
            int current_idx = *it;
            int non_dominated_count = 0;

            // Compare point i with current Pareto point
            if (strict_mode == 1) { // Strict Pareto dominance
                for (int j = 0; j < num_objectives; j++) {
                    if (all_points[i * num_objectives + j] >= all_points[current_idx * num_objectives + j]) {
                        non_dominated_count++;
                    }
                }
            } else { // Weak Pareto dominance
                for (int j = 0; j < num_objectives; j++) {
                    if (all_points[i * num_objectives + j] > all_points[current_idx * num_objectives + j]) {
                        non_dominated_count++;
                    }
                }
            }

            if (non_dominated_count == 0) {
                dominated = true; // Point i is dominated by current_idx
                break;
            } else if (non_dominated_count == num_objectives) {
                // Point i dominates current_idx, mark for removal
                erase_list.push_back(current_idx);
            }
        }

        if (!dominated) {
            // Remove dominated points
            for (int idx : erase_list) {
                pareto_indices.erase(idx);
            }
            // Add point i to Pareto set
            pareto_indices.insert(i);
        }
    }

    // Additional pass to ensure no dominated points remain in the Pareto set
    vector<int> final_erase;
    for (auto it1 = pareto_indices.begin(); it1 != pareto_indices.end(); it1++) {
        bool dominated = false;

        for (auto it2 = pareto_indices.begin(); it2 != pareto_indices.end(); it2++) {
            if (*it1 == *it2) continue;

            int non_dominated_count = 0;

            if (strict_mode == 1) { // Strict Pareto dominance
                for (int j = 0; j < num_objectives; j++) {
                    if (all_points[*it1 * num_objectives + j] >= all_points[*it2 * num_objectives + j]) {
                        non_dominated_count++;
                    }
                }
            } else { // Weak Pareto dominance
                for (int j = 0; j < num_objectives; j++) {
                    if (all_points[*it1 * num_objectives + j] > all_points[*it2 * num_objectives + j]) {
                        non_dominated_count++;
                    }
                }
            }

            if (non_dominated_count == 0) {
                dominated = true; // *it1 is dominated by *it2
                break;
            }
        }

        if (dominated) {
            final_erase.push_back(*it1);
        }
    }

    // Remove any remaining dominated points
    for (int idx : final_erase) {
        pareto_indices.erase(idx);
    }

    // Check if we have enough space for output
    int pareto_count = pareto_indices.size();
    if (pareto_count > max_pareto_points) {
        return -1; // Not enough space in output array
    }

    // Copy Pareto points to output array
    int output_idx = 0;
    for (auto it = pareto_indices.begin(); it != pareto_indices.end(); it++) {
        int point_idx = *it;
        for (int j = 0; j < num_objectives; j++) {
            pareto_points[output_idx * num_objectives + j] = all_points[point_idx * num_objectives + j];
        }
        output_idx++;
    }

    return pareto_count;
}

/**
 * Extract Pareto front indices from a set of points
 * @param all_points Array of points (flattened: [x1_obj1, x1_obj2, ..., x2_obj1, x2_obj2, ...])
 * @param num_points Number of points
 * @param num_objectives Number of objectives per point
 * @param pareto_indices Output array for Pareto front point indices (must be pre-allocated)
 * @param max_pareto_points Maximum number of Pareto indices that can be stored
 * @param strict_mode 0 for weak Pareto dominance, 1 for strict Pareto dominance
 * @return Number of Pareto front points found, or -1 on error
 */
int extract_pareto_indices(double* all_points, int num_points, int num_objectives,
                          int* pareto_indices, int max_pareto_points, int strict_mode) {

    if (!all_points || !pareto_indices || num_points <= 0 || num_objectives <= 0) {
        return -1;
    }

    // Use a set to store indices of Pareto optimal points
    set<int> pareto_set;
    vector<int> erase_list;

    // Main Pareto dominance check loop
    for (int i = 0; i < num_points; i++) {
        bool dominated = false;
        erase_list.clear();

        // Check against all current Pareto points
        for (auto it = pareto_set.begin(); it != pareto_set.end(); it++) {
            int current_idx = *it;
            int non_dominated_count = 0;

            // Compare point i with current Pareto point
            if (strict_mode == 1) { // Strict Pareto dominance
                for (int j = 0; j < num_objectives; j++) {
                    if (all_points[i * num_objectives + j] >= all_points[current_idx * num_objectives + j]) {
                        non_dominated_count++;
                    }
                }
            } else { // Weak Pareto dominance
                for (int j = 0; j < num_objectives; j++) {
                    if (all_points[i * num_objectives + j] > all_points[current_idx * num_objectives + j]) {
                        non_dominated_count++;
                    }
                }
            }

            if (non_dominated_count == 0) {
                dominated = true; // Point i is dominated by current_idx
                break;
            } else if (non_dominated_count == num_objectives) {
                // Point i dominates current_idx, mark for removal
                erase_list.push_back(current_idx);
            }
        }

        if (!dominated) {
            // Remove dominated points
            for (int idx : erase_list) {
                pareto_set.erase(idx);
            }
            // Add point i to Pareto set
            pareto_set.insert(i);
        }
    }

    // Check if we have enough space for output
    int pareto_count = pareto_set.size();
    if (pareto_count > max_pareto_points) {
        return -1; // Not enough space in output array
    }

    // Copy Pareto indices to output array
    int output_idx = 0;
    for (auto it = pareto_set.begin(); it != pareto_set.end(); it++) {
        pareto_indices[output_idx] = *it;
        output_idx++;
    }

    return pareto_count;
}

} // extern "C"
