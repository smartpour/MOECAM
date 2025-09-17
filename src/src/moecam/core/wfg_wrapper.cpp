/**
 * WFG Hypervolume CFFI Wrapper
 * ============================
 *
 * C++ wrapper for WFG hypervolume calculation that can be compiled
 * as a shared library using CFFI during package installation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

extern "C" {

// Include WFG header first
#include "wfg.h"

// We need to declare the functions from wfg.c since we can't include it directly
extern double hv(FRONT ps);
extern int n;
extern POINT ref;

// Global state for WFG
static bool wfg_initialized = false;
static FRONT current_front;
static int current_max_points = 0;
static int current_objectives = 0;

/**
 * Initialize WFG hypervolume calculation
 * @param num_objectives Number of objectives
 * @param max_points Maximum number of points to handle
 * @return 0 on success, -1 on error
 */
int wfg_init(int num_objectives, int max_points) {
    if (num_objectives < 2 || max_points < 1) {
        return -1;
    }

    // Clean up if already initialized
    if (wfg_initialized) {
        wfg_cleanup();
    }

    // Set global variables for WFG
    n = num_objectives;
    current_objectives = num_objectives;
    current_max_points = max_points;

    // Allocate memory for front
    current_front.points = (POINT*)malloc(max_points * sizeof(POINT));
    if (!current_front.points) {
        return -1;
    }

    // Allocate objectives for each point
    for (int i = 0; i < max_points; i++) {
        current_front.points[i].objectives = (double*)malloc(num_objectives * sizeof(double));
        if (!current_front.points[i].objectives) {
            // Clean up on error
            for (int j = 0; j < i; j++) {
                free(current_front.points[j].objectives);
            }
            free(current_front.points);
            return -1;
        }
    }

    // Allocate reference point
    ref.objectives = (double*)malloc(num_objectives * sizeof(double));
    if (!ref.objectives) {
        for (int i = 0; i < max_points; i++) {
            free(current_front.points[i].objectives);
        }
        free(current_front.points);
        return -1;
    }

    wfg_initialized = true;
    return 0;
}

/**
 * Calculate hypervolume of a set of points
 * @param points Array of points (flattened: [x1_obj1, x1_obj2, ..., x2_obj1, x2_obj2, ...])
 * @param reference_point Reference point for hypervolume calculation
 * @param num_points Number of points
 * @param num_objectives Number of objectives
 * @return Hypervolume value, or -1.0 on error
 */
double wfg_calculate_hypervolume(double* points, double* reference_point,
                                int num_points, int num_objectives) {
    if (!wfg_initialized || !points || !reference_point) {
        return -1.0;
    }

    if (num_points > current_max_points || num_objectives != current_objectives) {
        return -1.0;
    }

    // Set reference point
    for (int i = 0; i < num_objectives; i++) {
        ref.objectives[i] = reference_point[i];
    }

    // Copy points to WFG format
    current_front.nPoints = num_points;
    for (int i = 0; i < num_points; i++) {
        for (int j = 0; j < num_objectives; j++) {
            current_front.points[i].objectives[j] = points[i * num_objectives + j];
        }
    }

    // Calculate hypervolume
    double volume = hv(current_front);
    return volume;
}

/**
 * Clean up WFG resources
 */
void wfg_cleanup() {
    if (!wfg_initialized) {
        return;
    }

    if (current_front.points) {
        for (int i = 0; i < current_max_points; i++) {
            if (current_front.points[i].objectives) {
                free(current_front.points[i].objectives);
            }
        }
        free(current_front.points);
        current_front.points = NULL;
    }

    if (ref.objectives) {
        free(ref.objectives);
        ref.objectives = NULL;
    }

    wfg_initialized = false;
    current_max_points = 0;
    current_objectives = 0;
}

} // extern "C"
