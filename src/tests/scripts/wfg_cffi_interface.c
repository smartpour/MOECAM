
#include "wfg_cffi_interface.h"
#include "/Users/OTomar/Downloads/Create Project Using Multi-Objective Optimization Resources/csources/wfg/WFG_1.15/wfg.h"
#include <stdlib.h>
#include <string.h>

// External declaration of WFG function
extern double hv(FRONT);

static int initialized = 0;
static int current_objectives = 0;

int wfg_init(int num_objectives, int max_points) {
    current_objectives = num_objectives;
    initialized = 1;
    return 0;
}

double wfg_calculate_hypervolume(double* points, double* reference_point, int num_points, int num_objectives) {
    if (!initialized || num_objectives != current_objectives) {
        return -1.0;
    }

    // Create FRONT structure
    FRONT front;
    front.nPoints = num_points;
    front.n = num_objectives;
    front.points = (POINT*)malloc(num_points * sizeof(POINT));

    if (!front.points) {
        return -1.0;
    }

    // Fill points from input array
    for (int i = 0; i < num_points; i++) {
        front.points[i].objectives = (OBJECTIVE*)malloc(num_objectives * sizeof(OBJECTIVE));
        if (!front.points[i].objectives) {
            // Cleanup allocated memory
            for (int j = 0; j < i; j++) {
                free(front.points[j].objectives);
            }
            free(front.points);
            return -1.0;
        }

        for (int j = 0; j < num_objectives; j++) {
            // WFG expects minimization, so transform points relative to reference
            front.points[i].objectives[j] = reference_point[j] - points[i * num_objectives + j];
        }
    }

    // Calculate hypervolume using WFG function
    double volume = hv(front);

    // Cleanup
    for (int i = 0; i < num_points; i++) {
        free(front.points[i].objectives);
    }
    free(front.points);

    return volume;
}

void wfg_cleanup() {
    initialized = 0;
}
