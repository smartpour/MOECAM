#include "wfg_simple_interface.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Include WFG header
extern "C" {
#include "../../../csources/wfg/WFG_1.15/wfg.h"
}

// Forward declarations of WFG functions
extern "C" {
    double hv(FRONT ps);
}

// Global state for this interface
static int g_initialized = 0;
static int g_max_points = 0;
static int g_max_dimensions = 0;

int wfg_simple_init(int dimensions, int max_points) {
    g_max_dimensions = dimensions;
    g_max_points = max_points;
    g_initialized = 1;
    return 1;
}

double wfg_simple_calculate(const double* pts, const double* nadir, 
                           int nPoints, int dimensions) {
    if (!g_initialized) {
        return 0.0;
    }
    
    if (nPoints <= 0 || dimensions <= 0) {
        return 0.0;
    }
    
    // Create FRONT structure
    FRONT front;
    front.n = dimensions;
    front.nPoints = nPoints;
    front.points = (POINT*)malloc(sizeof(POINT) * nPoints);
    
    // Fill the front with transformed points
    for (int j = 0; j < nPoints; j++) {
        front.points[j].objectives = (OBJECTIVE*)malloc(sizeof(OBJECTIVE) * dimensions);
        for (int k = 0; k < dimensions; k++) {
            // Transform: fabs(objective - reference_point)
            double objective_value = pts[j * dimensions + k];
            double reference_value = nadir[k];
            front.points[j].objectives[k] = fabs(objective_value - reference_value);
        }
    }
    
    // Calculate hypervolume using WFG algorithm
    double hypervolume = hv(front);
    
    // Cleanup
    for (int j = 0; j < nPoints; j++) {
        free(front.points[j].objectives);
    }
    free(front.points);
    
    return hypervolume;
}

void wfg_simple_cleanup() {
    g_initialized = 0;
    g_max_points = 0;
    g_max_dimensions = 0;
}
