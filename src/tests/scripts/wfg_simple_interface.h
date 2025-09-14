#ifndef WFG_SIMPLE_INTERFACE_H
#define WFG_SIMPLE_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

// Simple C interface that calls the WFG algorithm directly
// This avoids conflicts with global variables by using a separate compilation unit

// Initialize WFG for given dimensions and max points
int wfg_simple_init(int dimensions, int max_points);

// Calculate hypervolume
// pts: flattened array of points (nPoints * dimensions)
// nadir: reference point (dimensions elements)
// nPoints: number of points
// dimensions: number of objectives
double wfg_simple_calculate(const double* pts, const double* nadir,
                           int nPoints, int dimensions);

// Cleanup
void wfg_simple_cleanup();

#ifdef __cplusplus
}
#endif

#endif // WFG_SIMPLE_INTERFACE_H
