
#ifndef WFG_CFFI_INTERFACE_H
#define WFG_CFFI_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

// Initialize WFG hypervolume calculation
int wfg_init(int num_objectives, int max_points);

// Calculate hypervolume for a set of points
double wfg_calculate_hypervolume(double* points, double* reference_point, int num_points, int num_objectives);

// Cleanup WFG resources
void wfg_cleanup();

#ifdef __cplusplus
}
#endif

#endif
