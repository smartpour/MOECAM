
#ifndef PARETO_CFFI_INTERFACE_H
#define PARETO_CFFI_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

// Extract Pareto front from a set of points
// Returns number of Pareto points found
int extract_pareto_front(double* all_points, int num_points, int num_objectives,
                        double* pareto_points, int max_pareto_points);

#ifdef __cplusplus
}
#endif

#endif
