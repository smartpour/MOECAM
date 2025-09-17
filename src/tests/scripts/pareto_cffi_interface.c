
#include "pareto_cffi_interface.h"
#include <stdlib.h>
#include <string.h>

// Simple Pareto dominance check
int dominates(double* point1, double* point2, int num_objectives) {
    int better_in_any = 0;
    for (int i = 0; i < num_objectives; i++) {
        if (point1[i] > point2[i]) {
            return 0; // point1 is worse in objective i
        }
        if (point1[i] < point2[i]) {
            better_in_any = 1;
        }
    }
    return better_in_any;
}

int extract_pareto_front(double* all_points, int num_points, int num_objectives,
                        double* pareto_points, int max_pareto_points) {
    int pareto_count = 0;

    for (int i = 0; i < num_points && pareto_count < max_pareto_points; i++) {
        int is_dominated = 0;
        double* current_point = &all_points[i * num_objectives];

        // Check if current point is dominated by any other point
        for (int j = 0; j < num_points; j++) {
            if (i != j) {
                double* other_point = &all_points[j * num_objectives];
                if (dominates(other_point, current_point, num_objectives)) {
                    is_dominated = 1;
                    break;
                }
            }
        }

        // If not dominated, add to Pareto front
        if (!is_dominated) {
            memcpy(&pareto_points[pareto_count * num_objectives], current_point,
                   num_objectives * sizeof(double));
            pareto_count++;
        }
    }

    return pareto_count;
}
