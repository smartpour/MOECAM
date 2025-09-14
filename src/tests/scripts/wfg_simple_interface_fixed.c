#include "wfg_simple_interface.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

// External WFG variables and function
extern int n;
extern POINT ref;
extern FRONT *fs;
extern int fr;
extern int maxm;
extern int maxn;
extern int safe;
extern int** nextp;
extern int** prevp;
extern int* firstp;
extern int* lastp;
extern int* psize;
extern double hv(FRONT front);

static int wfg_initialized = 0;

int wfg_simple_init(int max_points, int max_objectives) {
    if (wfg_initialized) {
        wfg_simple_cleanup();
    }

    maxm = max_points;
    maxn = max_objectives;

    // Allocate memory for algorithm workspace
    int maxdepth = maxn - 2;
    if (maxdepth <= 0) maxdepth = 1;

    fs = malloc(sizeof(FRONT) * maxdepth);
    if (!fs) return 0;

    for (int i = 0; i < maxdepth; i++) {
        fs[i].points = malloc(sizeof(POINT) * maxm);
        if (!fs[i].points) return 0;
        for (int j = 0; j < maxm; j++) {
            int objectives_count = maxn - i - 1;
            if (objectives_count <= 0) objectives_count = 1;
            fs[i].points[j].objectives = malloc(sizeof(OBJECTIVE) * objectives_count);
            if (!fs[i].points[j].objectives) return 0;
        }
    }

    nextp = malloc(sizeof(int*) * maxdepth);
    if (!nextp) return 0;
    for (int i = 0; i < maxdepth; i++) {
        nextp[i] = malloc(sizeof(int) * maxm);
        if (!nextp[i]) return 0;
    }

    prevp = malloc(sizeof(int*) * maxdepth);
    if (!prevp) return 0;
    for (int i = 0; i < maxdepth; i++) {
        prevp[i] = malloc(sizeof(int) * maxm);
        if (!prevp[i]) return 0;
    }

    firstp = malloc(sizeof(int) * maxdepth);
    if (!firstp) return 0;

    lastp = malloc(sizeof(int) * maxdepth);
    if (!lastp) return 0;

    psize = malloc(sizeof(int) * maxdepth);
    if (!psize) return 0;

    // Initialize reference point
    ref.objectives = malloc(sizeof(OBJECTIVE) * maxn);
    if (!ref.objectives) return 0;

    wfg_initialized = 1;
    return 1;
}

double wfg_simple_calculate(double** points_data, int num_points, int num_objectives,
                           double* reference_point) {
    if (!wfg_initialized) {
        if (!wfg_simple_init(num_points + 100, num_objectives)) {
            return -1.0;
        }
    }

    // Create FRONT structure
    FRONT front;
    front.nPoints = num_points;
    front.n = num_objectives;
    front.points = malloc(sizeof(POINT) * num_points);
    if (!front.points) return -1.0;

    // Copy points and set up objectives
    for (int i = 0; i < num_points; i++) {
        front.points[i].objectives = malloc(sizeof(OBJECTIVE) * num_objectives);
        if (!front.points[i].objectives) {
            // Cleanup on failure
            for (int j = 0; j < i; j++) {
                free(front.points[j].objectives);
            }
            free(front.points);
            return -1.0;
        }

        for (int j = 0; j < num_objectives; j++) {
            front.points[i].objectives[j] = points_data[i][j];
        }
    }

    // Set reference point
    for (int i = 0; i < num_objectives; i++) {
        ref.objectives[i] = reference_point[i];
    }

    // Transform objectives relative to reference point (as done in original WFG)
    for (int i = 0; i < front.nPoints; i++) {
        for (int j = 0; j < front.n; j++) {
            front.points[i].objectives[j] = fabs(front.points[i].objectives[j] - ref.objectives[j]);
        }
    }

    // Set global variables for algorithm
    n = num_objectives;
    fr = 0;
    safe = 0;

    // Calculate hypervolume
    double result = hv(front);

    // Cleanup
    for (int i = 0; i < num_points; i++) {
        free(front.points[i].objectives);
    }
    free(front.points);

    return result;
}

void wfg_simple_cleanup() {
    if (!wfg_initialized) return;

    int maxdepth = maxn - 2;
    if (maxdepth <= 0) maxdepth = 1;

    if (fs) {
        for (int i = 0; i < maxdepth; i++) {
            if (fs[i].points) {
                for (int j = 0; j < maxm; j++) {
                    if (fs[i].points[j].objectives) {
                        free(fs[i].points[j].objectives);
                    }
                }
                free(fs[i].points);
            }
        }
        free(fs);
        fs = NULL;
    }

    if (nextp) {
        for (int i = 0; i < maxdepth; i++) {
            if (nextp[i]) free(nextp[i]);
        }
        free(nextp);
        nextp = NULL;
    }

    if (prevp) {
        for (int i = 0; i < maxdepth; i++) {
            if (prevp[i]) free(prevp[i]);
        }
        free(prevp);
        prevp = NULL;
    }

    if (firstp) {
        free(firstp);
        firstp = NULL;
    }

    if (lastp) {
        free(lastp);
        lastp = NULL;
    }

    if (psize) {
        free(psize);
        psize = NULL;
    }

    if (ref.objectives) {
        free(ref.objectives);
        ref.objectives = NULL;
    }

    wfg_initialized = 0;
}
