#include "wfg_cpp_interface.h"
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>

// We need to include the WFG C code and handle the global variables
extern "C" {
    // WFG global variables that need to be accessible
    extern int n;
    extern POINT ref;
    extern FRONT* fs;
    extern int fr;
    extern int maxm;
    extern int maxn;
    extern int safe;
    extern int** nextp;
    extern int** prevp;
    extern int* firstp;
    extern int* lastp;
    extern int* psize;

    // WFG functions we need
    extern double hv(FRONT ps);
}

WFGHypervolumeCalculator::WFGHypervolumeCalculator()
    : m_maxPoints(0), m_maxDim(0), m_initialized(false), fs(nullptr),
      nextp(nullptr), prevp(nullptr), firstp(nullptr), lastp(nullptr), psize(nullptr) {
}

WFGHypervolumeCalculator::~WFGHypervolumeCalculator() {
    Cleanup();
}

bool WFGHypervolumeCalculator::PrepareWFG(int dimensions, int maxPoints) {
    if (m_initialized) {
        Cleanup();
    }

    m_maxDim = dimensions;
    m_maxPoints = maxPoints;

    // Set global variables
    maxm = maxPoints;
    maxn = dimensions;

    try {
        // Allocate memory like in the original WFG code
        int maxdepth = maxn - 2;
        if (maxdepth <= 0) maxdepth = 1; // Handle edge case

        fs = (FRONT*)malloc(sizeof(FRONT) * maxdepth);
        for (int i = 0; i < maxdepth; i++) {
            fs[i].points = (POINT*)malloc(sizeof(POINT) * maxm);
            for (int j = 0; j < maxm; j++) {
                fs[i].points[j].objectives = (OBJECTIVE*)malloc(sizeof(OBJECTIVE) * (maxn - i));
            }
        }

        nextp = (int**)malloc(sizeof(int*) * maxdepth);
        for (int i = 0; i < maxdepth; i++) {
            nextp[i] = (int*)malloc(sizeof(int) * maxm);
        }

        prevp = (int**)malloc(sizeof(int*) * maxdepth);
        for (int i = 0; i < maxdepth; i++) {
            prevp[i] = (int*)malloc(sizeof(int) * maxm);
        }

        firstp = (int*)malloc(sizeof(int) * maxdepth);
        lastp = (int*)malloc(sizeof(int) * maxdepth);
        psize = (int*)malloc(sizeof(int) * maxdepth);

        // Initialize reference point
        ref.objectives = (OBJECTIVE*)malloc(sizeof(OBJECTIVE) * maxn);

        m_initialized = true;
        return true;

    } catch (...) {
        Cleanup();
        return false;
    }
}

void WFGHypervolumeCalculator::Cleanup() {
    if (!m_initialized) return;

    if (fs) {
        int maxdepth = std::max(1, maxn - 2);
        for (int i = 0; i < maxdepth; i++) {
            if (fs[i].points) {
                for (int j = 0; j < maxm; j++) {
                    free(fs[i].points[j].objectives);
                }
                free(fs[i].points);
            }
        }
        free(fs);
        fs = nullptr;
    }

    if (nextp) {
        int maxdepth = std::max(1, maxn - 2);
        for (int i = 0; i < maxdepth; i++) {
            free(nextp[i]);
        }
        free(nextp);
        nextp = nullptr;
    }

    if (prevp) {
        int maxdepth = std::max(1, maxn - 2);
        for (int i = 0; i < maxdepth; i++) {
            free(prevp[i]);
        }
        free(prevp);
        prevp = nullptr;
    }

    free(firstp); firstp = nullptr;
    free(lastp); lastp = nullptr;
    free(psize); psize = nullptr;
    free(ref.objectives); ref.objectives = nullptr;

    m_initialized = false;
}

FRONT WFGHypervolumeCalculator::CreateFrontFromArray(const double* pts, const double* nadir,
                                                    int nPoints, int dimensions) {
    FRONT front;
    front.n = dimensions;
    front.nPoints = nPoints;
    front.points = (POINT*)malloc(sizeof(POINT) * nPoints);

    for (int j = 0; j < nPoints; j++) {
        front.points[j].objectives = (OBJECTIVE*)malloc(sizeof(OBJECTIVE) * dimensions);
        for (int k = 0; k < dimensions; k++) {
            // Apply the same transformation as in the original WFG code:
            // fabs(objective - reference_point)
            double objective_value = pts[j * dimensions + k];
            double reference_value = nadir[k];
            front.points[j].objectives[k] = fabs(objective_value - reference_value);
        }
    }

    return front;
}

void WFGHypervolumeCalculator::FreeFront(FRONT& front) {
    if (front.points) {
        for (int i = 0; i < front.nPoints; i++) {
            free(front.points[i].objectives);
        }
        free(front.points);
        front.points = nullptr;
    }
}

double WFGHypervolumeCalculator::CalculateHypervolume(const double* pts, const double* nadir,
                                                     int nPoints, int dimensions) {
    if (!m_initialized) {
        std::cerr << "WFG not initialized. Call PrepareWFG first." << std::endl;
        return 0.0;
    }

    if (nPoints <= 0 || dimensions <= 0) {
        return 0.0;
    }

    // Create front from input data
    FRONT front = CreateFrontFromArray(pts, nadir, nPoints, dimensions);

    // Set global variables for WFG algorithm
    n = dimensions;
    fr = 0;
    safe = 0;

    // Copy reference point
    for (int i = 0; i < dimensions; i++) {
        ref.objectives[i] = nadir[i];
    }

    // Calculate hypervolume
    double hypervolume = hv(front);

    // Cleanup
    FreeFront(front);

    return hypervolume;
}

double WFGHypervolumeCalculator::CalculateHypervolume(const std::vector<std::vector<double>>& points,
                                                     const std::vector<double>& nadir) {
    if (points.empty() || nadir.empty()) {
        return 0.0;
    }

    int nPoints = points.size();
    int dimensions = nadir.size();

    // Flatten points array
    std::vector<double> pts_flat;
    pts_flat.reserve(nPoints * dimensions);

    for (const auto& point : points) {
        if (point.size() != dimensions) {
            std::cerr << "Point dimension mismatch" << std::endl;
            return 0.0;
        }
        for (double val : point) {
            pts_flat.push_back(val);
        }
    }

    return CalculateHypervolume(pts_flat.data(), nadir.data(), nPoints, dimensions);
}

double WFGHypervolumeCalculator::MainEntryWFG(int nPoints, int dimensions, FRONT* front) {
    if (!m_initialized) {
        return 0.0;
    }

    // Set global variables
    n = dimensions;
    fr = 0;
    safe = 0;

    return hv(*front);
}

double WFGHypervolumeCalculator::MainEntryWFG(int nPoints, int dimensions,
                                             const double* nadir, const double* pts) {
    return CalculateHypervolume(pts, nadir, nPoints, dimensions);
}

// C interface for Python binding
extern "C" {
    WFGHypervolumeCalculator* wfg_create() {
        return new WFGHypervolumeCalculator();
    }

    void wfg_destroy(WFGHypervolumeCalculator* calc) {
        delete calc;
    }

    int wfg_prepare(WFGHypervolumeCalculator* calc, int dimensions, int maxPoints) {
        if (!calc) return 0;
        return calc->PrepareWFG(dimensions, maxPoints) ? 1 : 0;
    }

    double wfg_calculate(WFGHypervolumeCalculator* calc,
                        const double* pts, const double* nadir,
                        int nPoints, int dimensions) {
        if (!calc) return 0.0;
        return calc->CalculateHypervolume(pts, nadir, nPoints, dimensions);
    }
}
