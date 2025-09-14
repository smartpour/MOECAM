#ifndef WFG_CPP_INTERFACE_H
#define WFG_CPP_INTERFACE_H

#include <vector>
#include <memory>

extern "C" {
    #include "wfg.h"
}

class WFGHypervolumeCalculator {
private:
    int m_maxPoints;
    int m_maxDim;
    bool m_initialized;

    // WFG internal structures
    FRONT* fs;
    int** nextp;
    int** prevp;
    int* firstp;
    int* lastp;
    int* psize;
    POINT ref;

    // Global variables needed by WFG algorithm
    extern int n;
    extern POINT ref_global;
    extern FRONT* fs_global;
    extern int fr;
    extern int maxm;
    extern int maxn;
    extern int safe;
    extern int** nextp_global;
    extern int** prevp_global;
    extern int* firstp_global;
    extern int* lastp_global;
    extern int* psize_global;

public:
    WFGHypervolumeCalculator();
    ~WFGHypervolumeCalculator();

    // Initialize WFG with maximum expected dimensions and points
    bool PrepareWFG(int dimensions, int maxPoints);

    // Calculate hypervolume for a single front
    // pts: flattened array of points (nPoints * dimensions)
    // nadir: reference point (dimensions elements)
    // nPoints: number of points
    // dimensions: number of objectives
    double CalculateHypervolume(const double* pts, const double* nadir,
                               int nPoints, int dimensions);

    // Convenience method for std::vector
    double CalculateHypervolume(const std::vector<std::vector<double>>& points,
                               const std::vector<double>& nadir);

    // Main entry point similar to your code
    double MainEntryWFG(int nPoints, int dimensions, FRONT* front);

    // Alternative simple interface
    double MainEntryWFG(int nPoints, int dimensions, const double* nadir, const double* pts);

private:
    void Cleanup();
    void SetupGlobalVariables();
    FRONT CreateFrontFromArray(const double* pts, const double* nadir,
                              int nPoints, int dimensions);
    void FreeFront(FRONT& front);
};

// C-style interface for easier Python binding
extern "C" {
    // Create/destroy calculator instance
    WFGHypervolumeCalculator* wfg_create();
    void wfg_destroy(WFGHypervolumeCalculator* calc);

    // Initialize
    int wfg_prepare(WFGHypervolumeCalculator* calc, int dimensions, int maxPoints);

    // Calculate hypervolume
    double wfg_calculate(WFGHypervolumeCalculator* calc,
                        const double* pts, const double* nadir,
                        int nPoints, int dimensions);
}

#endif // WFG_CPP_INTERFACE_H
