#include "wfg_cpp_interface.h"
#include <iostream>
#include <vector>
#include <cmath>

int main() {
    std::cout << "Testing WFG C++ Interface" << std::endl;

    // Create calculator
    WFGHypervolumeCalculator calc;

    // Test case 1: Simple 2D example
    std::cout << "\n=== Test 1: Simple 2D Example ===" << std::endl;

    // Initialize for 2D, max 100 points
    if (!calc.PrepareWFG(2, 100)) {
        std::cerr << "Failed to initialize WFG" << std::endl;
        return 1;
    }

    // Create test points (minimization problem)
    std::vector<std::vector<double>> points = {
        {1.0, 5.0},
        {2.0, 4.0},
        {3.0, 3.0},
        {4.0, 2.0},
        {5.0, 1.0}
    };

    // Reference point (nadir)
    std::vector<double> nadir = {6.0, 6.0};

    double hv = calc.CalculateHypervolume(points, nadir);
    std::cout << "Hypervolume: " << hv << std::endl;

    // Test case 2: Using C-style interface
    std::cout << "\n=== Test 2: C-style Interface ===" << std::endl;

    // Flatten points
    double pts_flat[] = {1.0, 5.0, 2.0, 4.0, 3.0, 3.0, 4.0, 2.0, 5.0, 1.0};
    double nadir_array[] = {6.0, 6.0};

    double hv2 = calc.CalculateHypervolume(pts_flat, nadir_array, 5, 2);
    std::cout << "Hypervolume (C-style): " << hv2 << std::endl;

    // Test case 3: Test with actual Kursawe-like data (negative values)
    std::cout << "\n=== Test 3: Kursawe-like Data ===" << std::endl;

    std::vector<std::vector<double>> kur_points = {
        {-7.67, -5.41},
        {-7.80, -4.47},
        {-8.76, -0.68},
        {-7.56, -6.22},
        {-7.37, -7.36}
    };

    std::vector<double> kur_nadir = {0.0, 0.0}; // Origin as reference

    double hv3 = calc.CalculateHypervolume(kur_points, kur_nadir);
    std::cout << "Hypervolume (Kursawe-like): " << hv3 << std::endl;

    // Test case 4: C interface
    std::cout << "\n=== Test 4: Pure C Interface ===" << std::endl;

    WFGHypervolumeCalculator* c_calc = wfg_create();
    if (wfg_prepare(c_calc, 2, 100)) {
        double hv4 = wfg_calculate(c_calc, pts_flat, nadir_array, 5, 2);
        std::cout << "Hypervolume (C interface): " << hv4 << std::endl;
    }
    wfg_destroy(c_calc);

    std::cout << "\nTest completed!" << std::endl;
    return 0;
}
