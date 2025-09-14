#include "wfg_simple_interface.h"
#include <iostream>
#include <vector>

int main() {
    std::cout << "Testing Simple WFG Interface" << std::endl;

    // Initialize
    if (!wfg_simple_init(2, 100)) {
        std::cerr << "Failed to initialize" << std::endl;
        return 1;
    }

    // Test case 1: Simple 2D points
    std::cout << "\n=== Test 1: Simple 2D Points ===" << std::endl;

    // Points (flattened): {1,5}, {2,4}, {3,3}, {4,2}, {5,1}
    double points[] = {1.0, 5.0, 2.0, 4.0, 3.0, 3.0, 4.0, 2.0, 5.0, 1.0};
    double nadir[] = {6.0, 6.0};

    double hv1 = wfg_simple_calculate(points, nadir, 5, 2);
    std::cout << "Hypervolume: " << hv1 << std::endl;

    // Test case 2: Kursawe-like data (negative values)
    std::cout << "\n=== Test 2: Kursawe-like Data ===" << std::endl;

    double kur_points[] = {
        -7.67, -5.41,
        -7.80, -4.47,
        -8.76, -0.68,
        -7.56, -6.22,
        -7.37, -7.36
    };
    double kur_nadir[] = {0.0, 0.0};

    double hv2 = wfg_simple_calculate(kur_points, kur_nadir, 5, 2);
    std::cout << "Hypervolume (Kursawe-like): " << hv2 << std::endl;

    // Test case 3: Single point
    std::cout << "\n=== Test 3: Single Point ===" << std::endl;

    double single_point[] = {-8.0, -6.0};
    double hv3 = wfg_simple_calculate(single_point, kur_nadir, 1, 2);
    std::cout << "Hypervolume (single): " << hv3 << std::endl;

    wfg_simple_cleanup();

    std::cout << "\nAll tests completed!" << std::endl;
    return 0;
}
