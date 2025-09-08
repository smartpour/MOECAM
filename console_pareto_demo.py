#!/usr/bin/env python
"""Console-only Pareto front demonstration - no GUI dependencies."""

import numpy as np
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), 'src', 'src'))
from moecam.core.cffi_interface import CppOptimizer

def pareto_front_2d(objectives):
    """Find Pareto front for 2D objectives (minimization)."""
    n_points = objectives.shape[0]
    is_pareto = np.ones(n_points, dtype=bool)

    for i in range(n_points):
        if is_pareto[i]:
            dominated = np.all(objectives <= objectives[i], axis=1) & \
                       np.any(objectives < objectives[i], axis=1)
            is_pareto[dominated] = False

    return is_pareto

# Test objectives with point recording
all_points = []

def sphere_with_recording(x):
    """Sphere objectives that records every evaluation."""
    f = np.array([np.sum(x**2), np.sum((x-1)**2)])
    all_points.append((x.copy(), f.copy()))  # Record point
    if len(all_points) % 100 == 0:  # Progress indicator
        print(f"  Evaluated {len(all_points)} points...")
    return f

def main():
    global all_points
    print("=== Console Pareto Front Analysis ===")
    print("Recording all objective evaluations...")

    # Run optimization with recording
    opt = CppOptimizer(3, 2, [-2, -2, -2], [2, 2, 2], sphere_with_recording)
    X, F = opt.sample_random(500)

    # Extract recorded data
    recorded_X = np.array([point[0] for point in all_points])
    recorded_F = np.array([point[1] for point in all_points])

    print(f"\nRecorded {len(recorded_F)} objective evaluations")

    # Find Pareto front
    pareto_mask = pareto_front_2d(recorded_F)
    pareto_X = recorded_X[pareto_mask]
    pareto_F = recorded_F[pareto_mask]

    # Print results
    print(f"\n=== RESULTS ===")
    print(f"Total points evaluated: {len(recorded_F)}")
    print(f"Pareto optimal points: {len(pareto_F)}")
    print(f"Pareto efficiency: {len(pareto_F)/len(recorded_F)*100:.1f}%")

    # Show some example points
    print(f"\n=== SAMPLE POINTS ===")
    print("Decision vars (x1, x2, x3) -> Objectives (f1, f2)")
    for i in range(min(10, len(recorded_F))):
        x = recorded_X[i]
        f = recorded_F[i]
        print(f"({x[0]:6.3f}, {x[1]:6.3f}, {x[2]:6.3f}) -> ({f[0]:8.3f}, {f[1]:8.3f})")

    # Show Pareto points
    print(f"\n=== PARETO OPTIMAL POINTS ===")
    if len(pareto_F) > 0:
        # Sort by f1 for cleaner display
        sorted_idx = np.argsort(pareto_F[:, 0])
        for i in sorted_idx:
            x = pareto_X[i]
            f = pareto_F[i]
            print(f"({x[0]:6.3f}, {x[1]:6.3f}, {x[2]:6.3f}) -> ({f[0]:8.3f}, {f[1]:8.3f}) *PARETO*")
    else:
        print("No Pareto points found!")

    # Save to file
    filename = "recorded_objectives.txt"
    header = f"x1 x2 x3 f1 f2 pareto\nRecorded from sphere optimization ({len(recorded_F)} points)"
    data = np.column_stack([recorded_X, recorded_F, pareto_mask.astype(int)])
    np.savetxt(filename, data, header=header,
               fmt=['%.6f', '%.6f', '%.6f', '%.6f', '%.6f', '%d'])
    print(f"\nSaved all points to {filename}")
    print("Format: x1 x2 x3 f1 f2 pareto(1=yes,0=no)")

    return recorded_X, recorded_F, pareto_mask

if __name__ == "__main__":
    X, F, pareto = main()

    # Additional analysis
    print(f"\n=== OBJECTIVE SPACE STATISTICS ===")
    print(f"f1 range: [{F[:,0].min():.3f}, {F[:,0].max():.3f}]")
    print(f"f2 range: [{F[:,1].min():.3f}, {F[:,1].max():.3f}]")

    if np.sum(pareto) > 0:
        pf = F[pareto]
        print(f"Pareto f1 range: [{pf[:,0].min():.3f}, {pf[:,0].max():.3f}]")
        print(f"Pareto f2 range: [{pf[:,1].min():.3f}, {pf[:,1].max():.3f}]")

    print("\nDone! Check recorded_objectives.txt for detailed results.")
