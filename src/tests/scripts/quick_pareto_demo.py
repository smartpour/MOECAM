#!/usr/bin/env python
"""Quick Pareto front demonstration using random sampling and objective recording."""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add package path from new location
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from moecam.core.cffi_interface import CppOptimizer

def pareto_front_2d(objectives):
    """Extract Pareto front from 2D objectives (minimization).

    Args:
        objectives: array of shape (n_points, 2)

    Returns:
        pareto_mask: boolean array indicating Pareto optimal points
    """
    n_points = objectives.shape[0]
    is_pareto = np.ones(n_points, dtype=bool)

    for i in range(n_points):
        if is_pareto[i]:
            # Point i dominates point j if it's better in all objectives
            dominated = np.all(objectives <= objectives[i], axis=1) & \
                       np.any(objectives < objectives[i], axis=1)
            is_pareto[dominated] = False

    return is_pareto

def plot_pareto_analysis(X, F, title="Pareto Front Analysis"):
    """Plot objective space with Pareto front highlighted."""
    if F.shape[1] != 2:
        print(f"Can only plot 2D objectives, got {F.shape[1]}")
        return

    # Find Pareto front
    pareto_mask = pareto_front_2d(F)
    pareto_points = F[pareto_mask]

    # Sort Pareto points for nice line plot
    sorted_idx = np.argsort(pareto_points[:, 0])
    pareto_sorted = pareto_points[sorted_idx]

    plt.figure(figsize=(10, 4))

    # Plot 1: Objective space
    plt.subplot(1, 2, 1)
    plt.scatter(F[:, 0], F[:, 1], alpha=0.6, s=20, label=f'All points ({len(F)})')
    plt.scatter(pareto_points[:, 0], pareto_points[:, 1],
                color='red', s=40, label=f'Pareto front ({len(pareto_points)})')
    plt.plot(pareto_sorted[:, 0], pareto_sorted[:, 1], 'r-', alpha=0.7, linewidth=2)
    plt.xlabel('f1 (minimize)')
    plt.ylabel('f2 (minimize)')
    plt.title('Objective Space')
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Plot 2: Decision space (first 2 variables if available)
    plt.subplot(1, 2, 2)
    if X.shape[1] >= 2:
        plt.scatter(X[:, 0], X[:, 1], alpha=0.6, s=20, label='All points')
        plt.scatter(X[pareto_mask, 0], X[pareto_mask, 1],
                    color='red', s=40, label='Pareto optimal')
        plt.xlabel('x1')
        plt.ylabel('x2')
        plt.title('Decision Space (x1, x2)')
        plt.legend()
        plt.grid(True, alpha=0.3)
    else:
        plt.text(0.5, 0.5, 'Need ≥2 variables\nfor decision plot',
                ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('Decision Space')

    plt.suptitle(title)
    plt.tight_layout()

    # Print statistics
    print(f"\n=== {title} ===")
    print(f"Total points evaluated: {len(F)}")
    print(f"Pareto optimal points: {len(pareto_points)}")
    print(f"Pareto efficiency: {len(pareto_points)/len(F)*100:.1f}%")

    if len(pareto_points) > 0:
        print(f"f1 range: [{pareto_points[:,0].min():.3f}, {pareto_points[:,0].max():.3f}]")
        print(f"f2 range: [{pareto_points[:,1].min():.3f}, {pareto_points[:,1].max():.3f}]")

    return pareto_mask, pareto_points

# Test problems with known Pareto characteristics
def sphere_objectives(x):
    """Bi-objective: sphere and shifted sphere (convex Pareto front)."""
    return np.array([np.sum(x**2), np.sum((x-1)**2)])

def zdt1_objectives(x):
    """ZDT1: classic test problem with known Pareto front."""
    f1 = x[0]
    g = 1 + 9 * np.sum(x[1:]) / (len(x) - 1)
    f2 = g * (1 - np.sqrt(f1/g))
    return np.array([f1, f2])

def kursawe_objectives(x):
    """Kursawe: multi-modal test problem."""
    f1 = np.sum(-10 * np.exp(-0.2 * np.sqrt(x[:-1]**2 + x[1:]**2)))
    f2 = np.sum(np.abs(x)**0.8 + 5 * np.sin(x**3))
    return np.array([f1, f2])

if __name__ == "__main__":
    print("=== Quick Pareto Front Analysis Using Random Sampling ===")

    # Test 1: Simple sphere problem
    print("\n1. Testing sphere objectives...")
    opt1 = CppOptimizer(3, 2, [-2, -2, -2], [2, 2, 2], sphere_objectives)
    X1, F1 = opt1.sample_random(500)
    pareto1, pf1 = plot_pareto_analysis(X1, F1, "Sphere Objectives")

    # Test 2: ZDT1 problem
    print("\n2. Testing ZDT1 objectives...")
    opt2 = CppOptimizer(5, 2, [0, 0, 0, 0, 0], [1, 1, 1, 1, 1], zdt1_objectives)
    X2, F2 = opt2.sample_random(1000)
    pareto2, pf2 = plot_pareto_analysis(X2, F2, "ZDT1 Test Problem")

    # Test 3: Kursawe problem
    print("\n3. Testing Kursawe objectives...")
    opt3 = CppOptimizer(3, 2, [-5, -5, -5], [5, 5, 5], kursawe_objectives)
    X3, F3 = opt3.sample_random(800)
    pareto3, pf3 = plot_pareto_analysis(X3, F3, "Kursawe Problem")

    # Save all objective points to file for further analysis
    results_dir = os.path.join(os.path.dirname(__file__), '..', 'test_results')
    os.makedirs(results_dir, exist_ok=True)

    np.savetxt(os.path.join(results_dir, 'sphere_objectives.txt'), F1, header='f1 f2 (sphere)')
    np.savetxt(os.path.join(results_dir, 'zdt1_objectives.txt'), F2, header='f1 f2 (ZDT1)')
    np.savetxt(os.path.join(results_dir, 'kursawe_objectives.txt'), F3, header='f1 f2 (Kursawe)')

    print(f"\n=== Saved objective values to test_results/ ===")
    print("Files: sphere_objectives.txt, zdt1_objectives.txt, kursawe_objectives.txt")

    plt.show()
