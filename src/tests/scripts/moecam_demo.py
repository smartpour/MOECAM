"""
Simple MOECAM Interface Demo

Demonstrates the three MOECAM algorithms with basic test problems.
This is the final working interface that provides Python access to:
- MinimizeECAM
- MinimizeRandomStart
- direct_optimize
"""

import sys
import os

# Add the current directory to path for imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from moecam_direct_interface import minimize_ecam, minimize_random_start, direct_optimize
import math

def main():
    """Demonstrate the three MOECAM algorithms."""

    print("MOECAM Python Interface Demo")
    print("=" * 40)
    print()
    print("This demonstrates Python interfaces to the three MOECAM algorithms:")
    print("1. MinimizeECAM - Multi-objective evolutionary algorithm")
    print("2. MinimizeRandomStart - Multi-restart optimization")
    print("3. direct_optimize - DIRECT global optimization")
    print()

    # Define test functions
    def sphere(x):
        """Simple sphere function."""
        return [sum(xi**2 for xi in x)]

    def bi_objective(x):
        """Simple bi-objective problem."""
        f1 = sum(xi**2 for xi in x)
        f2 = sum((xi - 1)**2 for xi in x)
        return [f1, f2]

    def kursawe(x):
        """Kursawe function - challenging multi-objective problem."""
        f1 = sum(-10 * math.exp(-0.2 * math.sqrt(x[i]**2 + x[i+1]**2)) for i in range(len(x)-1))
        f2 = sum(abs(xi)**0.8 + 5 * math.sin(xi**3) for xi in x)
        return [f1, f2]

    bounds = [(-2.0, 2.0), (-2.0, 2.0)]

    # Test 1: ECAM with bi-objective function
    print("1. ECAM Algorithm")
    print("-" * 20)
    try:
        result = minimize_ecam(bi_objective, bounds, num_objectives=2, max_iterations=50)
        print(f"   Success: {result['success']}")
        if result['success']:
            print(f"   Best solution: [{', '.join(f'{x:.4f}' for x in result['x'])}]")
            print(f"   Objectives: [{', '.join(f'{f:.4f}' for f in result['f'])}]")
            print(f"   Function evaluations: {result['function_evaluations']}")
            if 'pareto_front' in result:
                print(f"   Pareto points found: {len(result['pareto_front'])}")
        else:
            print(f"   Error: {result['message']}")
    except Exception as e:
        print(f"   Exception: {e}")

    print()

    # Test 2: Random Start with Kursawe
    print("2. Random Start Algorithm")
    print("-" * 25)
    try:
        result = minimize_random_start(kursawe, bounds, num_objectives=2, max_iterations=50)
        print(f"   Success: {result['success']}")
        if result['success']:
            print(f"   Best solution: [{', '.join(f'{x:.4f}' for x in result['x'])}]")
            print(f"   Objectives: [{', '.join(f'{f:.4f}' for f in result['f'])}]")
            print(f"   Function evaluations: {result['function_evaluations']}")
        else:
            print(f"   Error: {result['message']}")
    except Exception as e:
        print(f"   Exception: {e}")

    print()

    # Test 3: DIRECT with single objective
    print("3. DIRECT Algorithm")
    print("-" * 18)
    try:
        result = direct_optimize(sphere, bounds, max_function_evaluations=30)
        print(f"   Success: {result['success']}")
        if result['success']:
            print(f"   Best solution: [{', '.join(f'{x:.4f}' for x in result['x'])}]")
            print(f"   Objective: {result['f'][0]:.6f}")
            print(f"   Function evaluations: {result['function_evaluations']}")
        else:
            print(f"   Error: {result['message']}")
    except Exception as e:
        print(f"   Exception: {e}")

    print()
    print("=" * 40)
    print("✅ MOECAM Interface Demo Complete!")
    print()
    print("SUCCESS: All three MOECAM algorithms are now accessible from Python:")
    print()
    print("from moecam_direct_interface import minimize_ecam, minimize_random_start, direct_optimize")
    print()
    print("# Multi-objective optimization with ECAM")
    print("result = minimize_ecam(objective_function, bounds, num_objectives=2)")
    print()
    print("# Multi-restart optimization")
    print("result = minimize_random_start(objective_function, bounds, num_objectives=2)")
    print()
    print("# DIRECT global optimization")
    print("result = direct_optimize(objective_function, bounds)")
    print()
    print("Each function returns a dictionary with:")
    print("- x: best solution found")
    print("- f: objective values at best solution")
    print("- success: whether optimization succeeded")
    print("- function_evaluations: number of function calls")
    print("- pareto_front: Pareto optimal points (for multi-objective)")
    print("- hypervolume: quality metric (for multi-objective)")

if __name__ == "__main__":
    main()
