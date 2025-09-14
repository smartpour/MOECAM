"""
Example usage of MOECAM algorithms from the proper package structure.

This demonstrates how to use the MOECAM algorithms after moving them
to the correct package location.
"""

import sys
import os

# Add the src directory to path to import moecam
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

try:
    # Import from the proper package structure
    from moecam import minimize_ecam, minimize_random_start, direct_optimize
    from moecam import extract_pareto_front, calculate_hypervolume
    print("✅ Successfully imported MOECAM from package structure!")
except ImportError as e:
    print(f"❌ Import failed: {e}")
    sys.exit(1)

def main():
    """Demo the MOECAM algorithms from proper package structure."""

    print("\nMOECAM Package Structure Demo")
    print("=" * 40)

    # Define test functions
    def bi_objective(x):
        """Simple bi-objective problem."""
        f1 = sum(xi**2 for xi in x)
        f2 = sum((xi - 1)**2 for xi in x)
        return [f1, f2]

    def sphere(x):
        """Single objective sphere function."""
        return [sum(xi**2 for xi in x)]

    bounds = [(-2.0, 2.0), (-2.0, 2.0)]

    print("\n1. Testing ECAM algorithm from moecam.algorithms")
    try:
        result = minimize_ecam(bi_objective, bounds, num_objectives=2, max_iterations=50)
        print(f"   Success: {result['success']}")
        if result['success']:
            print(f"   Solution: [{', '.join(f'{x:.4f}' for x in result['x'])}]")
            print(f"   Objectives: [{', '.join(f'{f:.4f}' for f in result['f'])}]")
            print(f"   Function evaluations: {result['function_evaluations']}")
    except Exception as e:
        print(f"   Error: {e}")

    print("\n2. Testing DIRECT algorithm from moecam.algorithms")
    try:
        result = direct_optimize(sphere, bounds, max_function_evaluations=30)
        print(f"   Success: {result['success']}")
        if result['success']:
            print(f"   Solution: [{', '.join(f'{x:.4f}' for x in result['x'])}]")
            print(f"   Objective: {result['f'][0]:.6f}")
            print(f"   Function evaluations: {result['function_evaluations']}")
    except Exception as e:
        print(f"   Error: {e}")

    print("\n3. Testing Pareto tools from moecam.core")
    try:
        # Generate some test points
        test_objectives = [
            [1.0, 4.0],
            [2.0, 3.0],
            [3.0, 2.0],
            [4.0, 1.0],
            [2.5, 2.5]  # This should be dominated
        ]

        pareto_front = extract_pareto_front(test_objectives)
        print(f"   Input points: {len(test_objectives)}")
        print(f"   Pareto points: {len(pareto_front)}")

        if len(pareto_front) > 0:
            hypervolume = calculate_hypervolume(pareto_front, [5.0, 5.0])
            print(f"   Hypervolume: {hypervolume:.6f}")
    except Exception as e:
        print(f"   Error: {e}")

    print("\n" + "=" * 40)
    print("✅ MOECAM Package Structure Working!")
    print("\nNow you can import MOECAM algorithms like this:")
    print("from moecam import minimize_ecam, minimize_random_start, direct_optimize")
    print("from moecam import extract_pareto_front, calculate_hypervolume")
    print("\nOr from specific subpackages:")
    print("from moecam.algorithms import minimize_ecam")
    print("from moecam.core import extract_pareto_front")

if __name__ == "__main__":
    main()
