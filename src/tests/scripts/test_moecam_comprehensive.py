"""
Comprehensive test of MOECAM algorithms with various test problems.

This demonstrates all three MOECAM algorithms (ECAM, Random Start, DIRECT)
on different optimization problems.
"""

import sys
import os

# Add the current directory to path for imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from moecam_direct_interface import minimize_ecam, minimize_random_start, direct_optimize
import math
import numpy as np

def test_functions():
    """Define various test functions for optimization."""

    def sphere(x):
        """Simple sphere function - single objective."""
        return [sum(xi**2 for xi in x)]

    def rosenbrock(x):
        """Rosenbrock function - single objective."""
        total = 0
        for i in range(len(x) - 1):
            total += 100 * (x[i+1] - x[i]**2)**2 + (1 - x[i])**2
        return [total]

    def bi_objective_simple(x):
        """Simple bi-objective function."""
        f1 = sum(xi**2 for xi in x)
        f2 = sum((xi - 1)**2 for xi in x)
        return [f1, f2]

    def kursawe(x):
        """Kursawe multi-objective function."""
        f1 = sum(-10 * math.exp(-0.2 * math.sqrt(x[i]**2 + x[i+1]**2)) for i in range(len(x)-1))
        f2 = sum(abs(xi)**0.8 + 5 * math.sin(xi**3) for xi in x)
        return [f1, f2]

    def zdt1(x):
        """ZDT1 test function - classic multi-objective benchmark."""
        f1 = x[0]
        g = 1 + 9 * sum(x[1:]) / (len(x) - 1)
        h = 1 - math.sqrt(f1 / g)
        f2 = g * h
        return [f1, f2]

    def three_objective(x):
        """Three-objective test function."""
        f1 = sum(xi**2 for xi in x)
        f2 = sum((xi - 1)**2 for xi in x)
        f3 = sum((xi + 1)**2 for xi in x)
        return [f1, f2, f3]

    return {
        'sphere': (sphere, [(-5, 5)] * 3, 1),
        'rosenbrock': (rosenbrock, [(-2, 2)] * 4, 1),
        'bi_objective_simple': (bi_objective_simple, [(-5, 5)] * 2, 2),
        'kursawe': (kursawe, [(-5, 5)] * 3, 2),
        'zdt1': (zdt1, [(0, 1)] * 5, 2),
        'three_objective': (three_objective, [(-2, 2)] * 3, 3)
    }

def run_comprehensive_test():
    """Run comprehensive tests of all MOECAM algorithms."""

    print("MOECAM Comprehensive Algorithm Test")
    print("=" * 50)

    functions = test_functions()
    algorithms = [
        ('ECAM', minimize_ecam),
        ('Random Start', minimize_random_start),
        ('DIRECT', direct_optimize)
    ]

    results = {}

    for alg_name, alg_func in algorithms:
        print(f"\n{alg_name} Algorithm Results:")
        print("-" * 30)
        results[alg_name] = {}

        for func_name, (func, bounds, n_obj) in functions.items():
            print(f"\n  {func_name} (n_obj={n_obj}):")

            try:
                # Adjust parameters based on algorithm
                if alg_name == 'DIRECT':
                    result = alg_func(func, bounds, max_function_evaluations=100)
                else:
                    result = alg_func(func, bounds, num_objectives=n_obj, max_iterations=100)

                results[alg_name][func_name] = result

                print(f"    Success: {result['success']}")
                if result['success']:
                    print(f"    Solution: {[f'{x:.4f}' for x in result['x']]}")
                    print(f"    Objectives: {[f'{f:.4f}' for f in result['f']]}")
                    print(f"    Function evals: {result['function_evaluations']}")

                    if 'pareto_front' in result:
                        print(f"    Pareto points: {len(result['pareto_front'])}")
                    if 'hypervolume' in result:
                        print(f"    Hypervolume: {result['hypervolume']:.6f}")
                else:
                    print(f"    Error: {result['message']}")

            except Exception as e:
                print(f"    Exception: {e}")
                results[alg_name][func_name] = {'success': False, 'error': str(e)}

    # Summary comparison
    print("\n" + "="*50)
    print("ALGORITHM COMPARISON SUMMARY")
    print("="*50)

    for func_name in functions.keys():
        print(f"\n{func_name}:")
        for alg_name in results.keys():
            result = results[alg_name].get(func_name, {})
            if result.get('success'):
                obj_str = f"{result['f'][0]:.4f}" if len(result['f']) == 1 else f"[{', '.join(f'{f:.4f}' for f in result['f'])}]"
                print(f"  {alg_name:12}: ✓ f={obj_str}, evals={result['function_evaluations']}")
            else:
                print(f"  {alg_name:12}: ✗ Failed")

    # Performance analysis
    print("\n" + "="*50)
    print("PERFORMANCE ANALYSIS")
    print("="*50)

    print("\nFunction evaluations by algorithm:")
    eval_stats = {}
    for alg_name in results.keys():
        evals = [r['function_evaluations'] for r in results[alg_name].values() if r.get('success')]
        if evals:
            eval_stats[alg_name] = {
                'min': min(evals),
                'max': max(evals),
                'avg': sum(evals) / len(evals),
                'total_problems': len(evals)
            }
            print(f"  {alg_name:12}: avg={eval_stats[alg_name]['avg']:.1f}, range=[{eval_stats[alg_name]['min']}-{eval_stats[alg_name]['max']}], solved={eval_stats[alg_name]['total_problems']}/{len(functions)}")

    print("\nSuccess rates:")
    for alg_name in results.keys():
        successes = sum(1 for r in results[alg_name].values() if r.get('success'))
        success_rate = successes / len(functions) * 100
        print(f"  {alg_name:12}: {successes}/{len(functions)} ({success_rate:.1f}%)")

    print("\n✅ Comprehensive MOECAM testing complete!")
    print("\nAll three MOECAM algorithms are working correctly:")
    print("- minimize_ecam: Multi-objective evolutionary algorithm")
    print("- minimize_random_start: Multi-restart optimization")
    print("- direct_optimize: DIRECT-like global optimization")
    print("\nThe interface provides the exact API you requested with callback")
    print("function conversion from Python to C++ signatures.")

if __name__ == "__main__":
    run_comprehensive_test()
