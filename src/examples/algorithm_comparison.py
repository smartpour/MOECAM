#!/usr/bin/env python3
"""
Algorithm comparison example for MOECAM package.

This example demonstrates how to compare different algorithms
on the same problem and evaluate their performance.
"""

import numpy as np
import matplotlib.pyplot as plt
from src.moecam.problems.test_functions import ZDT1
from src.moecam.algorithms.moea_algorithms import NSGAII, MOEAD
from src.moecam.metrics.performance_metrics import hypervolume, EvaluationCounter

def run_algorithm_comparison():
    print("MOECAM Algorithm Comparison Example")
    print("=" * 50)
    
    # Define the problem
    problem = ZDT1(n_dim=30)
    print(f"Problem: ZDT1 with {problem.n_dim} variables and {problem.n_obj} objectives")
    
    # Define algorithms to compare
    algorithms = {
        'NSGA-II': NSGAII(problem, pop_size=100, num_generations=250),
        'MOEA/D': MOEAD(problem, pop_size=100, num_generations=250)
    }
    
    results = {}
    reference_point = [1.1, 1.1]
    
    # Run each algorithm
    for name, algorithm in algorithms.items():
        print(f"\nRunning {name}...")
        
        # Initialize evaluation counter
        counter = EvaluationCounter()
        counter.start()
        
        # Run optimization
        pareto_front_result = algorithm.optimize()
        
        counter.stop()
        
        # Calculate metrics
        hv = hypervolume(pareto_front_result, reference_point)
        
        results[name] = {
            'pareto_front': pareto_front_result,
            'hypervolume': hv,
            'evaluations': counter.get_count(),
            'time': counter.get_time()
        }
        
        print(f"{name} Results:")
        print(f"  - Pareto solutions: {len(pareto_front_result)}")
        print(f"  - Hypervolume: {hv:.4f}")
        print(f"  - Evaluations: {counter.get_count()}")
        print(f"  - Time: {counter.get_time():.2f}s")
    
    # Visualize comparison
    plt.figure(figsize=(15, 5))
    
    # Plot Pareto fronts
    plt.subplot(1, 3, 1)
    for name, result in results.items():
        pf = result['pareto_front']
        plt.scatter(pf[:, 0], pf[:, 1], alpha=0.7, s=30, label=name)
    plt.xlabel('Objective 1')
    plt.ylabel('Objective 2')
    plt.title('Pareto Fronts Comparison')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot hypervolume comparison
    plt.subplot(1, 3, 2)
    names = list(results.keys())
    hvs = [results[name]['hypervolume'] for name in names]
    plt.bar(names, hvs, alpha=0.7)
    plt.ylabel('Hypervolume')
    plt.title('Hypervolume Comparison')
    plt.grid(True, alpha=0.3)
    
    # Plot time comparison
    plt.subplot(1, 3, 3)
    times = [results[name]['time'] for name in names]
    plt.bar(names, times, alpha=0.7)
    plt.ylabel('Time (seconds)')
    plt.title('Execution Time Comparison')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('algorithm_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("\nComparison results saved to 'algorithm_comparison.png'")

if __name__ == "__main__":
    run_algorithm_comparison()

