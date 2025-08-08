#!/usr/bin/env python3
"""
Basic usage example for MOECAM package.

This example demonstrates how to:
1. Define a multi-objective optimization problem
2. Run an optimization algorithm
3. Evaluate the results using performance metrics
4. Visualize the Pareto front
"""

import numpy as np
import matplotlib.pyplot as plt
from src.moecam.problems.test_functions import ZDT1, ZDT2
from src.moecam.algorithms.moea_algorithms import NSGAII
from src.moecam.metrics.performance_metrics import hypervolume, pareto_front

def main():
    print("MOECAM Basic Usage Example")
    print("=" * 40)
    
    # Define the problem
    problem = ZDT1(n_dim=30)
    print(f"Problem: ZDT1 with {problem.n_dim} variables and {problem.n_obj} objectives")
    
    # Initialize the algorithm
    algorithm = NSGAII(
        problem=problem,
        pop_size=100,
        num_generations=250,
        crossover_rate=0.9,
        mutation_rate=0.01
    )
    print(f"Algorithm: NSGA-II with population size {algorithm.pop_size}")
    
    # Run optimization
    print("Running optimization...")
    result_pareto_front = algorithm.optimize()
    print(f"Optimization completed. Found {len(result_pareto_front)} Pareto optimal solutions.")
    
    # Calculate performance metrics
    reference_point = [1.1, 1.1]  # Reference point for hypervolume calculation
    hv = hypervolume(result_pareto_front, reference_point)
    print(f"Hypervolume: {hv:.4f}")
    
    # Visualize results
    plt.figure(figsize=(10, 6))
    plt.scatter(result_pareto_front[:, 0], result_pareto_front[:, 1], 
                alpha=0.7, s=50, label='Pareto Front')
    plt.xlabel('Objective 1')
    plt.ylabel('Objective 2')
    plt.title('ZDT1 Pareto Front - NSGA-II')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('zdt1_pareto_front.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("Results saved to 'zdt1_pareto_front.png'")

if __name__ == "__main__":
    main()

