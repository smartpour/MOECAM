#!/usr/bin/env python3
"""
Test script for MOECAM-based Pareto visualization
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

def test_moecam_visualization():
    """Test the MOECAM algorithms and create basic visualization"""
    
    try:
        # Import MOECAM modules
        from moecam.problems.test_functions import ZDT1
        from moecam.algorithms.moea_algorithms import NSGAII
        from moecam.metrics.performance_metrics import pareto_front, hypervolume
        
        print("✓ MOECAM modules imported successfully")
        
        # Create problem
        problem = ZDT1(n_dim=10)
        print(f"✓ Created {problem.__class__.__name__} problem with {problem.n_dim} dimensions")
        
        # Test with different evaluation budgets
        budgets = [100, 200, 500]
        results = {}
        
        for budget in budgets:
            print(f"\nTesting with {budget} evaluations...")
            
            # Calculate pop_size and generations
            if budget <= 100:
                pop_size = 20
                num_generations = budget // pop_size
            else:
                pop_size = 50
                num_generations = budget // pop_size
            
            num_generations = max(num_generations, 1)
            
            print(f"  Population size: {pop_size}, Generations: {num_generations}")
            
            # Run NSGA-II
            algorithm = NSGAII(problem, pop_size=pop_size, num_generations=num_generations)
            
            # Track all evaluated solutions
            all_solutions = []
            original_evaluate = algorithm._evaluate_population
            
            def tracking_evaluate(population):
                objectives = original_evaluate(population)
                all_solutions.extend(objectives)
                return objectives
            
            algorithm._evaluate_population = tracking_evaluate
            
            # Run optimization
            pareto_solutions = algorithm.optimize()
            
            # Calculate metrics
            reference_point = [1.1, 1.1]
            hv = hypervolume(pareto_solutions, reference_point) if len(pareto_solutions) > 0 else 0.0
            
            results[budget] = {
                'pareto_solutions': pareto_solutions,
                'all_solutions': np.array(all_solutions),
                'hypervolume': hv,
                'n_pareto': len(pareto_solutions),
                'n_total': len(all_solutions)
            }
            
            print(f"  Results: {len(pareto_solutions)} Pareto solutions, "
                  f"{len(all_solutions)} total evaluations, HV = {hv:.4f}")
        
        # Create visualization
        plt.figure(figsize=(12, 8))
        
        # Generate true Pareto front
        f1_true = np.linspace(0, 1, 100)
        f2_true = 1 - np.sqrt(f1_true)
        
        # Plot true Pareto front
        plt.plot(f1_true, f2_true, 'k-', linewidth=2, label='True Pareto Front')
        
        colors = ['red', 'orange', 'blue']
        
        # Plot results for each budget
        for i, budget in enumerate(budgets):
            result = results[budget]
            all_sols = result['all_solutions']
            pareto_sols = result['pareto_solutions']
            
            if len(all_sols) > 0:
                # Plot all evaluated solutions
                plt.scatter(all_sols[:, 0], all_sols[:, 1], 
                           c=colors[i], alpha=0.5, s=20, 
                           label=f'All solutions ({budget} evals)')
                
                # Highlight Pareto solutions
                if len(pareto_sols) > 0:
                    plt.scatter(pareto_sols[:, 0], pareto_sols[:, 1], 
                               c=colors[i], s=50, marker='s', edgecolors='black',
                               label=f'Pareto solutions ({budget} evals)')
        
        plt.xlabel('Objective 1 (f₁)', fontsize=12)
        plt.ylabel('Objective 2 (f₂)', fontsize=12)
        plt.title('MOECAM: ZDT1 Pareto Front Analysis\nNSGA-II with Different Evaluation Budgets', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        # Save plot
        plt.savefig('moecam_pareto_analysis.png', dpi=300, bbox_inches='tight')
        print(f"\n✓ Visualization saved to: moecam_pareto_analysis.png")
        
        plt.show()
        
        # Print summary
        print(f"\n" + "="*50)
        print("MOECAM Analysis Summary:")
        print("="*50)
        for budget in budgets:
            result = results[budget]
            print(f"Budget {budget:4d}: {result['n_pareto']:2d} Pareto, "
                  f"{result['n_total']:3d} total, HV = {result['hypervolume']:.4f}")
        
        return True
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("MOECAM Visualization Test")
    print("=" * 30)
    
    success = test_moecam_visualization()
    
    if success:
        print("\n🎉 Test completed successfully!")
    else:
        print("\n❌ Test failed!")
