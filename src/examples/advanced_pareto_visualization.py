#!/usr/bin/env python3
"""
Advanced Pareto Front Visualization Module
Creates visualizations similar to Figure 5 showing Pareto front approximations
with different evaluation budgets and algorithm comparisons.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Tuple, Optional
import sys
import os

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

try:
    from moecam.problems.test_functions import ZDT1, ZDT2, ZDT3
    from moecam.algorithms.moea_algorithms import NSGAII, MOEAD
    from moecam.metrics.performance_metrics import pareto_front, hypervolume
except ImportError as e:
    print(f"Import error: {e}")
    print("Make sure the MOECAM package is properly installed")

class ParetoFrontAnalyzer:
    """
    Advanced analyzer for Pareto front approximations with multiple evaluation budgets
    """

    def __init__(self):
        self.colors = {
            'pareto_true': 'black',
            'eval_200': 'red',
            'eval_500': 'orange',
            'eval_1000': 'blue',
            'eval_1800': 'green',
            'eval_3000': 'purple'
        }

        self.markers = {
            'pareto_true': 's',  # squares
            'eval_200': '.',     # dots
            'eval_500': '.',
            'eval_1000': '.',
            'eval_1800': '.',
            'eval_3000': '.'
        }

        self.sizes = {
            'pareto_true': 25,
            'eval_200': 15,
            'eval_500': 15,
            'eval_1000': 15,
            'eval_1800': 15,
            'eval_3000': 15
        }

    def generate_true_pareto_front(self, problem_name: str, n_points: int = 100) -> np.ndarray:
        """Generate the true Pareto front for benchmark problems"""

        if problem_name.upper() == 'ZDT1':
            # True Pareto front: f2 = 1 - sqrt(f1)
            f1 = np.linspace(0, 1, n_points)
            f2 = 1 - np.sqrt(f1)
            return np.column_stack([f1, f2])

        elif problem_name.upper() == 'ZDT2':
            # True Pareto front: f2 = 1 - f1^2
            f1 = np.linspace(0, 1, n_points)
            f2 = 1 - f1**2
            return np.column_stack([f1, f2])

        elif problem_name.upper() == 'ZDT3':
            # Disconnected Pareto front
            f1_segments = []
            f2_segments = []

            # Define the disconnected segments
            segments = [(0.0, 0.0830015349), (0.1822287280, 0.2577623634),
                       (0.4093136748, 0.4538821041), (0.6183967944, 0.6525117038),
                       (0.8233317983, 0.8518328654)]

            for start, end in segments:
                f1_seg = np.linspace(start, end, n_points//len(segments))
                f2_seg = 1 - np.sqrt(f1_seg) - f1_seg * np.sin(10 * np.pi * f1_seg)
                f1_segments.extend(f1_seg)
                f2_segments.extend(f2_seg)

            return np.column_stack([f1_segments, f2_segments])

        else:
            raise ValueError(f"Unknown problem: {problem_name}")

    def run_algorithm_with_budget(self, algorithm_class, problem, evaluation_budget: int,
                                run_id: int = 0) -> Tuple[np.ndarray, Dict, np.ndarray]:
        """Run MOECAM algorithm with specific evaluation budget"""

        # Calculate population size and generations based on budget
        if evaluation_budget <= 200:
            pop_size = 20
            num_generations = evaluation_budget // pop_size
        elif evaluation_budget <= 1000:
            pop_size = 50
            num_generations = evaluation_budget // pop_size
        else:
            pop_size = 100
            num_generations = evaluation_budget // pop_size

        # Ensure minimum generations
        num_generations = max(num_generations, 1)

        print(f"  Run {run_id+1}: Budget={evaluation_budget}, Pop={pop_size}, Gen={num_generations}")

        # Initialize MOECAM algorithm
        algorithm = algorithm_class(problem, pop_size=pop_size, num_generations=num_generations)

        # Store all evaluated solutions during optimization
        original_evaluate = algorithm._evaluate_population
        all_evaluated_solutions = []

        def tracking_evaluate(population):
            objectives = original_evaluate(population)
            # Store all solutions for visualization
            for i, (solution, objective) in enumerate(zip(population, objectives)):
                all_evaluated_solutions.append(objective.copy())
            return objectives

        algorithm._evaluate_population = tracking_evaluate

        # Run optimization
        pareto_solutions = algorithm.optimize()

        # Get all evaluated solutions as numpy array
        all_solutions = np.array(all_evaluated_solutions) if all_evaluated_solutions else np.array([]).reshape(0, 2)

        # Calculate metrics using actual algorithm performance tracking
        actual_evaluations = algorithm.evaluation_count if hasattr(algorithm, 'evaluation_count') else len(all_solutions)
        execution_time = algorithm.execution_time if hasattr(algorithm, 'execution_time') else 0.0

        # Calculate hypervolume
        reference_point = [2.0, 2.0]  # Conservative reference point
        hv = hypervolume(pareto_solutions, reference_point) if len(pareto_solutions) > 0 else 0.0

        metrics = {
            'hypervolume': hv,
            'n_solutions': len(pareto_solutions),
            'evaluations_used': actual_evaluations,
            'execution_time': execution_time,
            'total_evaluated': len(all_solutions)
        }

        return pareto_solutions, metrics, all_solutions

    def multi_budget_analysis(self, problem_name: str = 'ZDT1',
                            algorithm_name: str = 'NSGAII',
                            evaluation_budgets: List[int] = [200, 500, 1000, 1800],
                            n_runs: int = 5) -> Dict:
        """
        Perform multi-budget analysis similar to Figure 5
        """

        print(f"Starting multi-budget analysis for {problem_name} using {algorithm_name}")
        print(f"Evaluation budgets: {evaluation_budgets}")
        print(f"Number of runs per budget: {n_runs}")

        # Initialize problem
        if problem_name.upper() == 'ZDT1':
            problem = ZDT1(n_dim=30)
        elif problem_name.upper() == 'ZDT2':
            problem = ZDT2(n_dim=30)
        elif problem_name.upper() == 'ZDT3':
            problem = ZDT3(n_dim=30)
        else:
            raise ValueError(f"Unknown problem: {problem_name}")

        # Initialize algorithm class
        if algorithm_name.upper() == 'NSGAII':
            algorithm_class = NSGAII
        elif algorithm_name.upper() == 'MOEAD':
            algorithm_class = MOEAD
        else:
            raise ValueError(f"Unknown algorithm: {algorithm_name}")

        # Generate true Pareto front
        true_pareto = self.generate_true_pareto_front(problem_name, n_points=200)

        # Storage for results
        results = {
            'true_pareto': true_pareto,
            'budgets': {},
            'summary_stats': {}
        }

        # Run experiments for each budget
        for budget in evaluation_budgets:
            print(f"\nRunning MOECAM experiments with {budget} evaluations...")

            budget_results = {
                'pareto_fronts': [],
                'metrics': [],
                'all_solutions': []  # All evaluated solutions from all runs
            }

            for run in range(n_runs):
                try:
                    pareto_solutions, metrics, all_evaluated = self.run_algorithm_with_budget(
                        algorithm_class, problem, budget, run
                    )

                    budget_results['pareto_fronts'].append(pareto_solutions)
                    budget_results['metrics'].append(metrics)
                    # Add all evaluated solutions from this run
                    if len(all_evaluated) > 0:
                        budget_results['all_solutions'].extend(all_evaluated)

                    print(f"    Run {run+1}: {len(pareto_solutions)} Pareto solutions, "
                          f"{len(all_evaluated)} total evaluations, HV = {metrics['hypervolume']:.4f}")

                except Exception as e:
                    print(f"    Error in run {run+1}: {e}")
                    continue

            # Convert all solutions to numpy array
            if budget_results['all_solutions']:
                budget_results['all_solutions'] = np.array(budget_results['all_solutions'])
            else:
                budget_results['all_solutions'] = np.array([]).reshape(0, 2)

            # Calculate summary statistics
            if budget_results['metrics']:
                hvs = [m['hypervolume'] for m in budget_results['metrics']]
                n_sols = [m['n_solutions'] for m in budget_results['metrics']]
                exec_times = [m['execution_time'] for m in budget_results['metrics']]
                actual_evals = [m['evaluations_used'] for m in budget_results['metrics']]

                summary = {
                    'mean_hypervolume': np.mean(hvs),
                    'std_hypervolume': np.std(hvs),
                    'mean_solutions': np.mean(n_sols),
                    'std_solutions': np.std(n_sols),
                    'mean_execution_time': np.mean(exec_times),
                    'std_execution_time': np.std(exec_times),
                    'mean_evaluations': np.mean(actual_evals),
                    'total_solutions': len(budget_results['all_solutions'])
                }
            else:
                summary = {
                    'mean_hypervolume': 0.0,
                    'std_hypervolume': 0.0,
                    'mean_solutions': 0.0,
                    'std_solutions': 0.0,
                    'mean_execution_time': 0.0,
                    'std_execution_time': 0.0,
                    'mean_evaluations': 0.0,
                    'total_solutions': 0
                }

            results['budgets'][budget] = budget_results
            results['summary_stats'][budget] = summary

            print(f"  Completed: {summary['total_solutions']} total solutions evaluated, "
                  f"HV = {summary['mean_hypervolume']:.4f} ± {summary['std_hypervolume']:.4f}, "
                  f"Avg time = {summary['mean_execution_time']:.3f}s")

        return results

    def create_figure5_style_plot(self, results: Dict, problem_name: str,
                                algorithm_name: str, save_path: str = None) -> None:
        """
        Create a visualization similar to Figure 5
        """

        plt.figure(figsize=(12, 8))

        # Set style
        plt.style.use('seaborn-v0_8')

        # Plot all candidate solutions for each budget
        budget_keys = sorted(results['budgets'].keys())

        for i, budget in enumerate(budget_keys):
            budget_data = results['budgets'][budget]
            solutions = budget_data['all_solutions']

            if len(solutions) > 0:
                color = list(self.colors.values())[i + 1]  # Skip 'pareto_true' color
                marker = self.markers[f'eval_{budget}'] if f'eval_{budget}' in self.markers else '.'
                size = self.sizes[f'eval_{budget}'] if f'eval_{budget}' in self.sizes else 15

                plt.scatter(solutions[:, 0], solutions[:, 1],
                          c=color, marker=marker, s=size, alpha=0.6,
                          label=f'using {budget} evaluations')

        # Plot true Pareto front (benchmarking Pareto)
        true_pareto = results['true_pareto']
        plt.plot(true_pareto[:, 0], true_pareto[:, 1],
                color=self.colors['pareto_true'],
                marker=self.markers['pareto_true'],
                markersize=4, linewidth=2,
                label='Benchmarking Pareto', markerfacecolor='none')

        # Formatting
        if problem_name.upper() == 'ZDT1':
            plt.xlabel('Objective 1 (f₁)', fontsize=14)
            plt.ylabel('Objective 2 (f₂)', fontsize=14)
            plt.title(f'MOECAM Analysis: {problem_name} using {algorithm_name}\n'
                     f'Population size varies, Multiple runs per budget', fontsize=16)
        else:
            plt.xlabel('Primary Energy [kWh/m2]', fontsize=14)
            plt.ylabel('LCC [Eur/m2]', fontsize=14)
            plt.title(f'{problem_name}: Population size = 64, FIT = 1, ε = 6%', fontsize=16)

        # Legend
        plt.legend(loc='upper right', fontsize=12)

        # Grid
        plt.grid(True, alpha=0.3)

        # Adjust layout
        plt.tight_layout()

        # Save or show
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved to: {save_path}")

        plt.show()

    def create_performance_summary_plot(self, results: Dict, problem_name: str,
                                      algorithm_name: str, save_path: str = None) -> None:
        """
        Create a summary plot showing hypervolume vs evaluation budget
        """

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

        budgets = sorted(results['summary_stats'].keys())
        hvs_mean = [results['summary_stats'][b]['mean_hypervolume'] for b in budgets]
        hvs_std = [results['summary_stats'][b]['std_hypervolume'] for b in budgets]
        n_sols_mean = [results['summary_stats'][b]['mean_solutions'] for b in budgets]
        n_sols_std = [results['summary_stats'][b]['std_solutions'] for b in budgets]

        # Hypervolume plot
        ax1.errorbar(budgets, hvs_mean, yerr=hvs_std, marker='o', capsize=5, linewidth=2)
        ax1.set_xlabel('Evaluation Budget', fontsize=12)
        ax1.set_ylabel('Hypervolume', fontsize=12)
        ax1.set_title(f'Hypervolume vs Budget\n{problem_name} - {algorithm_name}', fontsize=14)
        ax1.grid(True, alpha=0.3)

        # Number of solutions plot
        ax2.errorbar(budgets, n_sols_mean, yerr=n_sols_std, marker='s', capsize=5,
                    linewidth=2, color='orange')
        ax2.set_xlabel('Evaluation Budget', fontsize=12)
        ax2.set_ylabel('Number of Pareto Solutions', fontsize=12)
        ax2.set_title(f'Solution Count vs Budget\n{problem_name} - {algorithm_name}', fontsize=14)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Performance summary saved to: {save_path}")

        plt.show()

def main():
    """
    Main function to generate Figure 5 style visualization
    """

    # Initialize analyzer
    analyzer = ParetoFrontAnalyzer()

    # Configuration
    problem_name = 'ZDT1'
    algorithm_name = 'NSGAII'
    evaluation_budgets = [200, 500, 1000, 1800]
    n_runs = 5

    print("=" * 60)
    print("MOECAM: Advanced Pareto Front Analysis")
    print("=" * 60)

    try:
        # Run multi-budget analysis
        results = analyzer.multi_budget_analysis(
            problem_name=problem_name,
            algorithm_name=algorithm_name,
            evaluation_budgets=evaluation_budgets,
            n_runs=n_runs
        )

        print(f"\nAnalysis completed successfully!")
        print(f"Results summary:")
        for budget in sorted(results['summary_stats'].keys()):
            stats = results['summary_stats'][budget]
            print(f"  Budget {budget}: {stats['total_solutions']} solutions, "
                  f"HV = {stats['mean_hypervolume']:.4f} ± {stats['std_hypervolume']:.4f}")

        # Create visualizations
        print(f"\nGenerating visualizations...")

        # Figure 5 style plot
        analyzer.create_figure5_style_plot(
            results, problem_name, algorithm_name,
            save_path='pareto_front_budget_analysis.png'
        )

        # Performance summary plot
        analyzer.create_performance_summary_plot(
            results, problem_name, algorithm_name,
            save_path='performance_vs_budget.png'
        )

        print(f"\nVisualization complete! Check the generated PNG files.")

    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
