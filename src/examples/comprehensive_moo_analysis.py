#!/usr/bin/env python3
"""
Comprehensive Multi-Objective Optimization Analysis Script

This script demonstrates:
1. Theoretical differences between MOO and SOO
2. Enhanced NSGA-II with C++ integration
3. Advanced visualization and performance analysis
4. Comprehensive benchmarking
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import os

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from src.moecam.problems.test_functions import ZDT1, ZDT2, ZDT3, DTLZ1, DTLZ2
from src.moecam.algorithms.moea_algorithms import NSGAII, MOEAD
from src.moecam.metrics.performance_metrics import hypervolume, pareto_front
from src.moecam.visualization import MOOVisualizer, create_comprehensive_analysis

class MOOAnalysis:
    """Comprehensive analysis of Multi-Objective Optimization"""

    def __init__(self):
        self.visualizer = MOOVisualizer()

    def demonstrate_moo_vs_soo_differences(self):
        """
        Demonstrate fundamental differences between MOO and SOO
        """
        print("=" * 80)
        print("MULTI-OBJECTIVE OPTIMIZATION (MOO) vs SINGLE-OBJECTIVE OPTIMIZATION (SOO)")
        print("=" * 80)

        print("\n1. FUNDAMENTAL DIFFERENCES:")
        print("-" * 40)
        print("SOO: Seeks ONE optimal solution")
        print("MOO: Seeks SET of optimal trade-off solutions (Pareto Set)")
        print()
        print("SOO: Total ordering of solutions (better/worse)")
        print("MOO: Partial ordering via dominance relations")
        print()
        print("SOO: Convergence to single point")
        print("MOO: Convergence to Pareto front + diversity maintenance")

        print("\n2. ALGORITHMIC IMPLICATIONS:")
        print("-" * 40)
        print("SOO algorithms (gradient-based, simulated annealing) cannot work for MOO because:")
        print("• No single 'best' solution exists")
        print("• Need to maintain population diversity")
        print("• Require dominance-based selection")
        print("• Must balance convergence AND diversity")

        print("\n3. ALGORITHM DESIGN REQUIREMENTS:")
        print("-" * 40)
        print("MOO algorithms must implement:")
        print("• Non-dominated sorting (rank assignment)")
        print("• Diversity preservation mechanisms (crowding distance)")
        print("• Multi-objective selection operators")
        print("• Population-based search")

        print("\n4. PERFORMANCE METRICS:")
        print("-" * 40)
        print("SOO: Function value, convergence rate")
        print("MOO: Hypervolume, IGD, spread, convergence metrics")

        # Demonstrate with actual example
        self._demonstrate_pareto_concept()

    def _demonstrate_pareto_concept(self):
        """Demonstrate Pareto dominance concept with concrete example"""
        print("\n5. PARETO DOMINANCE DEMONSTRATION:")
        print("-" * 40)

        # Create sample solutions for bi-objective minimization
        solutions = np.array([
            [1.0, 5.0],  # A
            [2.0, 4.0],  # B
            [3.0, 3.0],  # C
            [4.0, 2.0],  # D
            [5.0, 1.0],  # E
            [2.5, 4.5],  # F (dominated)
            [1.5, 3.5],  # G (non-dominated)
        ])

        labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

        # Find Pareto front
        pf = pareto_front(solutions)

        print(f"Solutions: {dict(zip(labels, solutions.tolist()))}")
        print(f"Pareto optimal solutions: {len(pf)} out of {len(solutions)}")

        # Visualize
        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot all solutions
        for i, (sol, label) in enumerate(zip(solutions, labels)):
            is_pareto = any(np.array_equal(sol, pf_sol) for pf_sol in pf)
            color = 'red' if is_pareto else 'blue'
            marker = 'o' if is_pareto else 's'
            ax.scatter(sol[0], sol[1], color=color, s=100, marker=marker)
            ax.annotate(label, (sol[0], sol[1]), xytext=(5, 5),
                       textcoords='offset points', fontweight='bold')

        # Draw Pareto front
        pf_sorted = pf[pf[:, 0].argsort()]
        ax.plot(pf_sorted[:, 0], pf_sorted[:, 1], 'r--', linewidth=2, alpha=0.7)

        ax.set_xlabel('Objective 1 (minimize)', fontweight='bold')
        ax.set_ylabel('Objective 2 (minimize)', fontweight='bold')
        ax.set_title('Pareto Dominance Concept', fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(['Pareto Front', 'Pareto Optimal', 'Dominated'], loc='upper right')

        plt.tight_layout()
        plt.savefig('pareto_concept_demonstration.png', dpi=300, bbox_inches='tight')
        plt.show()

        print("→ Red circles: Pareto optimal (non-dominated)")
        print("→ Blue squares: Dominated solutions")
        print("→ Red dashed line: Pareto front")

    def run_comprehensive_benchmark(self):
        """Run comprehensive benchmarking study"""
        print("\n" + "=" * 80)
        print("COMPREHENSIVE BENCHMARKING STUDY")
        print("=" * 80)

        # Test problems
        problems = {
            'ZDT1': ZDT1(n_dim=30),
            'ZDT2': ZDT2(n_dim=30),
            'ZDT3': ZDT3(n_dim=30),
            'DTLZ1': DTLZ1(n_dim=7, n_obj=3),
            'DTLZ2': DTLZ2(n_dim=12, n_obj=3)
        }

        # Algorithm configurations
        algorithms = {
            'NSGA-II (Python)': lambda prob: NSGAII(prob, pop_size=100, num_generations=250, use_cpp=False),
            'NSGA-II (C++)': lambda prob: NSGAII(prob, pop_size=100, num_generations=250, use_cpp=True),
            'MOEA/D': lambda prob: MOEAD(prob, pop_size=100, num_generations=250)
        }

        all_results = {}

        for prob_name, problem in problems.items():
            print(f"\n--- Testing {prob_name} ---")

            if problem.n_obj > 2:
                print(f"Skipping visualization for {prob_name} (>2 objectives)")
                continue

            prob_results = {}

            for alg_name, alg_factory in algorithms.items():
                if 'C++' in alg_name and not hasattr(NSGAII, 'cpp_optimizer'):
                    print(f"Skipping {alg_name} (C++ not available)")
                    continue

                print(f"  Running {alg_name}...")

                try:
                    # Create and run algorithm
                    algorithm = alg_factory(problem)

                    start_time = time.time()
                    pareto_front_result = algorithm.optimize()
                    end_time = time.time()

                    # Calculate metrics
                    reference_point = np.array([1.1, 1.1]) if problem.n_obj == 2 else np.array([1.1] * problem.n_obj)
                    hv = hypervolume(pareto_front_result, reference_point)

                    # Store results
                    result = {
                        'pareto_front': pareto_front_result,
                        'hypervolume': hv,
                        'execution_time': end_time - start_time,
                        'algorithm': algorithm
                    }

                    # Add performance metrics if available
                    if hasattr(algorithm, 'get_performance_metrics'):
                        result['performance_metrics'] = algorithm.get_performance_metrics()
                        result['convergence_history'] = algorithm.convergence_history

                    prob_results[alg_name] = result

                    print(f"    ✓ Solutions: {len(pareto_front_result)}, HV: {hv:.6f}, Time: {end_time-start_time:.2f}s")

                except Exception as e:
                    print(f"    ✗ Error: {e}")
                    continue

            all_results[prob_name] = prob_results

            # Create visualizations for this problem
            if prob_results:
                self._create_problem_analysis(prob_name, prob_results)

        return all_results

    def _create_problem_analysis(self, problem_name, results):
        """Create comprehensive analysis for a specific problem"""
        print(f"  Creating analysis for {problem_name}...")

        # 1. Pareto front comparison
        fig1, ax1 = self.visualizer.plot_pareto_comparison(
            results, problem_name,
            save_path=f"{problem_name}_pareto_comparison.png",
            show_all_points=True
        )

        # 2. Performance comparison
        fig2, axes2 = self.visualizer.plot_performance_comparison(
            results,
            save_path=f"{problem_name}_performance_comparison.png"
        )

        # 3. Convergence analysis (if data available)
        has_convergence_data = any('convergence_history' in result for result in results.values())
        if has_convergence_data:
            fig3, axes3 = self.visualizer.plot_convergence_analysis(
                results,
                save_path=f"{problem_name}_convergence_analysis.png"
            )

        plt.show()

    def demonstrate_efficiency_metrics(self, results):
        """Demonstrate how to measure MOO algorithm efficiency"""
        print("\n" + "=" * 80)
        print("EFFICIENCY MEASUREMENT IN MULTI-OBJECTIVE OPTIMIZATION")
        print("=" * 80)

        print("\n1. SPEED METRICS:")
        print("-" * 20)
        print("• CPU Time: Wall-clock execution time")
        print("• Function Evaluations: Number of objective function calls")
        print("• Convergence Rate: Hypervolume improvement per generation")

        print("\n2. QUALITY METRICS:")
        print("-" * 20)
        print("• Hypervolume: Volume dominated by Pareto set")
        print("• IGD (Inverted Generational Distance): Distance to true front")
        print("• Spread: Distribution of solutions along front")
        print("• Convergence: Proximity to true Pareto front")

        print("\n3. EFFICIENCY ANALYSIS:")
        print("-" * 20)

        for prob_name, prob_results in results.items():
            if not prob_results:
                continue

            print(f"\n{prob_name}:")

            for alg_name, result in prob_results.items():
                perf = result.get('performance_metrics', {})

                print(f"  {alg_name}:")
                print(f"    • Evaluations: {perf.get('total_evaluations', 'N/A')}")
                print(f"    • Time: {result.get('execution_time', 'N/A'):.3f}s")
                print(f"    • Final HV: {result.get('hypervolume', 'N/A'):.6f}")
                print(f"    • Solutions: {len(result.get('pareto_front', []))}")

                # Calculate efficiency ratio
                if 'execution_time' in result and result['execution_time'] > 0:
                    hv_per_second = result.get('hypervolume', 0) / result['execution_time']
                    print(f"    • HV/second: {hv_per_second:.6f}")

    def create_publication_plots(self, results):
        """Create publication-quality plots similar to the reference figures"""
        print("\n" + "=" * 80)
        print("CREATING PUBLICATION-QUALITY VISUALIZATIONS")
        print("=" * 80)

        # Create plots similar to Figures 5-8 in the reference
        for prob_name, prob_results in results.items():
            if not prob_results or len(prob_results) == 0:
                continue

            print(f"\nCreating plots for {prob_name}...")

            # Figure 5 style: Comprehensive Pareto front comparison
            fig, ax = plt.subplots(figsize=(12, 8), dpi=300)

            colors = ['green', 'red', 'blue', 'orange', 'purple']
            markers = ['o', 's', '^', 'v', 'D']

            # Plot with different evaluation counts (like 200 vs 1800 in Fig 5)
            for i, (alg_name, result) in enumerate(prob_results.items()):
                color = colors[i % len(colors)]
                marker = markers[i % len(markers)]

                pareto_front_data = result['pareto_front']

                # Plot Pareto front points
                if len(pareto_front_data) > 0:
                    ax.scatter(pareto_front_data[:, 0], pareto_front_data[:, 1],
                             color=color, marker=marker, s=50, alpha=0.8,
                             label=f'{alg_name} ({len(pareto_front_data)} solutions)')

                # Plot all evaluated points if available
                if hasattr(result.get('algorithm'), 'all_evaluated_solutions'):
                    all_points = result['algorithm'].all_evaluated_solutions
                    ax.scatter(all_points[:, 0], all_points[:, 1],
                             color=color, marker='.', s=10, alpha=0.3)

            # Add true Pareto front for ZDT problems
            if prob_name.startswith('ZDT1'):
                f1_true = np.linspace(0, 1, 100)
                f2_true = 1 - np.sqrt(f1_true)
                ax.plot(f1_true, f2_true, 'k-', linewidth=2,
                       label='True Pareto Front', alpha=0.8)

            ax.set_xlabel('Primary Energy [kWh/m2]', fontsize=12, fontweight='bold')
            ax.set_ylabel('LCC [EUR/m2]', fontsize=12, fontweight='bold')
            ax.set_title(f'{prob_name}: Benchmarking Pareto and Candidate Solutions',
                        fontsize=14, fontweight='bold')
            ax.legend(loc='best')
            ax.grid(True, alpha=0.3)

            plt.tight_layout()
            plt.savefig(f'{prob_name}_publication_style.png', dpi=300, bbox_inches='tight')

        plt.show()

def main():
    """Main analysis function"""
    print("MOECAM: Comprehensive Multi-Objective Optimization Analysis")
    print("=" * 60)

    # Initialize analysis
    analysis = MOOAnalysis()

    # 1. Demonstrate theoretical differences
    analysis.demonstrate_moo_vs_soo_differences()

    # 2. Run comprehensive benchmarks
    results = analysis.run_comprehensive_benchmark()

    # 3. Demonstrate efficiency metrics
    analysis.demonstrate_efficiency_metrics(results)

    # 4. Create publication plots
    analysis.create_publication_plots(results)

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE!")
    print("Generated files:")
    print("• pareto_concept_demonstration.png")
    print("• [Problem]_pareto_comparison.png")
    print("• [Problem]_performance_comparison.png")
    print("• [Problem]_convergence_analysis.png")
    print("• [Problem]_publication_style.png")
    print("=" * 60)

if __name__ == "__main__":
    main()
