"""
Advanced Visualization Module for Multi-Objective Optimization Analysis
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches

# Set publication-quality plotting style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

class MOOVisualizer:
    """Comprehensive visualization tools for multi-objective optimization"""
    
    def __init__(self, figsize=(12, 8), dpi=300):
        self.figsize = figsize
        self.dpi = dpi
        self.colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                      '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    
    def plot_pareto_comparison(self, results_dict, problem_name="Test Problem", 
                              save_path=None, show_all_points=True):
        """
        Plot Pareto fronts comparison similar to Figure 5 in the reference
        
        Args:
            results_dict: Dict with algorithm names as keys and results as values
            problem_name: Name of the test problem
            save_path: Path to save the figure
            show_all_points: Whether to show all evaluated points
        """
        fig, ax = plt.subplots(figsize=self.figsize, dpi=self.dpi)
        
        # Plot all evaluated points if available (like red and green dots in Fig 5)
        for i, (alg_name, result) in enumerate(results_dict.items()):
            color = self.colors[i % len(self.colors)]
            
            # Plot all evaluated solutions if available
            if show_all_points and hasattr(result, 'all_evaluated_solutions'):
                all_points = result.all_evaluated_solutions
                ax.scatter(all_points[:, 0], all_points[:, 1], 
                          s=8, alpha=0.3, color=color, 
                          label=f'{alg_name} (all evaluations)')
            
            # Plot Pareto front
            pareto_front = result['pareto_front'] if isinstance(result, dict) else result
            if len(pareto_front) > 0:
                # Sort for better line plotting
                sorted_pf = pareto_front[pareto_front[:, 0].argsort()]
                ax.plot(sorted_pf[:, 0], sorted_pf[:, 1], 
                       'o-', color=color, linewidth=2, markersize=6,
                       label=f'{alg_name} Pareto Front', alpha=0.8)
        
        # Add theoretical Pareto front if known (black squares like in Fig 5)
        if problem_name.startswith('ZDT1'):
            f1_theory = np.linspace(0, 1, 100)
            f2_theory = 1 - np.sqrt(f1_theory)
            ax.plot(f1_theory, f2_theory, 'ks-', markersize=3, linewidth=1,
                   label='True Pareto Front', alpha=0.7)
        
        ax.set_xlabel('Objective 1 (f₁)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Objective 2 (f₂)', fontsize=12, fontweight='bold')
        ax.set_title(f'{problem_name} - Pareto Front Comparison', 
                    fontsize=14, fontweight='bold')
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=self.dpi, bbox_inches='tight')
        
        return fig, ax
    
    def plot_disconnected_pareto(self, pareto_front, problem_name="R2", save_path=None):
        """
        Plot disconnected Pareto fronts like Figures 6-8 in the reference
        """
        fig, ax = plt.subplots(figsize=(8, 6), dpi=self.dpi)
        
        if len(pareto_front) > 0:
            # For disconnected fronts, don't connect with lines
            ax.scatter(pareto_front[:, 0], pareto_front[:, 1], 
                      s=50, color='darkblue', alpha=0.7, edgecolors='black')
            
            # Add theoretical curves for specific problems
            if problem_name == "R2":
                # R2 has multiple disconnected segments
                x_segments = [
                    np.linspace(-1, -0.5, 50),
                    np.linspace(-0.3, 0.2, 50),
                    np.linspace(0.4, 0.7, 50)
                ]
                for x_seg in x_segments:
                    y_seg = self._r2_theoretical_front(x_seg)
                    ax.plot(x_seg, y_seg, 'k-', linewidth=2, alpha=0.8)
            
            elif problem_name == "R3":
                # R3 theoretical front
                x_theory = np.linspace(0.3, 1, 100)
                y_theory = 1 - x_theory
                ax.plot(x_theory, y_theory, 'k-', linewidth=2, alpha=0.8)
            
            elif problem_name == "R4":
                # R4 theoretical front (curved)
                x_theory = np.linspace(-1, 1, 100)
                y_theory = 1 - x_theory**2
                ax.plot(x_theory, y_theory, 'k-', linewidth=2, alpha=0.8)
        
        ax.set_xlabel('f₁', fontsize=12, fontweight='bold')
        ax.set_ylabel('f₂', fontsize=12, fontweight='bold')
        ax.set_title(f'{problem_name} Pareto-optimal front and feasible region', 
                    fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=self.dpi, bbox_inches='tight')
        
        return fig, ax
    
    def _r2_theoretical_front(self, x):
        """Theoretical R2 Pareto front calculation"""
        # Simplified R2 function for demonstration
        return 4 - x**2 - 2*np.abs(x)
    
    def plot_convergence_analysis(self, algorithms_results, save_path=None):
        """Plot convergence analysis showing hypervolume evolution"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10), dpi=self.dpi)
        
        for i, (alg_name, result) in enumerate(algorithms_results.items()):
            color = self.colors[i % len(self.colors)]
            
            if 'convergence_history' in result:
                history = result['convergence_history']
                generations = [h['generation'] for h in history]
                hypervolumes = [h['hypervolume'] for h in history]
                pareto_sizes = [h['pareto_size'] for h in history]
                evaluations = [h['evaluations'] for h in history]
                
                # Hypervolume vs Generation
                ax1.plot(generations, hypervolumes, 'o-', color=color, 
                        label=alg_name, alpha=0.8)
                
                # Pareto set size vs Generation
                ax2.plot(generations, pareto_sizes, 's-', color=color, 
                        label=alg_name, alpha=0.8)
                
                # Hypervolume vs Evaluations
                ax3.plot(evaluations, hypervolumes, '^-', color=color, 
                        label=alg_name, alpha=0.8)
                
                # Convergence rate
                if len(hypervolumes) > 1:
                    conv_rate = np.diff(hypervolumes)
                    ax4.plot(generations[1:], conv_rate, 'v-', color=color, 
                            label=alg_name, alpha=0.8)
        
        ax1.set_xlabel('Generation')
        ax1.set_ylabel('Hypervolume')
        ax1.set_title('Hypervolume Convergence')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        ax2.set_xlabel('Generation')
        ax2.set_ylabel('Pareto Set Size')
        ax2.set_title('Pareto Set Size Evolution')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        ax3.set_xlabel('Function Evaluations')
        ax3.set_ylabel('Hypervolume')
        ax3.set_title('Efficiency Analysis')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        ax4.set_xlabel('Generation')
        ax4.set_ylabel('HV Improvement Rate')
        ax4.set_title('Convergence Rate')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=self.dpi, bbox_inches='tight')
        
        return fig, (ax1, ax2, ax3, ax4)
    
    def plot_performance_comparison(self, algorithms_results, save_path=None):
        """Create comprehensive performance comparison charts"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10), dpi=self.dpi)
        
        alg_names = list(algorithms_results.keys())
        n_algs = len(alg_names)
        
        # Extract metrics
        hypervolumes = []
        execution_times = []
        evaluations = []
        pareto_sizes = []
        
        for alg_name in alg_names:
            result = algorithms_results[alg_name]
            
            # Final hypervolume
            if 'convergence_history' in result and result['convergence_history']:
                hypervolumes.append(result['convergence_history'][-1]['hypervolume'])
                pareto_sizes.append(result['convergence_history'][-1]['pareto_size'])
            else:
                hypervolumes.append(0)
                pareto_sizes.append(0)
            
            # Performance metrics
            perf_metrics = result.get('performance_metrics', {})
            execution_times.append(perf_metrics.get('execution_time', 0))
            evaluations.append(perf_metrics.get('total_evaluations', 0))
        
        # Bar charts
        x_pos = np.arange(n_algs)
        
        # Hypervolume comparison
        bars1 = ax1.bar(x_pos, hypervolumes, color=self.colors[:n_algs], alpha=0.7)
        ax1.set_xlabel('Algorithm')
        ax1.set_ylabel('Final Hypervolume')
        ax1.set_title('Hypervolume Comparison')
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels(alg_names, rotation=45)
        
        # Execution time comparison
        bars2 = ax2.bar(x_pos, execution_times, color=self.colors[:n_algs], alpha=0.7)
        ax2.set_xlabel('Algorithm')
        ax2.set_ylabel('Execution Time (seconds)')
        ax2.set_title('Runtime Comparison')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(alg_names, rotation=45)
        
        # Function evaluations
        bars3 = ax3.bar(x_pos, evaluations, color=self.colors[:n_algs], alpha=0.7)
        ax3.set_xlabel('Algorithm')
        ax3.set_ylabel('Function Evaluations')
        ax3.set_title('Computational Effort')
        ax3.set_xticks(x_pos)
        ax3.set_xticklabels(alg_names, rotation=45)
        
        # Efficiency scatter plot (HV vs Time)
        ax4.scatter(execution_times, hypervolumes, s=100, 
                   c=self.colors[:n_algs], alpha=0.7)
        for i, alg_name in enumerate(alg_names):
            ax4.annotate(alg_name, (execution_times[i], hypervolumes[i]),
                        xytext=(5, 5), textcoords='offset points')
        ax4.set_xlabel('Execution Time (seconds)')
        ax4.set_ylabel('Final Hypervolume')
        ax4.set_title('Efficiency Analysis')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=self.dpi, bbox_inches='tight')
        
        return fig, ((ax1, ax2), (ax3, ax4))
    
    def plot_3d_pareto_front(self, pareto_front, algorithm_name="Algorithm", save_path=None):
        """Plot 3D Pareto front for problems with 3 objectives"""
        if pareto_front.shape[1] != 3:
            raise ValueError("3D plotting requires exactly 3 objectives")
        
        fig = plt.figure(figsize=self.figsize, dpi=self.dpi)
        ax = fig.add_subplot(111, projection='3d')
        
        # Create 3D scatter plot
        scatter = ax.scatter(pareto_front[:, 0], pareto_front[:, 1], pareto_front[:, 2],
                           c=pareto_front[:, 2], cmap='viridis', s=50, alpha=0.7)
        
        ax.set_xlabel('Objective 1 (f₁)', fontweight='bold')
        ax.set_ylabel('Objective 2 (f₂)', fontweight='bold')
        ax.set_zlabel('Objective 3 (f₃)', fontweight='bold')
        ax.set_title(f'3D Pareto Front - {algorithm_name}', fontweight='bold')
        
        # Add colorbar
        plt.colorbar(scatter, ax=ax, shrink=0.5, aspect=5)
        
        if save_path:
            plt.savefig(save_path, dpi=self.dpi, bbox_inches='tight')
        
        return fig, ax

def create_comprehensive_analysis(algorithms_results, problem_name="Test Problem", 
                                output_dir="./analysis_output/"):
    """
    Create a comprehensive analysis report with all visualizations
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    visualizer = MOOVisualizer()
    
    # 1. Pareto front comparison
    fig1, _ = visualizer.plot_pareto_comparison(
        algorithms_results, problem_name, 
        save_path=os.path.join(output_dir, f"{problem_name}_pareto_comparison.png")
    )
    
    # 2. Convergence analysis
    fig2, _ = visualizer.plot_convergence_analysis(
        algorithms_results,
        save_path=os.path.join(output_dir, f"{problem_name}_convergence_analysis.png")
    )
    
    # 3. Performance comparison
    fig3, _ = visualizer.plot_performance_comparison(
        algorithms_results,
        save_path=os.path.join(output_dir, f"{problem_name}_performance_comparison.png")
    )
    
    plt.show()
    
    return fig1, fig2, fig3
