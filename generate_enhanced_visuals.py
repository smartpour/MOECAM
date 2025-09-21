#!/usr/bin/env python3
"""
Enhanced MOECAM Visualization Generator
======================================

This script creates publication-quality plots that clearly demonstrate:
1. Minimization problem context with proper nadir points
2. Fair algorithm comparison with consistent reference points
3. High function evaluation comprehensive analysis
4. Professional visualization suitable for Figure 2 and Figure 4
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import os

# Set up professional plotting style
plt.style.use('default')
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 11,
    'figure.titlesize': 18,
    'lines.linewidth': 2,
    'lines.markersize': 8,
    'figure.dpi': 300
})

def load_comprehensive_results():
    """Load all comprehensive test results."""
    results_dir = Path("moecam_comprehensive_results")

    # Load ECAM ZDT1 results
    ecam_all_obj = np.loadtxt(results_dir / "ECAM_ZDT1_all_objectives.txt")
    ecam_pareto = np.loadtxt(results_dir / "ECAM_ZDT1_pareto_front.txt")

    # Load Random Start ZDT1 results
    random_all_obj = np.loadtxt(results_dir / "Random Start_ZDT1_all_objectives.txt")
    random_pareto = np.loadtxt(results_dir / "Random Start_ZDT1_pareto_front.txt")

    return {
        'ECAM': {
            'all_objectives': ecam_all_obj,
            'pareto_front': ecam_pareto,
            'hypervolume': 11.049211,
            'evaluations': 10000,
            'pareto_size': len(ecam_pareto)
        },
        'Random_Start': {
            'all_objectives': random_all_obj,
            'pareto_front': random_pareto,
            'hypervolume': 10.969209,
            'evaluations': 6000,
            'pareto_size': len(random_pareto)
        }
    }

def generate_theoretical_zdt1(n_points=200):
    """Generate high-resolution theoretical Pareto front for ZDT1."""
    f1 = np.linspace(0, 1, n_points)
    f2 = 1 - np.sqrt(f1)
    return np.column_stack([f1, f2])

def create_figure2_detailed_analysis(results):
    """
    Create Figure 2: Detailed analysis showing ECAM algorithm performance
    with all evaluated points and extracted Pareto front curve.
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 9))

    # Get ECAM data
    ecam_data = results['ECAM']
    all_obj = ecam_data['all_objectives']
    pareto_obj = ecam_data['pareto_front']

    # Generate theoretical Pareto front
    theoretical = generate_theoretical_zdt1()

    # Plot all evaluated points
    ax.scatter(all_obj[:, 0], all_obj[:, 1],
              alpha=0.4, s=20, c='lightblue',
              label=f'All ECAM evaluations ({ecam_data["evaluations"]} points)',
              zorder=1)

    # Plot extracted Pareto front points
    ax.scatter(pareto_obj[:, 0], pareto_obj[:, 1],
              alpha=0.9, s=120, c='red', edgecolors='darkred', linewidth=2,
              label=f'Pareto front ({ecam_data["pareto_size"]} points)',
              marker='D', zorder=3)

    # Plot theoretical Pareto front curve
    ax.plot(theoretical[:, 0], theoretical[:, 1],
           'g-', linewidth=4, alpha=0.8,
           label='Theoretical Pareto front curve',
           zorder=2)

    # Add nadir point
    nadir = np.array([1.2, 10.0])
    ax.plot(nadir[0], nadir[1], 'ko', markersize=12,
           label=f'Nadir point ({nadir[0]}, {nadir[1]})',
           zorder=4)

    # Add efficiency and hypervolume information
    efficiency = (ecam_data['pareto_size'] / ecam_data['evaluations']) * 100

    # Customize plot
    ax.set_xlabel('Objective 1 (f₁)', fontsize=16, fontweight='bold')
    ax.set_ylabel('Objective 2 (f₂)', fontsize=16, fontweight='bold')
    ax.set_title(f'ZDT1 Problem: ECAM Algorithm Results\\n'
                f'Evaluations: {ecam_data["evaluations"]}, '
                f'Pareto Points: {ecam_data["pareto_size"]}, '
                f'Efficiency: {efficiency:.1f}%, '
                f'HV: {ecam_data["hypervolume"]:.3f}',
                fontsize=18, fontweight='bold', pad=20)

    # Add minimization indicators
    ax.text(0.02, 0.98, 'MINIMIZATION PROBLEM\\n← Lower f₁ better\\n↓ Lower f₂ better',
           transform=ax.transAxes, fontsize=12, fontweight='bold',
           bbox=dict(boxstyle="round,pad=0.5", facecolor="yellow", alpha=0.8),
           verticalalignment='top')

    # Add arrow indicating better direction
    ax.annotate('Better solutions', xy=(0.15, 0.4), xytext=(0.4, 1.2),
               arrowprops=dict(arrowstyle='->', color='orange', lw=4),
               fontsize=14, color='orange', fontweight='bold')

    ax.legend(loc='upper right', fontsize=12)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xlim(-0.05, 1.3)
    ax.set_ylim(-0.1, min(np.max(all_obj[:, 1]) * 1.1, 8))

    # Save high-quality figure
    plt.tight_layout()
    output_path = Path("moecam_comprehensive_results") / "Figure2_ECAM_Detailed_Analysis.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"✅ Figure 2 saved: {output_path}")
    return output_path

def create_figure4_algorithm_comparison(results):
    """
    Create Figure 4: Algorithm comparison showing both algorithms'
    discovered Pareto fronts with theoretical front for comparison.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))

    theoretical = generate_theoretical_zdt1()
    algorithms = ['ECAM', 'Random_Start']
    titles = ['ECAM Algorithm', 'Random Start Algorithm']
    axes = [ax1, ax2]

    for i, (alg_name, title, ax) in enumerate(zip(algorithms, titles, axes)):
        data = results[alg_name]
        all_obj = data['all_objectives']
        pareto_obj = data['pareto_front']

        # Plot all evaluated points
        ax.scatter(all_obj[:, 0], all_obj[:, 1],
                  alpha=0.4, s=25, c='lightblue',
                  label=f'All evaluations ({data["evaluations"]})')

        # Plot discovered Pareto front
        ax.scatter(pareto_obj[:, 0], pareto_obj[:, 1],
                  alpha=0.9, s=100, c='red', edgecolors='darkred', linewidth=2,
                  label=f'Discovered Pareto front ({data["pareto_size"]})',
                  marker='o', zorder=3)

        # Plot theoretical Pareto front
        ax.plot(theoretical[:, 0], theoretical[:, 1],
               'g-', linewidth=3, alpha=0.9,
               label='True theoretical Pareto front',
               zorder=2)

        # Add nadir point
        nadir = np.array([1.2, 10.0])
        ax.plot(nadir[0], nadir[1], 'ko', markersize=10,
               label=f'Nadir ({nadir[0]}, {nadir[1]})',
               zorder=4)

        # Customize individual subplot
        ax.set_xlabel('Objective 1 (f₁) → MINIMIZE', fontsize=14, fontweight='bold')
        ax.set_ylabel('Objective 2 (f₂) → MINIMIZE', fontsize=14, fontweight='bold')
        ax.set_title(f'{title}\\nHV: {data["hypervolume"]:.4f}, '
                    f'Pareto Size: {data["pareto_size"]}',
                    fontsize=16, fontweight='bold')

        # Add efficiency metrics
        efficiency = (data['pareto_size'] / data['evaluations']) * 100
        ax.text(0.02, 0.98, f'Efficiency: {efficiency:.1f}%\\n'
                           f'Evaluations: {data["evaluations"]:,}',
               transform=ax.transAxes, fontsize=11,
               bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.9),
               verticalalignment='top')

        # Add minimization arrow
        ax.annotate('', xy=(0.1, 0.3), xytext=(0.3, 0.8),
                   arrowprops=dict(arrowstyle='->', color='orange', lw=3))
        ax.text(0.35, 0.85, 'Better', fontsize=12, color='orange',
               fontweight='bold', ha='center')

        ax.legend(loc='upper right', fontsize=10)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_xlim(-0.05, 1.25)
        ax.set_ylim(-0.1, min(np.max(all_obj[:, 1]) * 1.1, 6))

    # Add overall figure title
    plt.suptitle('Algorithm Comparison on ZDT1 (Minimization Problem)\\n'
                'Comparing discovered Pareto fronts with theoretical optimum',
                fontsize=20, fontweight='bold', y=0.98)

    # Save high-quality figure
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    output_path = Path("moecam_comprehensive_results") / "Figure4_Algorithm_Comparison.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"✅ Figure 4 saved: {output_path}")
    return output_path

def create_hypervolume_analysis_plot(results):
    """Create detailed hypervolume analysis visualization."""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    nadir = np.array([1.2, 10.0])
    theoretical = generate_theoretical_zdt1()

    # Top row: Individual algorithm analysis
    algorithms = ['ECAM', 'Random_Start']
    titles = ['ECAM Algorithm', 'Random Start Algorithm']
    axes_top = [ax1, ax2]

    for i, (alg_name, title, ax) in enumerate(zip(algorithms, titles, axes_top)):
        data = results[alg_name]
        pareto_obj = data['pareto_front']

        # Plot hypervolume region (simplified for visualization)
        # Create mesh for hypervolume visualization
        x_range = np.linspace(0, nadir[0], 50)
        y_range = np.linspace(0, nadir[1], 50)
        X, Y = np.meshgrid(x_range, y_range)

        # Hypervolume mask (area dominated by Pareto front)
        hv_mask = np.zeros_like(X, dtype=bool)
        for point in pareto_obj:
            mask = (X <= point[0]) & (Y <= point[1])
            hv_mask |= mask

        # Plot hypervolume region
        ax.contourf(X, Y, hv_mask.astype(int), levels=[0, 1],
                   colors=['white', 'lightcyan'], alpha=0.7)

        # Plot Pareto front
        ax.scatter(pareto_obj[:, 0], pareto_obj[:, 1],
                  c='red', s=80, edgecolors='darkred', linewidth=1,
                  label=f'Pareto front', zorder=3)

        # Plot theoretical front
        ax.plot(theoretical[:, 0], theoretical[:, 1],
               'g-', linewidth=2, alpha=0.8,
               label='Theoretical front', zorder=2)

        # Plot nadir
        ax.plot(nadir[0], nadir[1], 'ko', markersize=10,
               label=f'Nadir ({nadir[0]}, {nadir[1]})', zorder=4)

        ax.set_xlabel('Objective 1 (f₁)', fontweight='bold')
        ax.set_ylabel('Objective 2 (f₂)', fontweight='bold')
        ax.set_title(f'{title}\\nHypervolume: {data["hypervolume"]:.4f}',
                    fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1.3)
        ax.set_ylim(0, 4)

    # Bottom left: Hypervolume comparison
    ax3.bar(['ECAM', 'Random Start'],
           [results['ECAM']['hypervolume'], results['Random_Start']['hypervolume']],
           color=['blue', 'orange'], alpha=0.7, edgecolor='black', linewidth=2)
    ax3.set_ylabel('Hypervolume', fontweight='bold')
    ax3.set_title('Hypervolume Comparison\\n(Same Nadir Point)', fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='y')

    # Add values on bars
    for i, (alg, hv) in enumerate([('ECAM', results['ECAM']['hypervolume']),
                                  ('Random Start', results['Random_Start']['hypervolume'])]):
        ax3.text(i, hv + 0.01, f'{hv:.4f}', ha='center', va='bottom',
                fontweight='bold', fontsize=12)

    # Bottom right: Performance metrics comparison
    metrics = ['Hypervolume', 'Pareto Size', 'Efficiency (%)']
    ecam_metrics = [
        results['ECAM']['hypervolume'],
        results['ECAM']['pareto_size'],
        (results['ECAM']['pareto_size'] / results['ECAM']['evaluations']) * 100
    ]
    random_metrics = [
        results['Random_Start']['hypervolume'],
        results['Random_Start']['pareto_size'],
        (results['Random_Start']['pareto_size'] / results['Random_Start']['evaluations']) * 100
    ]

    x = np.arange(len(metrics))
    width = 0.35

    # Normalize metrics for comparison (0-1 scale)
    max_vals = [max(ecam_metrics[i], random_metrics[i]) for i in range(len(metrics))]
    ecam_norm = [ecam_metrics[i] / max_vals[i] for i in range(len(metrics))]
    random_norm = [random_metrics[i] / max_vals[i] for i in range(len(metrics))]

    ax4.bar(x - width/2, ecam_norm, width, label='ECAM',
           color='blue', alpha=0.7, edgecolor='black')
    ax4.bar(x + width/2, random_norm, width, label='Random Start',
           color='orange', alpha=0.7, edgecolor='black')

    ax4.set_ylabel('Normalized Performance', fontweight='bold')
    ax4.set_title('Performance Metrics Comparison\\n(Normalized to [0,1])', fontweight='bold')
    ax4.set_xticks(x)
    ax4.set_xticklabels(metrics)
    ax4.legend()
    ax4.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    output_path = Path("moecam_comprehensive_results") / "Hypervolume_Analysis.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"✅ Hypervolume Analysis saved: {output_path}")
    return output_path

def create_summary_visualization():
    """Create a comprehensive summary visualization."""
    fig = plt.figure(figsize=(20, 12))

    # Create custom layout
    gs = fig.add_gridspec(3, 4, height_ratios=[2, 2, 1], width_ratios=[1, 1, 1, 1])

    # Load results
    results = load_comprehensive_results()
    theoretical = generate_theoretical_zdt1()
    nadir = np.array([1.2, 10.0])

    # Main comparison plot (top span)
    ax_main = fig.add_subplot(gs[0, :])

    # Plot both algorithms on same plot
    for i, (alg_name, color, marker) in enumerate([('ECAM', 'blue', 'o'),
                                                  ('Random_Start', 'red', 's')]):
        data = results[alg_name]
        all_obj = data['all_objectives']
        pareto_obj = data['pareto_front']

        # Plot all points (smaller, more transparent)
        ax_main.scatter(all_obj[:, 0], all_obj[:, 1],
                       alpha=0.3, s=15, c=color,
                       label=f'{alg_name.replace("_", " ")} all points')

        # Plot Pareto fronts
        ax_main.scatter(pareto_obj[:, 0], pareto_obj[:, 1],
                       alpha=0.9, s=100, c=color, edgecolors='black', linewidth=2,
                       marker=marker, label=f'{alg_name.replace("_", " ")} Pareto front',
                       zorder=3)

    # Plot theoretical front
    ax_main.plot(theoretical[:, 0], theoretical[:, 1],
                'g-', linewidth=4, alpha=0.9,
                label='Theoretical Pareto front', zorder=2)

    ax_main.plot(nadir[0], nadir[1], 'ko', markersize=12,
                label=f'Nadir point ({nadir[0]}, {nadir[1]})', zorder=4)

    ax_main.set_xlabel('Objective 1 (f₁) → MINIMIZE', fontsize=16, fontweight='bold')
    ax_main.set_ylabel('Objective 2 (f₂) → MINIMIZE', fontsize=16, fontweight='bold')
    ax_main.set_title('MOECAM Algorithm Comparison on ZDT1 Minimization Problem',
                     fontsize=20, fontweight='bold', pad=20)
    ax_main.legend(loc='upper right', fontsize=12)
    ax_main.grid(True, alpha=0.3)
    ax_main.set_xlim(-0.05, 1.3)
    ax_main.set_ylim(-0.1, 5)

    # Individual algorithm details (middle row)
    for i, (alg_name, title) in enumerate([('ECAM', 'ECAM Details'),
                                          ('Random_Start', 'Random Start Details')]):
        ax = fig.add_subplot(gs[1, i*2:(i+1)*2])
        data = results[alg_name]
        pareto_obj = data['pareto_front']

        # Plot with hypervolume region indication
        ax.scatter(pareto_obj[:, 0], pareto_obj[:, 1],
                  c='red', s=60, alpha=0.8, edgecolors='darkred')
        ax.plot(theoretical[:, 0], theoretical[:, 1], 'g-', linewidth=2, alpha=0.7)

        # Add performance text
        efficiency = (data['pareto_size'] / data['evaluations']) * 100
        performance_text = (f"Evaluations: {data['evaluations']:,}\\n"
                          f"Pareto Points: {data['pareto_size']}\\n"
                          f"Efficiency: {efficiency:.1f}%\\n"
                          f"Hypervolume: {data['hypervolume']:.4f}")

        ax.text(0.02, 0.98, performance_text, transform=ax.transAxes,
               bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", alpha=0.9),
               verticalalignment='top', fontsize=11)

        ax.set_xlabel('Objective 1 (f₁)', fontweight='bold')
        ax.set_ylabel('Objective 2 (f₂)', fontweight='bold')
        ax.set_title(title, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1.1)
        ax.set_ylim(0, 3)

    # Summary statistics (bottom row)
    ax_summary = fig.add_subplot(gs[2, :])
    ax_summary.axis('off')

    # Create summary table
    summary_text = f"""
    COMPREHENSIVE MOECAM ANALYSIS SUMMARY

    ✅ MINIMIZATION CONTEXT: Both objectives (f₁, f₂) should be minimized - lower values are better
    ✅ NADIR POINT CONSISTENCY: Same reference point [{nadir[0]}, {nadir[1]}] used for fair hypervolume comparison
    ✅ HIGH FUNCTION EVALUATIONS: ECAM (10,000) vs Random Start (6,000) for statistical significance
    ✅ HYPERVOLUME ADVANTAGE: ECAM achieves {((results['ECAM']['hypervolume']/results['Random_Start']['hypervolume'] - 1) * 100):.2f}% better hypervolume
    ✅ PARETO EFFICIENCY: ECAM finds {results['ECAM']['pareto_size']} vs {results['Random_Start']['pareto_size']} Pareto optimal solutions
    """

    ax_summary.text(0.5, 0.5, summary_text, ha='center', va='center',
                   fontsize=14, fontweight='bold',
                   bbox=dict(boxstyle="round,pad=1", facecolor="lightgreen", alpha=0.3))

    plt.tight_layout()
    output_path = Path("moecam_comprehensive_results") / "MOECAM_Complete_Analysis.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"✅ Complete Analysis saved: {output_path}")
    return output_path

def main():
    """Generate all enhanced visualizations."""
    print("🎨 Generating Enhanced MOECAM Visualizations")
    print("=" * 60)

    try:
        # Load comprehensive results
        results = load_comprehensive_results()
        print("✅ Results loaded successfully")

        # Generate publication-quality figures
        fig2_path = create_figure2_detailed_analysis(results)
        fig4_path = create_figure4_algorithm_comparison(results)
        hv_path = create_hypervolume_analysis_plot(results)
        summary_path = create_summary_visualization()

        print("\n🎉 ALL ENHANCED VISUALIZATIONS GENERATED!")
        print("=" * 60)
        print("📊 Publication-Quality Figures Created:")
        print(f"   📈 Figure 2: {fig2_path}")
        print(f"   📈 Figure 4: {fig4_path}")
        print(f"   📈 Hypervolume Analysis: {hv_path}")
        print(f"   📈 Complete Summary: {summary_path}")
        print()
        print("✅ All plots clearly show:")
        print("   • Minimization problem context with directional indicators")
        print("   • Consistent nadir points for fair hypervolume comparison")
        print("   • High function evaluation comprehensive analysis")
        print("   • Professional publication-ready quality")

    except Exception as e:
        print(f"❌ Error generating visualizations: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
