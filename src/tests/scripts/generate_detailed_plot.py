#!/usr/bin/env python3
"""
Generate Individual Algorithm Plots for Clear Visualization
==========================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from moecam_direct_interface import MOECAMDirectInterface

def zdt1(x):
    f1 = x[0]
    g = 1 + 9 * np.sum(x[1:]) / (len(x) - 1)
    f2 = g * (1 - np.sqrt(f1 / g))
    return [f1, f2]

def main():
    print('Generating Individual High-Quality Algorithm Plots')
    print('=' * 60)

    moecam = MOECAMDirectInterface()

    # Run ECAM on ZDT1 with more iterations for better results
    print('Running ECAM on ZDT1 (150 iterations)...')
    result = moecam.minimize_ecam(zdt1, [(0, 1)] * 3, num_objectives=2, max_iterations=150)

    if result['success']:
        all_objectives = np.array(result['all_objectives'])
        pareto_front = result['pareto_front']

        print(f'Total evaluations: {len(all_objectives)}')
        print(f'Pareto points: {len(pareto_front)}')
        print(f'Hypervolume: {result["hypervolume"]:.6f}')

        # Create high-quality individual plot
        plt.figure(figsize=(12, 8))

        # Plot all evaluated points (light blue, smaller)
        plt.scatter(all_objectives[:, 0], all_objectives[:, 1],
                   alpha=0.5, s=30, c='lightblue',
                   label=f'All ECAM evaluations ({len(all_objectives)} points)')

        # Plot Pareto front (red, larger, outlined)
        plt.scatter(pareto_front[:, 0], pareto_front[:, 1],
                   c='red', s=120, marker='D',
                   label=f'Pareto front ({len(pareto_front)} points)',
                   edgecolors='darkred', linewidth=2, zorder=5)

        # Connect Pareto points with a line to show the front
        pareto_sorted = pareto_front[np.argsort(pareto_front[:, 0])]
        plt.plot(pareto_sorted[:, 0], pareto_sorted[:, 1],
                'r--', alpha=0.7, linewidth=2, label='Pareto front curve')

        plt.title('ZDT1 Problem: ECAM Algorithm Results\\n' +
                 f'Evaluations: {len(all_objectives)}, Pareto Points: {len(pareto_front)}, ' +
                 f'Efficiency: {len(pareto_front)/len(all_objectives):.1%}, HV: {result["hypervolume"]:.3f}',
                 fontsize=14, fontweight='bold')
        plt.xlabel('Objective 1 (f₁)', fontsize=12)
        plt.ylabel('Objective 2 (f₂)', fontsize=12)
        plt.legend(fontsize=11)
        plt.grid(True, alpha=0.3)

        # Add some annotation
        plt.text(0.02, 0.98,
                f'Algorithm: ECAM\\nProblem: ZDT1\\nDimensions: 3D → 2D objectives\\nPareto Efficiency: {len(pareto_front)/len(all_objectives):.1%}',
                transform=plt.gca().transAxes, fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.8))

        plt.tight_layout()
        plt.savefig('ZDT1_ECAM_Detailed_Results.png', dpi=300, bbox_inches='tight')
        plt.savefig('ZDT1_ECAM_Detailed_Results.pdf', bbox_inches='tight')

        print('✅ High-quality plot saved: ZDT1_ECAM_Detailed_Results.png')
        print('✅ PDF version saved: ZDT1_ECAM_Detailed_Results.pdf')

        # Save detailed output files
        np.savetxt('ZDT1_DETAILED_all_evaluations.csv', all_objectives,
                   header='f1,f2 - All objective values evaluated by ECAM on ZDT1',
                   fmt='%.8f', delimiter=',')

        np.savetxt('ZDT1_DETAILED_pareto_front.csv', pareto_front,
                   header='f1,f2 - Pareto front extracted from ECAM on ZDT1',
                   fmt='%.8f', delimiter=',')

        print('✅ CSV files saved: ZDT1_DETAILED_all_evaluations.csv and ZDT1_DETAILED_pareto_front.csv')

    else:
        print(f'❌ Algorithm failed: {result["message"]}')

if __name__ == '__main__':
    main()
