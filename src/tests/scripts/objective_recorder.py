#!/usr/bin/env python
"""Objective recording utility for tracking optimization progress."""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

class ObjectiveRecorder:
    """Records objective function evaluations and provides analysis tools."""

    def __init__(self, problem_name="problem"):
        self.problem_name = problem_name
        self.evaluations = []
        self.objectives = []
        self.eval_count = 0

    def record(self, x, f):
        """Record an evaluation: x (decision vars) -> f (objectives)."""
        self.evaluations.append(np.array(x).copy())
        self.objectives.append(np.array(f).copy())
        self.eval_count += 1

    def get_data(self):
        """Get recorded data as arrays."""
        if not self.evaluations:
            return np.array([]), np.array([])
        X = np.array(self.evaluations)
        F = np.array(self.objectives)
        return X, F

    def pareto_front(self):
        """Extract Pareto front from recorded objectives."""
        X, F = self.get_data()
        if len(F) == 0 or F.shape[1] != 2:
            return np.array([]), np.array([])

        n_points = F.shape[0]
        is_pareto = np.ones(n_points, dtype=bool)

        for i in range(n_points):
            if is_pareto[i]:
                dominated = np.all(F <= F[i], axis=1) & np.any(F < F[i], axis=1)
                is_pareto[dominated] = False

        return X[is_pareto], F[is_pareto]

    def save_data(self, filename=None):
        """Save objective data to file."""
        if filename is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            results_dir = os.path.join(os.path.dirname(__file__), '..', 'test_results')
            os.makedirs(results_dir, exist_ok=True)
            filename = os.path.join(results_dir, f"{self.problem_name}_objectives_{timestamp}.txt")

        X, F = self.get_data()
        if len(F) > 0:
            # Save objectives with header
            header = f"Objectives for {self.problem_name} ({len(F)} evaluations)"
            if F.shape[1] == 2:
                header += "\nf1 f2"
            np.savetxt(filename, F, header=header)
            print(f"Saved {len(F)} objective values to {filename}")
        return filename

    def plot(self, show_pareto=True, save_plot=False):
        """Plot objectives with optional Pareto front."""
        X, F = self.get_data()
        if len(F) == 0:
            print("No data to plot")
            return

        if F.shape[1] != 2:
            print(f"Can only plot 2D objectives, got {F.shape[1]}")
            return

        plt.figure(figsize=(8, 6))
        plt.scatter(F[:, 0], F[:, 1], alpha=0.6, s=30, label=f'All points ({len(F)})')

        if show_pareto:
            X_pf, F_pf = self.pareto_front()
            if len(F_pf) > 0:
                # Sort for line plot
                sorted_idx = np.argsort(F_pf[:, 0])
                F_sorted = F_pf[sorted_idx]

                plt.scatter(F_pf[:, 0], F_pf[:, 1], color='red', s=50,
                           label=f'Pareto front ({len(F_pf)})', zorder=5)
                plt.plot(F_sorted[:, 0], F_sorted[:, 1], 'r-', alpha=0.8, linewidth=2, zorder=4)

        plt.xlabel('f1')
        plt.ylabel('f2')
        plt.title(f'{self.problem_name} - Objective Space')
        plt.legend()
        plt.grid(True, alpha=0.3)

        if save_plot:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            results_dir = os.path.join(os.path.dirname(__file__), '..', 'test_results')
            os.makedirs(results_dir, exist_ok=True)
            plotname = os.path.join(results_dir, f"{self.problem_name}_pareto_{timestamp}.png")
            plt.savefig(plotname, dpi=150, bbox_inches='tight')
            print(f"Saved plot to {plotname}")

        return plt.gca()

    def summary(self):
        """Print summary statistics."""
        X, F = self.get_data()
        X_pf, F_pf = self.pareto_front()

        print(f"\n=== {self.problem_name} Summary ===")
        print(f"Total evaluations: {len(F)}")
        if len(F) > 0:
            print(f"Decision variables: {X.shape[1]}")
            print(f"Objectives: {F.shape[1]}")
            print(f"Pareto optimal points: {len(F_pf)}")
            if len(F) > 0:
                print(f"Pareto efficiency: {len(F_pf)/len(F)*100:.1f}%")

# Example usage with recorded objectives
def demo_objective_recording():
    """Demonstrate objective recording with sphere problem."""

    # Create recorder
    recorder = ObjectiveRecorder("bi_sphere")

    def sphere_with_recording(x):
        """Sphere objectives with automatic recording."""
        f = np.array([np.sum(x**2), np.sum((x-1)**2)])
        recorder.record(x, f)  # Record every evaluation
        return f

    # Use with optimizer
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
    from moecam.core.cffi_interface import CppOptimizer

    opt = CppOptimizer(3, 2, [-2, -2, -2], [2, 2, 2], sphere_with_recording)
    X, F = opt.sample_random(300)  # This will trigger recording

    # Analyze results
    recorder.summary()
    recorder.save_data()
    ax = recorder.plot(save_plot=True)

    return recorder

if __name__ == "__main__":
    import sys
    import os

    print("=== Objective Recording Demo ===")
    recorder = demo_objective_recording()

    # Show the plot
    plt.show()
