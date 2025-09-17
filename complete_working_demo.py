#!/usr/bin/env python3
"""
MOECAM Complete Working Test Program
===================================

Complete demonstration showing:
- Multi-objective test function evaluation
- Simple optimization (random sampling)
- Pareto front extraction
- Simple hypervolume calculation
- Visualization and file outputs

This demonstrates the core functionality requested in the requirements.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
import sys

# Add the src directory to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

try:
    from working_cffi_interface import extract_pareto_front_cffi, CFFI_AVAILABLE
except ImportError:
    CFFI_AVAILABLE = False


class MOECAMWorkingDemo:
    """Complete working demonstration of MOECAM functionality."""

    def __init__(self, output_dir: str = "moecam_working_demo_results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.results = {}

    def zdt1_function(self, x: np.ndarray) -> np.ndarray:
        """ZDT1 test function (2 objectives)."""
        x = np.asarray(x)
        if x.ndim == 1:
            x = x.reshape(1, -1)

        n_points, n_dim = x.shape
        objectives = np.zeros((n_points, 2))

        for i in range(n_points):
            f1 = x[i, 0]
            g = 1 + 9 * np.sum(x[i, 1:]) / (n_dim - 1) if n_dim > 1 else 1
            f2 = g * (1 - np.sqrt(f1 / g)) if f1 <= g else g
            objectives[i] = [f1, f2]

        return objectives.squeeze() if n_points == 1 else objectives

    def zdt2_function(self, x: np.ndarray) -> np.ndarray:
        """ZDT2 test function (2 objectives)."""
        x = np.asarray(x)
        if x.ndim == 1:
            x = x.reshape(1, -1)

        n_points, n_dim = x.shape
        objectives = np.zeros((n_points, 2))

        for i in range(n_points):
            f1 = x[i, 0]
            g = 1 + 9 * np.sum(x[i, 1:]) / (n_dim - 1) if n_dim > 1 else 1
            f2 = g * (1 - (f1 / g)**2) if f1 <= g else g
            objectives[i] = [f1, f2]

        return objectives.squeeze() if n_points == 1 else objectives

    def simple_random_optimization(self, objective_func, bounds, max_evaluations=1000):
        """Simple random sampling optimization algorithm."""
        print(f"   Running random optimization with {max_evaluations} evaluations...")

        n_dim = len(bounds)

        # Generate random solutions
        solutions = np.zeros((max_evaluations, n_dim))
        objectives = np.zeros((max_evaluations, 2))

        for i in range(max_evaluations):
            # Generate random solution within bounds
            solution = np.array([
                np.random.uniform(bounds[j][0], bounds[j][1])
                for j in range(n_dim)
            ])

            # Evaluate objectives
            obj = objective_func(solution)

            solutions[i] = solution
            objectives[i] = obj

        return {
            'solutions': solutions,
            'objectives': objectives,
            'num_evaluations': max_evaluations,
            'algorithm': 'Random Sampling'
        }

    def calculate_simple_hypervolume(self, points: np.ndarray, reference_point: np.ndarray) -> float:
        """Simple hypervolume calculation using Monte Carlo method."""
        # Simple 2D hypervolume calculation
        if points.shape[1] != 2:
            raise ValueError("Simple hypervolume only works for 2 objectives")

        # Sort points by first objective
        sorted_points = points[np.argsort(points[:, 0])]

        # Calculate hypervolume using rectangles
        hypervolume = 0.0
        prev_x = 0.0  # Assuming minimization starting from 0

        for point in sorted_points:
            x, y = point
            if x > prev_x:
                # Calculate rectangle area
                width = x - prev_x
                height = reference_point[1] - y if y < reference_point[1] else 0
                hypervolume += width * height
                prev_x = x

        # Add final rectangle
        if prev_x < reference_point[0]:
            width = reference_point[0] - prev_x
            height = reference_point[1] - sorted_points[-1, 1] if sorted_points[-1, 1] < reference_point[1] else 0
            hypervolume += width * height

        return hypervolume

    def save_results_to_file(self, filename: str, data: np.ndarray, header: str = ""):
        """Save numerical results to file."""
        filepath = self.output_dir / filename
        with open(filepath, 'w') as f:
            if header:
                f.write(f"# {header}\n")
            if data.ndim == 1:
                for val in data:
                    f.write(f"{val:.6f}\n")
            else:
                for row in data:
                    f.write(" ".join(f"{val:.6f}" for val in row) + "\n")
        print(f"   Saved to {filepath}")

    def plot_results(self, all_objectives: np.ndarray, pareto_objectives: np.ndarray,
                    true_pareto: np.ndarray, title: str, filename: str):
        """Create visualization plots."""
        plt.figure(figsize=(12, 8))

        # Plot all evaluated points
        plt.scatter(all_objectives[:, 0], all_objectives[:, 1],
                   alpha=0.6, s=20, c='lightblue', label='All Evaluated Points')

        # Plot Pareto front
        plt.scatter(pareto_objectives[:, 0], pareto_objectives[:, 1],
                   alpha=0.8, s=50, c='red', label='Extracted Pareto Front')

        # Plot true Pareto front
        if true_pareto is not None:
            plt.plot(true_pareto[:, 0], true_pareto[:, 1],
                    'g-', linewidth=2, label='True Pareto Front')

        plt.xlabel('Objective 1 (f1)')
        plt.ylabel('Objective 2 (f2)')
        plt.title(title)
        plt.legend()
        plt.grid(True, alpha=0.3)

        # Save plot
        filepath = self.output_dir / filename
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   Plot saved to {filepath}")

    def generate_true_pareto_front(self, func_name: str, n_points: int = 100) -> np.ndarray:
        """Generate true Pareto front for comparison."""
        if func_name == "ZDT1":
            f1 = np.linspace(0, 1, n_points)
            f2 = 1 - np.sqrt(f1)
            return np.column_stack([f1, f2])
        elif func_name == "ZDT2":
            f1 = np.linspace(0, 1, n_points)
            f2 = 1 - f1**2
            return np.column_stack([f1, f2])
        else:
            return None

    def test_function(self, func_name: str, objective_func, bounds, max_evaluations: int = 500):
        """Test a specific function with complete workflow."""
        print(f"\n=== Testing {func_name} Function ===")

        # 1. Run optimization
        print("📊 Running optimization...")
        result = self.simple_random_optimization(objective_func, bounds, max_evaluations)

        all_solutions = result['solutions']
        all_objectives = result['objectives']

        print(f"   ✅ Completed {result['num_evaluations']} evaluations")

        # 2. Save all points
        print("💾 Saving evaluation data...")
        self.save_results_to_file(
            f"{func_name}_all_solutions.txt",
            all_solutions,
            f"{func_name} - All evaluated solutions"
        )

        self.save_results_to_file(
            f"{func_name}_all_objectives.txt",
            all_objectives,
            f"{func_name} - All objective evaluations"
        )

        # 3. Extract Pareto front
        print("🔍 Extracting Pareto front...")
        if not CFFI_AVAILABLE:
            print("   ❌ CFFI not available, using simple Python implementation")
            # Simple Python Pareto extraction
            is_pareto = np.ones(len(all_objectives), dtype=bool)
            for i in range(len(all_objectives)):
                for j in range(len(all_objectives)):
                    if i != j:
                        if np.all(all_objectives[j] <= all_objectives[i]) and np.any(all_objectives[j] < all_objectives[i]):
                            is_pareto[i] = False
                            break
            pareto_objectives = all_objectives[is_pareto]
            pareto_indices = np.where(is_pareto)[0]
        else:
            pareto_objectives, pareto_indices = extract_pareto_front_cffi(all_objectives, strict_mode=False)

        print(f"   ✅ Found {len(pareto_objectives)} Pareto optimal points")
        print(f"   📈 Pareto ratio: {len(pareto_objectives)/len(all_objectives):.2%}")

        # 4. Save Pareto front
        self.save_results_to_file(
            f"{func_name}_pareto_front.txt",
            pareto_objectives,
            f"{func_name} - Pareto Front"
        )

        self.save_results_to_file(
            f"{func_name}_pareto_indices.txt",
            pareto_indices,
            f"{func_name} - Pareto Front Indices"
        )

        # 5. Calculate hypervolume
        print("📏 Calculating hypervolume...")
        reference_point = np.array([1.1, 1.1])  # Reference point
        try:
            hypervolume = self.calculate_simple_hypervolume(pareto_objectives, reference_point)
            print(f"   ✅ Hypervolume: {hypervolume:.6f}")

            # Save hypervolume
            with open(self.output_dir / f"{func_name}_hypervolume.txt", 'w') as f:
                f.write(f"# {func_name} - Hypervolume calculation\n")
                f.write(f"# Reference point: {reference_point}\n")
                f.write(f"{hypervolume:.6f}\n")

        except Exception as e:
            print(f"   ❌ Hypervolume calculation failed: {e}")
            hypervolume = 0.0

        # 6. Generate visualization
        print("📈 Creating visualization...")
        true_pareto = self.generate_true_pareto_front(func_name)

        self.plot_results(
            all_objectives, pareto_objectives, true_pareto,
            f"{func_name} - Complete Analysis (Random Sampling)",
            f"{func_name}_complete_analysis.png"
        )

        return {
            'function': func_name,
            'all_solutions': all_solutions,
            'all_objectives': all_objectives,
            'pareto_objectives': pareto_objectives,
            'pareto_indices': pareto_indices,
            'hypervolume': hypervolume,
            'reference_point': reference_point,
            'num_evaluations': result['num_evaluations']
        }

    def run_complete_demonstration(self):
        """Run the complete MOECAM working demonstration."""
        print("🚀 MOECAM Complete Working Demonstration")
        print("=" * 50)

        if CFFI_AVAILABLE:
            print("✅ CFFI library available - using optimized Pareto extraction")
        else:
            print("⚠️  CFFI library not available - using Python fallback")

        # Test configurations
        test_configs = [
            {
                'name': 'ZDT1',
                'function': self.zdt1_function,
                'bounds': [(0, 1)] * 5,  # 5-dimensional
                'evaluations': 1000
            },
            {
                'name': 'ZDT2',
                'function': self.zdt2_function,
                'bounds': [(0, 1)] * 3,  # 3-dimensional
                'evaluations': 500
            }
        ]

        # Run tests
        for config in test_configs:
            result = self.test_function(
                config['name'],
                config['function'],
                config['bounds'],
                config['evaluations']
            )

            self.results[config['name']] = result

        # Generate summary
        self.generate_summary_report()

        print(f"\n🎉 Complete demonstration finished!")
        print(f"📁 All results saved to: {self.output_dir.absolute()}")
        print("\n📋 Generated files:")
        for file in sorted(self.output_dir.glob("*")):
            print(f"   {file.name}")

    def generate_summary_report(self):
        """Generate summary report."""
        print("\n📋 Generating summary report...")

        summary_file = self.output_dir / "summary_report.txt"
        with open(summary_file, 'w') as f:
            f.write("MOECAM Complete Working Demonstration - Summary Report\n")
            f.write("=" * 60 + "\n\n")

            for func_name, result in self.results.items():
                f.write(f"Function: {func_name}\n")
                f.write(f"  Evaluations: {result['num_evaluations']}\n")
                f.write(f"  Total points: {len(result['all_objectives'])}\n")
                f.write(f"  Pareto points: {len(result['pareto_objectives'])}\n")
                f.write(f"  Pareto ratio: {len(result['pareto_objectives'])/len(result['all_objectives']):.2%}\n")
                f.write(f"  Hypervolume: {result['hypervolume']:.6f}\n")
                f.write(f"  Reference point: {result['reference_point']}\n")
                f.write("\n")

            f.write("Generated Files:\n")
            f.write("- *_all_solutions.txt: All evaluated decision variables\n")
            f.write("- *_all_objectives.txt: All objective function evaluations\n")
            f.write("- *_pareto_front.txt: Extracted Pareto optimal points\n")
            f.write("- *_pareto_indices.txt: Indices of Pareto optimal points\n")
            f.write("- *_hypervolume.txt: Hypervolume calculation results\n")
            f.write("- *_complete_analysis.png: Visualization plots\n")

        print(f"   Summary saved to {summary_file}")


def main():
    """Main function to run the complete demonstration."""
    print("MOECAM Complete Working Demonstration")
    print("=" * 40)

    # Set random seed for reproducible results
    np.random.seed(42)

    # Create and run demonstration
    demo = MOECAMWorkingDemo("moecam_complete_demo_results")
    demo.run_complete_demonstration()


if __name__ == "__main__":
    main()
