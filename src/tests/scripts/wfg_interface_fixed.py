#!/usr/bin/env python
"""Python interface for the WFG hypervolume calculation tool.

This module provides a Python wrapper around the WFG hypervolume calculator,
allowing computation of hypervolume indicators for Pareto fronts.
"""

import numpy as np
import subprocess
import tempfile
import os
from pathlib import Path

class WFGHypervolume:
    """Python interface for WFG hypervolume calculation."""

    def __init__(self, wfg_executable_path=None):
        """Initialize WFG hypervolume calculator.

        Args:
            wfg_executable_path: Path to compiled wfg executable.
                               If None, attempts to compile from source.
        """
        self.wfg_path = wfg_executable_path
        # Get the root directory of the project (3 levels up from scripts)
        root_dir = Path(__file__).parent.parent.parent.parent
        self.source_dir = root_dir / "csources" / "wfg" / "WFG_1.15"

        if self.wfg_path is None:
            self.wfg_path = self._compile_executable()

    def _compile_executable(self):
        """Compile the WFG hypervolume calculator."""
        exe_path = self.source_dir / "wfg"

        if exe_path.exists():
            print(f"✓ Using existing WFG executable: {exe_path}")
            return str(exe_path)

        print(f"Compiling WFG hypervolume calculator in {self.source_dir}")
        try:
            # Use the provided makefile
            result = subprocess.run(
                ["make"],
                cwd=str(self.source_dir),
                check=True,
                capture_output=True,
                text=True
            )
            print(f"✓ Compiled successfully: {exe_path}")
            return str(exe_path)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to compile WFG: {e.stderr}")

    def calculate_hypervolume(self, fronts, reference_point=None):
        """Calculate hypervolume for one or more Pareto fronts.

        Args:
            fronts: Single front (2D array) or list of fronts
            reference_point: Reference point for hypervolume calculation.
                           If None, uses origin (0, 0, ..., 0)

        Returns:
            dict: {
                'hypervolumes': list of hypervolume values,
                'total_hypervolume': sum of all hypervolumes,
                'reference_point': used reference point,
                'summary': dict with statistics
            }
        """
        # Normalize input to list of fronts
        if isinstance(fronts, np.ndarray) and fronts.ndim == 2:
            fronts = [fronts]

        if not fronts:
            return {'hypervolumes': [], 'total_hypervolume': 0.0, 'reference_point': None}

        # Determine dimensionality
        first_front = np.asarray(fronts[0])
        if first_front.ndim != 2:
            raise ValueError("Each front must be a 2D array (n_points, n_objectives)")

        n_objectives = first_front.shape[1]

        # Set default reference point
        if reference_point is None:
            reference_point = np.zeros(n_objectives)
        else:
            reference_point = np.asarray(reference_point)
            if len(reference_point) != n_objectives:
                raise ValueError(f"Reference point must have {n_objectives} objectives")

        # Create temporary input file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.wfg', delete=False) as input_file:
            input_path = input_file.name

            # Write fronts in WFG format
            for front in fronts:
                front = np.asarray(front)
                input_file.write('#\n')  # Start of front marker
                for point in front:
                    input_file.write(' '.join(f'{val:.12f}' for val in point) + '\n')

        try:
            # Build command: wfg FRONTFILE [r1 r2 ... rd]
            cmd = [str(self.wfg_path), input_path]
            if not np.allclose(reference_point, 0):
                cmd.extend(f'{val:.12f}' for val in reference_point)

            # Run WFG hypervolume calculator
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

            if result.returncode != 0:
                raise RuntimeError(f"WFG hypervolume calculation failed: {result.stderr}")

            # Parse hypervolume results from output
            hypervolumes = self._parse_wfg_output(result.stdout)

            # Calculate summary
            total_hv = sum(hypervolumes) if hypervolumes else 0.0
            summary = {
                'num_fronts': len(fronts),
                'total_points': sum(len(f) for f in fronts),
                'objectives': n_objectives,
                'avg_hypervolume': total_hv / len(hypervolumes) if hypervolumes else 0.0
            }

            return {
                'hypervolumes': hypervolumes,
                'total_hypervolume': total_hv,
                'reference_point': reference_point,
                'summary': summary
            }

        finally:
            # Cleanup
            if os.path.exists(input_path):
                os.unlink(input_path)

    def _parse_wfg_output(self, output_text):
        """Parse hypervolume values from WFG output."""
        hypervolumes = []

        for line in output_text.strip().split('\n'):
            line = line.strip()
            # Look for lines like "hv(1) = 12.345678"
            if line.startswith('hv(') and ' = ' in line:
                try:
                    # Extract the hypervolume value after " = "
                    hv_str = line.split(' = ')[1]
                    hv = float(hv_str)
                    hypervolumes.append(hv)
                except (ValueError, IndexError):
                    continue  # Skip malformed lines

        return hypervolumes

    def save_fronts_for_wfg(self, fronts, filename):
        """Save Pareto fronts in WFG format.

        Args:
            fronts: Single front (2D array) or list of fronts
            filename: Output filename
        """
        if isinstance(fronts, np.ndarray) and fronts.ndim == 2:
            fronts = [fronts]

        with open(filename, 'w') as f:
            for front in fronts:
                front = np.asarray(front)
                f.write('#\n')  # Front separator
                for point in front:
                    f.write(' '.join(f'{val:.12f}' for val in point) + '\n')

        print(f"Saved {len(fronts)} front(s) to {filename}")

    def compare_fronts(self, fronts_dict, reference_point=None):
        """Compare multiple Pareto fronts by hypervolume.

        Args:
            fronts_dict: Dict mapping names to Pareto fronts
            reference_point: Reference point for comparison

        Returns:
            dict: Comparison results with rankings
        """
        results = {}

        for name, front in fronts_dict.items():
            hv_result = self.calculate_hypervolume([front], reference_point)
            results[name] = {
                'hypervolume': hv_result['hypervolumes'][0] if hv_result['hypervolumes'] else 0.0,
                'num_points': len(front),
                'front': front
            }

        # Rank by hypervolume (higher is better)
        ranked = sorted(results.items(), key=lambda x: x[1]['hypervolume'], reverse=True)

        comparison = {
            'results': results,
            'ranking': ranked,
            'reference_point': reference_point
        }

        return comparison

def demo_hypervolume_calculation():
    """Demonstrate hypervolume calculation with test data."""
    print("=== WFG Hypervolume Calculation Demo ===")

    # Create test Pareto fronts
    np.random.seed(42)

    # Front 1: Simple convex front
    n_points1 = 50
    t = np.linspace(0, 1, n_points1)
    front1 = np.column_stack([t, (1-t)**2])

    # Front 2: Different trade-off
    n_points2 = 40
    t = np.linspace(0, 1, n_points2)
    front2 = np.column_stack([t**2, 1-t])

    # Front 3: Dominated front (worse)
    front3 = front1 + 0.5  # Shift entire front

    fronts_dict = {
        'Convex Front': front1,
        'Concave Front': front2,
        'Dominated Front': front3
    }

    # Calculate hypervolumes
    wfg = WFGHypervolume()

    # Set reference point (slightly worse than any front)
    ref_point = np.array([2.0, 2.0])

    comparison = wfg.compare_fronts(fronts_dict, ref_point)

    # Display results
    print(f"\nHypervolume Comparison (reference point: {ref_point}):")
    print(f"{'Rank':<4} {'Front Name':<15} {'Hypervolume':<12} {'Points':<8}")
    print("-" * 45)

    for rank, (name, data) in enumerate(comparison['ranking'], 1):
        hv = data['hypervolume']
        pts = data['num_points']
        print(f"{rank:<4} {name:<15} {hv:<12.6f} {pts:<8}")

    # Test individual hypervolume calculation
    print(f"\n=== Individual Front Analysis ===")
    for name, front in fronts_dict.items():
        result = wfg.calculate_hypervolume(front, ref_point)
        hv = result['hypervolumes'][0] if result['hypervolumes'] else 0.0
        print(f"{name}: HV = {hv:.6f} ({len(front)} points)")

    # Save test data
    test_dir = Path(__file__).parent / "test_results"
    test_dir.mkdir(exist_ok=True)

    wfg.save_fronts_for_wfg(
        list(fronts_dict.values()),
        test_dir / "demo_fronts_for_wfg.txt"
    )

    return comparison

if __name__ == "__main__":
    demo_hypervolume_calculation()
