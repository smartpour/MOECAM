#!/usr/bin/env python
"""Python interface for the Pareto front extraction tool.

This module provides a Python wrapper around the C++ paretofront executable,
allowing seamless integration with MOECAM optimization results.
"""

import numpy as np
import subprocess
import tempfile
import os
from pathlib import Path

class ParetoFrontExtractor:
    """Python interface for C++ Pareto front extraction tool."""

    def __init__(self, cpp_executable_path=None):
        """Initialize Pareto front extractor.

        Args:
            cpp_executable_path: Path to compiled paretofront executable.
                                If None, attempts to compile from source.
        """
        self.cpp_path = cpp_executable_path
        # Get the root directory of the project (4 levels up from src/src/moecam/core)
        root_dir = Path(__file__).parent.parent.parent.parent
        self.source_dir = root_dir / "csources" / "pareto"

        if self.cpp_path is None:
            self.cpp_path = self._compile_executable()

    def _compile_executable(self):
        """Compile the C++ Pareto front extractor."""
        # First check for local executable in the same directory
        local_exe = Path(__file__).parent / "paretofront"
        if local_exe.exists():
            print(f"✓ Using local executable: {local_exe}")
            return str(local_exe)

        exe_path = self.source_dir / "paretofront"
        source_path = self.source_dir / "paretofront.cpp"

        if exe_path.exists():
            print(f"✓ Using existing executable: {exe_path}")
            return str(exe_path)

        print(f"Compiling Pareto front extractor from {source_path}")
        try:
            subprocess.run([
                "g++", "-std=c++17", "-O3",
                str(source_path), "-o", str(exe_path)
            ], check=True, capture_output=True, text=True)
            print(f"✓ Compiled successfully: {exe_path}")
            return str(exe_path)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to compile Pareto extractor: {e.stderr}")

    def extract_pareto_front(self, objectives, strict=False, include_indices=False):
        """Extract Pareto front from objective data.

        Args:
            objectives: Array of shape (n_points, n_objectives)
            strict: If True, use strict Pareto dominance (default: weak)
            include_indices: If True, return original indices of Pareto points

        Returns:
            dict: {
                'pareto_front': np.array of Pareto optimal objectives,
                'indices': np.array of original indices (if include_indices=True),
                'summary': dict with statistics
            }
        """
        objectives = np.asarray(objectives)
        if objectives.ndim != 2:
            raise ValueError("Objectives must be 2D array (n_points, n_objectives)")

        n_points, n_objectives = objectives.shape

        # Create temporary files
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as input_file:
            input_path = input_file.name
            # Write objectives in required format (no indices, just objectives)
            for obj_vec in objectives:
                input_file.write(' '.join(f'{obj:.12f}' for obj in obj_vec) + '\n')

        output_path = input_path + '.pareto'

        try:
            # Run C++ Pareto extractor
            # Usage: paretofront inputfile outputfile dim ndata strict cols
            cmd = [
                str(self.cpp_path),
                input_path,
                output_path,
                str(n_objectives),  # dimensions
                str(n_points),      # number of data points
                str(1 if strict else 0),  # strict dominance
                str(0)              # cols (no index columns in input)
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

            if result.returncode != 0:
                raise RuntimeError(f"Pareto extractor failed: {result.stderr}")

            # Read results
            pareto_front = self._read_pareto_output(output_path, n_objectives)
            indices = None

            if include_indices:
                indices = self._read_pareto_indices(output_path + '.label')

            # Calculate statistics
            efficiency = len(pareto_front) / n_points * 100 if n_points > 0 else 0
            reduction_ratio = (n_points - len(pareto_front)) / n_points if n_points > 0 else 0
            summary = {
                'total_points': n_points,
                'pareto_points': len(pareto_front),
                'efficiency_percent': efficiency,
                'reduction_ratio': reduction_ratio,
                'strict_dominance': strict,
                'objectives': n_objectives
            }

            result_dict = {
                'pareto_front': pareto_front,
                'summary': summary
            }

            if include_indices:
                result_dict['indices'] = indices

            return result_dict

        finally:
            # Cleanup temporary files
            for path in [input_path, output_path, output_path + '.label', output_path + '.test']:
                if os.path.exists(path):
                    os.unlink(path)

    def _read_pareto_output(self, output_path, n_objectives):
        """Read Pareto front from output file."""
        if not os.path.exists(output_path):
            return np.array([]).reshape(0, n_objectives)

        pareto_points = []
        with open(output_path, 'r') as f:
            reading_points = False
            for line in f:
                line = line.strip()
                if line == '#':
                    if reading_points:
                        break  # End of Pareto points
                    else:
                        reading_points = True  # Start of Pareto points
                elif reading_points and line:
                    values = list(map(float, line.split()))
                    if len(values) >= n_objectives:
                        pareto_points.append(values[:n_objectives])

        return np.array(pareto_points) if pareto_points else np.array([]).reshape(0, n_objectives)

    def _read_pareto_indices(self, label_path):
        """Read original indices of Pareto points."""
        if not os.path.exists(label_path):
            return np.array([])

        indices = []
        with open(label_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    indices.append(int(line))

        return np.array(indices)

    def save_objectives_for_pareto(self, objectives, filename, include_indices=False):
        """Save objectives in format suitable for Pareto extractor.

        Args:
            objectives: Array of shape (n_points, n_objectives)
            filename: Output filename
            include_indices: If True, include point indices in first column
        """
        objectives = np.asarray(objectives)

        with open(filename, 'w') as f:
            for i, obj_vec in enumerate(objectives):
                if include_indices:
                    f.write(f'{i} ')
                f.write(' '.join(f'{obj:.12f}' for obj in obj_vec) + '\n')

        print(f"Saved {len(objectives)} objective vectors to {filename}")

def demo_pareto_extraction():
    """Demonstrate Pareto front extraction with test data."""
    print("=== Pareto Front Extraction Demo ===")

    # Create test data (bi-objective)
    np.random.seed(42)
    n_points = 1000

    # Generate some test objectives (sphere-like functions)
    x = np.random.uniform(-2, 2, (n_points, 3))
    objectives = np.column_stack([
        np.sum(x**2, axis=1),           # f1: sum of squares
        np.sum((x - 1)**2, axis=1)      # f2: sum of squares shifted
    ])

    # Extract Pareto front
    extractor = ParetoFrontExtractor()
    result = extractor.extract_pareto_front(objectives, include_indices=True)

    # Display results
    summary = result['summary']
    print(f"\nResults:")
    print(f"  Total points: {summary['total_points']}")
    print(f"  Pareto optimal: {summary['pareto_points']}")
    print(f"  Efficiency: {summary['efficiency_percent']:.2f}%")

    if len(result['pareto_front']) > 0:
        pf = result['pareto_front']
        print(f"\nPareto front range:")
        print(f"  f1: [{pf[:,0].min():.6f}, {pf[:,0].max():.6f}]")
        print(f"  f2: [{pf[:,1].min():.6f}, {pf[:,1].max():.6f}]")

        print(f"\nFirst few Pareto points:")
        for i in range(min(5, len(pf))):
            idx = result['indices'][i] if 'indices' in result else i
            print(f"  Point {idx}: f1={pf[i,0]:.6f}, f2={pf[i,1]:.6f}")

    # Save test data for external processing
    test_dir = Path(__file__).parent / "test_results"
    test_dir.mkdir(exist_ok=True)

    extractor.save_objectives_for_pareto(
        objectives,
        test_dir / "demo_objectives_for_pareto.txt"
    )

    return result

if __name__ == "__main__":
    demo_pareto_extraction()
