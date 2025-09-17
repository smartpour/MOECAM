#!/usr/bin/env python3
"""
Simple MOECAM Test - Works with installed package
================================================

This test demonstrates the working CFFI functionality after installation.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

try:
    import moecam
    print(f"✅ MOECAM imported successfully")
    print(f"✅ CFFI available: {moecam.CFFI_AVAILABLE}")

    if moecam.CFFI_AVAILABLE:
        # Test basic functionality
        print(moecam.hello_moecam())

        # Test Pareto front extraction
        print("\n🔍 Testing Pareto front extraction...")

        # Generate test data (ZDT1-like)
        np.random.seed(42)
        n_points = 100
        objectives = np.random.rand(n_points, 2)

        # Make some points clearly Pareto optimal
        for i in range(10):
            objectives[i] = [i/10.0, 1 - i/10.0]

        # Extract Pareto front
        pareto_points, pareto_indices = moecam.extract_pareto_front(objectives)

        print(f"   📊 Total points: {len(objectives)}")
        print(f"   🎯 Pareto points: {len(pareto_points)}")
        print(f"   📈 Pareto ratio: {len(pareto_points)/len(objectives):.1%}")

        # Create visualization
        plt.figure(figsize=(10, 6))
        plt.scatter(objectives[:, 0], objectives[:, 1], alpha=0.6, label='All Points')
        plt.scatter(pareto_points[:, 0], pareto_points[:, 1], color='red', s=50, label='Pareto Front')
        plt.xlabel('Objective 1')
        plt.ylabel('Objective 2')
        plt.title('MOECAM CFFI Test - Pareto Front Extraction')
        plt.legend()
        plt.grid(True, alpha=0.3)

        output_dir = Path("moecam_installation_test_results")
        output_dir.mkdir(exist_ok=True)

        plt.savefig(output_dir / "cffi_test.png", dpi=300, bbox_inches='tight')
        plt.close()

        # Save results
        np.savetxt(output_dir / "all_points.txt", objectives, header="All objective points")
        np.savetxt(output_dir / "pareto_points.txt", pareto_points, header="Pareto optimal points")
        np.savetxt(output_dir / "pareto_indices.txt", pareto_indices, fmt='%d', header="Pareto point indices")

        print(f"\n💾 Results saved to: {output_dir.absolute()}")
        print("🎉 MOECAM CFFI test completed successfully!")

    else:
        print("❌ CFFI not available - package not properly installed")

except ImportError as e:
    print(f"❌ Failed to import MOECAM: {e}")
    print("💡 Try running: pip install -e moecam_package/")
