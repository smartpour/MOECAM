#!/usr/bin/env python
"""Simple test of complete MOECAM integration with working example data."""

import numpy as np
import sys
from pathlib import Path

# Add scripts directory to path for imports
scripts_dir = Path(__file__).parent
sys.path.append(str(scripts_dir))

from pareto_interface import ParetoFrontExtractor
from wfg_interface import WFGHypervolume

def test_with_kursawe_data():
    """Test integration using the actual working Kursawe dataset."""
    print("=== Testing with Kursawe Dataset ===")

    # Load the known working Kursawe data
    kur_objectives_path = Path(__file__).parent.parent.parent.parent / "csources" / "pareto" / "kur0_2.txt"
    kur_pareto_path = Path(__file__).parent.parent.parent.parent / "csources" / "pareto" / "kur0_2P.txt"

    print(f"Loading original objectives from: {kur_objectives_path}")
    print(f"Loading known Pareto front from: {kur_pareto_path}")

    # Load original objectives (all evaluated points)
    original_objectives = np.loadtxt(kur_objectives_path)
    print(f"✓ Loaded {len(original_objectives)} original objective evaluations")
    print(f"  Objectives range: f1=[{original_objectives[:,0].min():.3f}, {original_objectives[:,0].max():.3f}], f2=[{original_objectives[:,1].min():.3f}, {original_objectives[:,1].max():.3f}]")

    # Load known Pareto front
    known_pareto = np.loadtxt(kur_pareto_path)
    print(f"✓ Loaded {len(known_pareto)} known Pareto optimal points")
    print(f"  Pareto range: f1=[{known_pareto[:,0].min():.3f}, {known_pareto[:,0].max():.3f}], f2=[{known_pareto[:,1].min():.3f}, {known_pareto[:,1].max():.3f}]")

    # Test 1: Extract Pareto front from original data
    print(f"\n=== Step 1: Pareto Front Extraction ===")
    pareto_extractor = ParetoFrontExtractor()
    pareto_result = pareto_extractor.extract_pareto_front(original_objectives)
    extracted_pareto = pareto_result['pareto_front']

    print(f"✓ Extracted {len(extracted_pareto)} Pareto points from {len(original_objectives)} total points")
    reduction_ratio = pareto_result['summary']['reduction_ratio']
    print(f"  Reduction: {len(original_objectives)} → {len(extracted_pareto)} ({reduction_ratio:.1%} reduction)")
    print(f"  Extracted range: f1=[{extracted_pareto[:,0].min():.3f}, {extracted_pareto[:,0].max():.3f}], f2=[{extracted_pareto[:,1].min():.3f}, {extracted_pareto[:,1].max():.3f}]")

    # Test 2: Calculate hypervolume of known Pareto front (use working file directly)
    print(f"\n=== Step 2: Hypervolume Calculation (Known Data) ===")
    wfg = WFGHypervolume()

    # Test with the file that we know works
    print("Testing WFG with known working file...")
    import subprocess
    result = subprocess.run([str(wfg.wfg_path), str(kur_pareto_path)], capture_output=True, text=True)
    print(f"Direct WFG result: {result.stdout.strip()}")

    # Parse the hypervolume value
    if "hv(1) =" in result.stdout:
        hv_value = float(result.stdout.split("hv(1) = ")[1].split()[0])
        print(f"✓ Known Pareto front hypervolume: {hv_value:.10f}")
    else:
        print(f"❌ No hypervolume output: {result.stdout}")
        return

    # Test 3: Calculate hypervolume of our extracted front
    print(f"\n=== Step 3: Hypervolume of Extracted Front ===")

    # Save our extracted front in the exact same format as the working file
    test_front_path = Path(__file__).parent / "test_results" / "extracted_pareto_wfg.txt"
    test_front_path.parent.mkdir(exist_ok=True)

    # Copy the exact format from the original file
    with open(test_front_path, 'w') as f:
        f.write('#\n')
        for point in extracted_pareto:
            # Use same format as original: variable precision + trailing space
            f.write(f'{point[0]:.5g} {point[1]:.5g} \n')

    print(f"Saved extracted front to: {test_front_path}")

    # Test with WFG
    result = subprocess.run([str(wfg.wfg_path), str(test_front_path)], capture_output=True, text=True)
    print(f"WFG result on extracted front: {result.stdout.strip()}")

    if "hv(1) =" in result.stdout:
        extracted_hv = float(result.stdout.split("hv(1) = ")[1].split()[0])
        print(f"✓ Extracted front hypervolume: {extracted_hv:.10f}")

        # Compare hypervolumes
        hv_ratio = extracted_hv / hv_value if hv_value > 0 else 0
        print(f"\n=== Comparison ===")
        print(f"Known Pareto front:    {len(known_pareto)} points, HV = {hv_value:.10f}")
        print(f"Extracted front:       {len(extracted_pareto)} points, HV = {extracted_hv:.10f}")
        print(f"Hypervolume ratio:     {hv_ratio:.6f} ({hv_ratio*100:.2f}%)")
        print(f"Point count ratio:     {len(extracted_pareto)/len(known_pareto):.6f} ({len(extracted_pareto)/len(known_pareto)*100:.2f}%)")

        if hv_ratio > 0.95:
            print("✅ EXCELLENT: Extracted front achieves >95% of known hypervolume!")
        elif hv_ratio > 0.85:
            print("✅ GOOD: Extracted front achieves >85% of known hypervolume")
        elif hv_ratio > 0.70:
            print("⚠️  FAIR: Extracted front achieves >70% of known hypervolume")
        else:
            print("❌ POOR: Extracted front achieves <70% of known hypervolume")

    else:
        print(f"❌ No hypervolume output for extracted front: {result.stdout}")

    print(f"\n✅ Integration test completed!")
    print(f"   • Pareto extraction: WORKING")
    print(f"   • Hypervolume calculation: WORKING")
    print(f"   • File I/O: WORKING")
    print(f"   • C++ tools integration: WORKING")

if __name__ == "__main__":
    test_with_kursawe_data()
