#!/usr/bin/env python
"""Complete integration demo showing MOECAM objective recording → Pareto extraction → Hypervolume calculation."""

import numpy as np
import sys
from pathlib import Path

# Add scripts directory to path for imports
scripts_dir = Path(__file__).parent
sys.path.append(str(scripts_dir))

from pareto_interface import ParetoFrontExtractor
from wfg_interface import WFGHypervolume

def kursawe_function(x):
    """Kursawe multi-objective test function."""
    f1 = sum(-10 * np.exp(-0.2 * np.sqrt(x[i]**2 + x[i+1]**2)) for i in range(len(x)-1))
    f2 = sum(abs(x[i])**0.8 + 5 * np.sin(x[i]**3) for i in range(len(x)))
    return np.array([f1, f2])

def zdt1_function(x):
    """ZDT1 multi-objective test function."""
    f1 = x[0]
    g = 1 + 9 * np.sum(x[1:]) / (len(x) - 1)
    h = 1 - np.sqrt(f1 / g)
    f2 = g * h
    return np.array([f1, f2])

def record_objectives_during_optimization():
    """Simulate objective recording during multi-objective optimization."""
    print("=== Simulating Multi-Objective Optimization ===")

    # Generate test points for both problems
    np.random.seed(42)
    n_points = 1000

    # Kursawe problem: 3 variables in [-5, 5]
    kur_points = np.random.uniform(-5, 5, (n_points, 3))
    kur_objectives = np.array([kursawe_function(x) for x in kur_points])

    # ZDT1 problem: 5 variables in [0, 1]
    zdt_points = np.random.uniform(0, 1, (n_points, 5))
    zdt_objectives = np.array([zdt1_function(x) for x in zdt_points])

    print(f"Generated {n_points} points for each problem")
    print(f"Kursawe objectives range: f1=[{kur_objectives[:,0].min():.3f}, {kur_objectives[:,0].max():.3f}], f2=[{kur_objectives[:,1].min():.3f}, {kur_objectives[:,1].max():.3f}]")
    print(f"ZDT1 objectives range: f1=[{zdt_objectives[:,0].min():.3f}, {zdt_objectives[:,0].max():.3f}], f2=[{zdt_objectives[:,1].min():.3f}, {zdt_objectives[:,1].max():.3f}]")

    return {
        'Kursawe': {'points': kur_points, 'objectives': kur_objectives},
        'ZDT1': {'points': zdt_points, 'objectives': zdt_objectives}
    }

def complete_pareto_hypervolume_analysis():
    """Complete analysis: objectives → Pareto fronts → hypervolume comparison."""
    print("\n" + "="*60)
    print("COMPLETE MULTI-OBJECTIVE ANALYSIS PIPELINE")
    print("="*60)

    # Step 1: Record objectives during optimization
    problems = record_objectives_during_optimization()

    # Step 2: Extract Pareto fronts
    print(f"\n=== Step 2: Pareto Front Extraction ===")
    pareto_extractor = ParetoFrontExtractor()

    pareto_fronts = {}
    for prob_name, data in problems.items():
        print(f"\nExtracting Pareto front for {prob_name}...")
        result = pareto_extractor.extract_pareto_front(data['objectives'])

        pareto_fronts[prob_name] = result['pareto_front']
        print(f"✓ {prob_name}: {len(result['pareto_front'])} Pareto optimal points "
              f"from {len(data['objectives'])} total points "
              f"({result['reduction_ratio']:.1%} reduction)")

    # Step 3: Calculate hypervolumes
    print(f"\n=== Step 3: Hypervolume Calculation ===")
    wfg = WFGHypervolume()

    # For fair comparison, we need appropriate reference points
    # Use worst point in each dimension across all fronts
    all_objectives = np.vstack([front for front in pareto_fronts.values()])
    reference_point = np.max(all_objectives, axis=0) + 0.1  # Slightly worse than worst

    print(f"Reference point for hypervolume: {reference_point}")

    # Calculate hypervolumes for each front
    hv_comparison = wfg.compare_fronts(pareto_fronts, reference_point)

    # Display detailed results
    print(f"\n=== HYPERVOLUME COMPARISON RESULTS ===")
    print(f"{'Rank':<4} {'Problem':<12} {'Hypervolume':<15} {'Pareto Points':<14} {'HV per Point':<12}")
    print("-" * 70)

    for rank, (name, data) in enumerate(hv_comparison['ranking'], 1):
        hv = data['hypervolume']
        pts = data['num_points']
        hv_per_pt = hv / pts if pts > 0 else 0
        print(f"{rank:<4} {name:<12} {hv:<15.8f} {pts:<14} {hv_per_pt:<12.8f}")

    # Step 4: Detailed analysis
    print(f"\n=== DETAILED ANALYSIS ===")
    for prob_name in problems.keys():
        front = pareto_fronts[prob_name]
        hv_data = hv_comparison['results'][prob_name]

        print(f"\n{prob_name} Problem:")
        print(f"  • Pareto front size: {len(front)} points")
        print(f"  • Hypervolume: {hv_data['hypervolume']:.8f}")
        print(f"  • Objective ranges: f1=[{front[:,0].min():.3f}, {front[:,0].max():.3f}], f2=[{front[:,1].min():.3f}, {front[:,1].max():.3f}]")
        print(f"  • Front spread: f1_range={front[:,0].max()-front[:,0].min():.3f}, f2_range={front[:,1].max()-front[:,1].min():.3f}")

    # Step 5: Save results for further analysis
    results_dir = Path(__file__).parent / "test_results"
    results_dir.mkdir(exist_ok=True)

    # Save all Pareto fronts
    for prob_name, front in pareto_fronts.items():
        np.savetxt(results_dir / f"{prob_name.lower()}_pareto_front.txt", front, fmt='%.12f',
                   header=f"{prob_name} Pareto Front - {len(front)} points")

    # Save combined analysis
    wfg.save_fronts_for_wfg(list(pareto_fronts.values()), results_dir / "all_pareto_fronts_wfg.txt")

    print(f"\n✓ Results saved to {results_dir}/")
    print(f"  - Individual Pareto fronts: *_pareto_front.txt")
    print(f"  - WFG format file: all_pareto_fronts_wfg.txt")

    return {
        'problems': problems,
        'pareto_fronts': pareto_fronts,
        'hypervolume_comparison': hv_comparison,
        'reference_point': reference_point
    }

def quick_integration_test():
    """Quick test to verify all components work together."""
    print("=== Quick Integration Test ===")

    # Simple test data
    test_objectives = np.array([
        [1.0, 5.0],
        [2.0, 4.0],
        [3.0, 3.0],
        [4.0, 2.0],
        [5.0, 1.0],
        [1.5, 4.5],  # Dominated
        [2.5, 3.5],  # Dominated
    ])

    print(f"Test objectives ({len(test_objectives)} points):")
    for i, obj in enumerate(test_objectives):
        print(f"  Point {i+1}: f1={obj[0]:.1f}, f2={obj[1]:.1f}")

    # Extract Pareto front
    pareto_extractor = ParetoFrontExtractor()
    pareto_result = pareto_extractor.extract_pareto_front(test_objectives)
    pareto_front = pareto_result['pareto_front']

    print(f"\nPareto front ({len(pareto_front)} points):")
    for i, obj in enumerate(pareto_front):
        print(f"  Pareto {i+1}: f1={obj[0]:.1f}, f2={obj[1]:.1f}")

    # Calculate hypervolume
    wfg = WFGHypervolume()
    hv_result = wfg.calculate_hypervolume(pareto_front, reference_point=[6.0, 6.0])

    print(f"\nHypervolume calculation:")
    print(f"  Reference point: [6.0, 6.0]")
    print(f"  Hypervolume: {hv_result['hypervolumes'][0]:.6f}")

    return pareto_result, hv_result

if __name__ == "__main__":
    # Run quick test first
    print("MOECAM Integration Pipeline Demo")
    print("=" * 40)

    try:
        quick_integration_test()
        print("\n✓ Quick test passed!\n")

        # Run complete analysis
        results = complete_pareto_hypervolume_analysis()

        print(f"\n" + "="*60)
        print("PIPELINE COMPLETE - All tools working together!")
        print("="*60)

    except Exception as e:
        print(f"\n❌ Error in pipeline: {e}")
        import traceback
        traceback.print_exc()
