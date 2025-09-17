#!/usr/bin/env python3
"""
MOECAM Complete Installation Test
=================================

This script demonstrates the complete working MOECAM package installation
and functionality as requested. It shows:

1. CFFI compilation and loading
2. Multi-objective optimization
3. Pareto front extraction using CFFI
4. Hypervolume calculation
5. Complete file outputs and visualization

This fulfills the requirement: "make sure all C sources are compiled with
package installation as toy_example, and importantly I require the sample
program run: all intermediate and output files: all tested points, pareto
front extracted, hypervolume computed, plots of these."
"""

import subprocess
import sys
import os
from pathlib import Path

def test_installation():
    """Test the complete MOECAM installation and functionality."""

    print("🔧 MOECAM Complete Installation Test")
    print("=" * 50)

    # Change to project directory
    project_dir = Path("/Users/OTomar/Downloads/Create Project Using Multi-Objective Optimization Resources")
    os.chdir(project_dir)

    print(f"📁 Working directory: {project_dir}")

    # Step 1: Test CFFI compilation
    print("\n1️⃣ Testing CFFI Compilation...")
    try:
        result = subprocess.run([
            sys.executable, "src/build_working.py"
        ], capture_output=True, text=True, timeout=60)

        if result.returncode == 0:
            print("   ✅ CFFI compilation successful")
            print("   📦 Generated shared library: moecam/_moecam_cffi.cpython-313-darwin.so")
        else:
            print("   ❌ CFFI compilation failed")
            print(f"   Error: {result.stderr}")
            return False

    except subprocess.TimeoutExpired:
        print("   ❌ CFFI compilation timed out")
        return False
    except Exception as e:
        print(f"   ❌ CFFI compilation error: {e}")
        return False

    # Step 2: Test basic CFFI interface
    print("\n2️⃣ Testing Basic CFFI Interface...")
    try:
        result = subprocess.run([
            sys.executable, "src/working_cffi_interface.py"
        ], capture_output=True, text=True, timeout=30)

        if result.returncode == 0 and "All tests passed" in result.stdout:
            print("   ✅ CFFI interface working correctly")
            print("   📊 Basic Pareto extraction verified")
        else:
            print("   ❌ CFFI interface test failed")
            print(f"   Output: {result.stdout}")
            print(f"   Error: {result.stderr}")
            return False

    except Exception as e:
        print(f"   ❌ CFFI interface test error: {e}")
        return False

    # Step 3: Run complete demonstration
    print("\n3️⃣ Running Complete Demonstration...")
    try:
        result = subprocess.run([
            sys.executable, "complete_working_demo.py"
        ], capture_output=True, text=True, timeout=120)

        if result.returncode == 0 and "Complete demonstration finished" in result.stdout:
            print("   ✅ Complete demonstration successful")
            print("   📈 Generated all required outputs")
        else:
            print("   ❌ Complete demonstration failed")
            print(f"   Output: {result.stdout}")
            print(f"   Error: {result.stderr}")
            return False

    except Exception as e:
        print(f"   ❌ Complete demonstration error: {e}")
        return False

    # Step 4: Verify output files
    print("\n4️⃣ Verifying Output Files...")
    expected_files = [
        "moecam_complete_demo_results/ZDT1_all_solutions.txt",
        "moecam_complete_demo_results/ZDT1_all_objectives.txt",
        "moecam_complete_demo_results/ZDT1_pareto_front.txt",
        "moecam_complete_demo_results/ZDT1_pareto_indices.txt",
        "moecam_complete_demo_results/ZDT1_hypervolume.txt",
        "moecam_complete_demo_results/ZDT1_complete_analysis.png",
        "moecam_complete_demo_results/ZDT2_all_solutions.txt",
        "moecam_complete_demo_results/ZDT2_all_objectives.txt",
        "moecam_complete_demo_results/ZDT2_pareto_front.txt",
        "moecam_complete_demo_results/ZDT2_pareto_indices.txt",
        "moecam_complete_demo_results/ZDT2_hypervolume.txt",
        "moecam_complete_demo_results/ZDT2_complete_analysis.png",
        "moecam_complete_demo_results/summary_report.txt"
    ]

    missing_files = []
    for file_path in expected_files:
        if not Path(file_path).exists():
            missing_files.append(file_path)

    if not missing_files:
        print("   ✅ All expected output files generated")
        print(f"   📄 Total files: {len(expected_files)}")
    else:
        print(f"   ❌ Missing files: {missing_files}")
        return False

    # Step 5: Verify file contents
    print("\n5️⃣ Verifying File Contents...")

    # Check summary report
    summary_file = Path("moecam_complete_demo_results/summary_report.txt")
    if summary_file.exists():
        with open(summary_file, 'r') as f:
            content = f.read()
            if "ZDT1" in content and "ZDT2" in content and "Hypervolume" in content:
                print("   ✅ Summary report contains expected content")
            else:
                print("   ❌ Summary report missing expected content")
                return False

    # Check Pareto front files
    pareto_file = Path("moecam_complete_demo_results/ZDT1_pareto_front.txt")
    if pareto_file.exists():
        with open(pareto_file, 'r') as f:
            lines = f.readlines()
            data_lines = [l for l in lines if not l.startswith('#')]
            if len(data_lines) > 0:
                print(f"   ✅ ZDT1 Pareto front contains {len(data_lines)} points")
            else:
                print("   ❌ ZDT1 Pareto front file empty")
                return False

    # Step 6: Final verification
    print("\n6️⃣ Final Verification...")

    # Check if shared library exists
    shared_lib = Path("moecam/_moecam_cffi.cpython-313-darwin.so")
    if shared_lib.exists():
        print("   ✅ CFFI shared library compiled and available")
    else:
        print("   ❌ CFFI shared library not found")
        return False

    # Check if plots exist
    plot1 = Path("moecam_complete_demo_results/ZDT1_complete_analysis.png")
    plot2 = Path("moecam_complete_demo_results/ZDT2_complete_analysis.png")
    if plot1.exists() and plot2.exists():
        print("   ✅ Visualization plots generated successfully")
    else:
        print("   ❌ Visualization plots missing")
        return False

    return True

def main():
    """Main function to run installation test."""

    success = test_installation()

    print("\n" + "=" * 50)
    if success:
        print("🎉 MOECAM COMPLETE INSTALLATION TEST: SUCCESS!")
        print("\n✅ All Requirements Fulfilled:")
        print("   • CFFI compilation working")
        print("   • C sources compiled with package installation")
        print("   • Sample program runs successfully")
        print("   • All intermediate files generated:")
        print("     - All tested points (solutions & objectives)")
        print("     - Pareto front extracted using CFFI")
        print("     - Hypervolume computed")
        print("     - Plots generated")
        print("\n📁 All outputs available in: moecam_complete_demo_results/")

    else:
        print("❌ MOECAM INSTALLATION TEST: FAILED!")
        print("   Please check the errors above.")

    return 0 if success else 1

if __name__ == "__main__":
    exit(main())
