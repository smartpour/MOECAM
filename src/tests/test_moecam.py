#!/usr/bin/env python3
"""
Test script to verify MOECAM installation and functionality
"""

import sys
import os

def test_basic_import():
    """Test basic imports"""
    try:
        import moecam
        print("✓ Successfully imported moecam")
        print(f"✓ Version: {moecam.__version__}")
        print(f"✓ Author: {moecam.__author__}")
        return True
    except ImportError as e:
        print(f"✗ Failed to import moecam: {e}")
        return False

def test_module_imports():
    """Test module imports"""
    modules_to_test = [
        'moecam.problems.test_functions',
        'moecam.algorithms.moea_algorithms',
        'moecam.metrics.performance_metrics',
        'moecam.core',
        'moecam.utils'
    ]

    successful = 0
    for module in modules_to_test:
        try:
            __import__(module)
            print(f"✓ Successfully imported {module}")
            successful += 1
        except ImportError as e:
            print(f"✗ Failed to import {module}: {e}")

    return successful == len(modules_to_test)

def test_problem_creation():
    """Test creating a test problem"""
    try:
        from moecam.problems.test_functions import ZDT1
        problem = ZDT1(n_dim=5)
        print(f"✓ Successfully created ZDT1 problem with {problem.n_dim} dimensions")

        # Test evaluation
        import numpy as np
        x = np.random.rand(5)
        objectives = problem.evaluate(x)
        print(f"✓ Problem evaluation works, got {len(objectives)} objectives")
        return True
    except Exception as e:
        print(f"✗ Failed to create/test problem: {e}")
        return False

def test_algorithm_creation():
    """Test creating an algorithm"""
    try:
        from moecam.problems.test_functions import ZDT1
        from moecam.algorithms.moea_algorithms import NSGAII

        problem = ZDT1(n_dim=5)
        algorithm = NSGAII(problem, pop_size=20, num_generations=5)
        print("✓ Successfully created NSGA-II algorithm")
        return True
    except Exception as e:
        print(f"✗ Failed to create algorithm: {e}")
        return False

def run_mini_optimization():
    """Run a very small optimization test"""
    try:
        from moecam.problems.test_functions import ZDT1
        from moecam.algorithms.moea_algorithms import NSGAII

        print("🔄 Running mini optimization test...")

        # Create problem and algorithm
        problem = ZDT1(n_dim=5)
        algorithm = NSGAII(problem, pop_size=10, num_generations=3)

        # Run optimization
        pareto_front = algorithm.optimize()

        print(f"✓ Optimization complete! Found {len(pareto_front)} solutions")
        return True

    except Exception as e:
        print(f"✗ Error during optimization: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests"""
    print("Testing MOECAM Installation and Functionality...")
    print("=" * 50)

    tests = [
        ("Basic Import", test_basic_import),
        ("Module Imports", test_module_imports),
        ("Problem Creation", test_problem_creation),
        ("Algorithm Creation", test_algorithm_creation),
        ("Mini Optimization", run_mini_optimization)
    ]

    passed = 0
    for test_name, test_func in tests:
        print(f"\n{test_name}:")
        print("-" * 20)
        if test_func():
            passed += 1

    print("\n" + "=" * 50)
    print(f"Tests passed: {passed}/{len(tests)}")

    if passed == len(tests):
        print("🎉 All tests passed! MOECAM is ready to use.")
    else:
        print("⚠️  Some tests failed. Check the error messages above.")

    return passed == len(tests)

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
