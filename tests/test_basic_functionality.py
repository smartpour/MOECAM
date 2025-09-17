#!/usr/bin/env python3
"""
Basic functionality tests for MOECAM package.
"""

import unittest
import numpy as np
import sys
import os

# Add src to path for testing
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'src'))

from moecam.problems.test_functions import ZDT1, ZDT2, DTLZ1
from moecam.algorithms.moea_algorithms import NSGAII, MOEAD
from moecam.metrics.performance_metrics import pareto_front, hypervolume, EvaluationCounter

class TestBasicFunctionality(unittest.TestCase):
    
    def test_zdt1_problem(self):
        """Test ZDT1 problem initialization and evaluation."""
        problem = ZDT1(n_dim=5)
        self.assertEqual(problem.n_dim, 5)
        self.assertEqual(problem.n_obj, 2)
        
        # Test evaluation
        x = np.array([0.5, 0.5, 0.5, 0.5, 0.5])
        objectives = problem.evaluate(x)
        self.assertEqual(len(objectives), 2)
        self.assertTrue(all(isinstance(obj, (int, float, np.number)) for obj in objectives))
    
    def test_pareto_front_calculation(self):
        """Test Pareto front extraction."""
        # Simple test case with known Pareto front
        points = np.array([
            [1.0, 3.0],
            [2.0, 2.0],
            [3.0, 1.0],
            [2.5, 2.5],  # Dominated point
            [1.5, 1.5]   # Pareto optimal
        ])
        
        pf = pareto_front(points)
        # Should contain points [1.0, 3.0], [2.0, 2.0], [3.0, 1.0], [1.5, 1.5]
        self.assertTrue(len(pf) >= 3)  # At least the non-dominated points
    
    def test_hypervolume_calculation(self):
        """Test hypervolume calculation."""
        # Simple 2D case
        points = np.array([
            [1.0, 2.0],
            [2.0, 1.0]
        ])
        reference_point = [3.0, 3.0]
        
        hv = hypervolume(points, reference_point)
        self.assertGreater(hv, 0)  # Should be positive
    
    def test_evaluation_counter(self):
        """Test evaluation counter functionality."""
        counter = EvaluationCounter()
        counter.start()
        
        # Simulate some evaluations
        for _ in range(10):
            counter.increment()
        
        counter.stop()
        
        self.assertEqual(counter.get_count(), 10)
        self.assertGreater(counter.get_time(), 0)
    
    def test_nsga2_algorithm(self):
        """Test NSGA-II algorithm initialization and basic functionality."""
        problem = ZDT1(n_dim=5)
        algorithm = NSGAII(problem, pop_size=20, num_generations=5)
        
        # Test that algorithm can be initialized
        self.assertEqual(algorithm.pop_size, 20)
        self.assertEqual(algorithm.num_generations, 5)
        
        # Test optimization (with very small parameters for speed)
        result = algorithm.optimize()
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(result.shape[1], 2)  # Should have 2 objectives
    
    def test_moead_algorithm(self):
        """Test MOEA/D algorithm initialization and basic functionality."""
        problem = ZDT1(n_dim=5)
        algorithm = MOEAD(problem, pop_size=20, num_generations=5)
        
        # Test that algorithm can be initialized
        self.assertEqual(algorithm.pop_size, 20)
        self.assertEqual(algorithm.num_generations, 5)
        
        # Test optimization (with very small parameters for speed)
        result = algorithm.optimize()
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(result.shape[1], 2)  # Should have 2 objectives

if __name__ == '__main__':
    unittest.main()

