import numpy as np
import time
from ..problems.test_functions import ZDT1
from ..metrics.performance_metrics import pareto_front

try:
    from ..core.cffi_interface import CppOptimizer
    CPP_AVAILABLE = True
except ImportError:
    CPP_AVAILABLE = False
    print("Warning: C++ interface not available, using Python-only implementation")

class NSGAII:
    def __init__(self, problem, pop_size=100, num_generations=250, 
                 crossover_rate=0.9, mutation_rate=0.01, use_cpp=False):
        self.problem = problem
        self.pop_size = pop_size
        self.num_generations = num_generations
        self.crossover_rate = crossover_rate
        self.mutation_rate = mutation_rate
        self.use_cpp = use_cpp and CPP_AVAILABLE
        
        # Initialize C++ optimizer if requested and available
        if self.use_cpp:
            self.cpp_optimizer = CppOptimizer(
                problem.n_dim, 
                problem.n_obj,
                problem.xl.tolist(), 
                problem.xu.tolist()
            )
        
        # Performance tracking
        self.evaluation_count = 0
        self.execution_time = 0
        self.convergence_history = []

    def _initialize_population(self):
        if self.use_cpp:
            # Use C++ for population initialization
            population = []
            for _ in range(self.pop_size):
                individual = self.cpp_optimizer.generate_random_solution()
                population.append(individual)
            return np.array(population)
        else:
            # Python implementation
            return np.random.rand(self.pop_size, self.problem.n_dim)

    def _evaluate_population(self, population):
        objectives = np.array([self.problem.evaluate(ind) for ind in population])
        self.evaluation_count += len(population)
        return objectives

    def _fast_non_dominated_sort(self, objectives):
        """Enhanced non-dominated sorting with proper dominance relations"""
        n = objectives.shape[0]
        fronts = [[] for _ in range(n)]
        domination_count = np.zeros(n)
        dominated_solutions = [[] for _ in range(n)]
        
        # Calculate dominance relationships
        for i in range(n):
            for j in range(n):
                if i != j:
                    if self._dominates(objectives[i], objectives[j]):
                        dominated_solutions[i].append(j)
                    elif self._dominates(objectives[j], objectives[i]):
                        domination_count[i] += 1
        
        # Find first front
        current_front = []
        for i in range(n):
            if domination_count[i] == 0:
                current_front.append(i)
                fronts[0] = current_front
        
        # Build subsequent fronts
        front_idx = 0
        while len(fronts[front_idx]) > 0:
            next_front = []
            for i in fronts[front_idx]:
                for j in dominated_solutions[i]:
                    domination_count[j] -= 1
                    if domination_count[j] == 0:
                        next_front.append(j)
            front_idx += 1
            if front_idx < n:
                fronts[front_idx] = next_front
            else:
                break
        
        # Remove empty fronts
        return [front for front in fronts if len(front) > 0]
    
    def _dominates(self, obj1, obj2):
        """Check if obj1 dominates obj2 (minimization)"""
        return np.all(obj1 <= obj2) and np.any(obj1 < obj2)

    def _crowding_distance_assignment(self, objectives, front_indices):
        """Calculate crowding distance for diversity preservation"""
        n = len(front_indices)
        if n <= 2:
            return np.full(n, float('inf'))
        
        distances = np.zeros(n)
        obj_front = objectives[front_indices]
        
        for m in range(objectives.shape[1]):  # For each objective
            # Sort by objective m
            sorted_indices = np.argsort(obj_front[:, m])
            
            # Boundary points get infinite distance
            distances[sorted_indices[0]] = float('inf')
            distances[sorted_indices[-1]] = float('inf')
            
            # Calculate distance for interior points
            obj_range = obj_front[sorted_indices[-1], m] - obj_front[sorted_indices[0], m]
            if obj_range > 0:
                for i in range(1, n - 1):
                    distances[sorted_indices[i]] += (
                        obj_front[sorted_indices[i + 1], m] - 
                        obj_front[sorted_indices[i - 1], m]
                    ) / obj_range
        
        return distances

    def _crossover(self, parent1, parent2):
        """Enhanced crossover with C++ integration and SBX"""
        if self.use_cpp:
            # Use C++ crossover
            return self.cpp_optimizer.crossover(parent1, parent2)
        else:
            # Python SBX (Simulated Binary Crossover)
            eta_c = 20  # Distribution index for crossover
            child1, child2 = parent1.copy(), parent2.copy()
            
            for i in range(len(parent1)):
                if np.random.random() <= self.crossover_rate:
                    if abs(parent1[i] - parent2[i]) > 1e-14:
                        # SBX crossover
                        u = np.random.random()
                        if u <= 0.5:
                            beta = (2 * u) ** (1.0 / (eta_c + 1))
                        else:
                            beta = (1.0 / (2 * (1 - u))) ** (1.0 / (eta_c + 1))
                        
                        child1[i] = 0.5 * ((1 + beta) * parent1[i] + (1 - beta) * parent2[i])
                        child2[i] = 0.5 * ((1 - beta) * parent1[i] + (1 + beta) * parent2[i])
                        
                        # Bound check
                        child1[i] = np.clip(child1[i], self.problem.xl[i], self.problem.xu[i])
                        child2[i] = np.clip(child2[i], self.problem.xl[i], self.problem.xu[i])
            
            return child1, child2

    def _mutate(self, individual):
        """Enhanced mutation with C++ integration and polynomial mutation"""
        if self.use_cpp:
            # Use C++ mutation
            return self.cpp_optimizer.mutate(individual, self.mutation_rate)
        else:
            # Python polynomial mutation
            eta_m = 20  # Distribution index for mutation
            mutated = individual.copy()
            
            for i in range(len(individual)):
                if np.random.random() <= self.mutation_rate:
                    u = np.random.random()
                    if u < 0.5:
                        delta = (2 * u) ** (1.0 / (eta_m + 1)) - 1
                    else:
                        delta = 1 - (2 * (1 - u)) ** (1.0 / (eta_m + 1))
                    
                    mutated[i] = individual[i] + delta * (self.problem.xu[i] - self.problem.xl[i])
                    mutated[i] = np.clip(mutated[i], self.problem.xl[i], self.problem.xu[i])
            
            return mutated
    
    def _tournament_selection(self, population, objectives, fronts, crowding_distances):
        """Tournament selection based on dominance and crowding distance"""
        tournament_size = 2
        selected = []
        
        for _ in range(self.pop_size):
            candidates = np.random.choice(len(population), tournament_size, replace=False)
            
            # Find which front each candidate belongs to
            candidate_fronts = []
            for candidate in candidates:
                for front_idx, front in enumerate(fronts):
                    if candidate in front:
                        candidate_fronts.append(front_idx)
                        break
            
            # Select based on front rank first, then crowding distance
            if candidate_fronts[0] < candidate_fronts[1]:
                selected.append(candidates[0])
            elif candidate_fronts[1] < candidate_fronts[0]:
                selected.append(candidates[1])
            else:
                # Same front, select based on crowding distance
                front_idx = candidate_fronts[0]
                front = fronts[front_idx]
                pos1 = front.index(candidates[0])
                pos2 = front.index(candidates[1])
                
                if crowding_distances[pos1] > crowding_distances[pos2]:
                    selected.append(candidates[0])
                else:
                    selected.append(candidates[1])
        
        return selected

    def optimize(self):
        """Enhanced NSGA-II optimization with performance tracking"""
        start_time = time.time()
        
        # Initialize population
        population = self._initialize_population()
        objectives = self._evaluate_population(population)
        
        # Track convergence
        self.convergence_history = []
        all_solutions = []  # Store all evaluated solutions for visualization
        
        for generation in range(self.num_generations):
            # Non-dominated sorting
            fronts = self._fast_non_dominated_sort(objectives)
            
            # Calculate crowding distances for each front
            all_crowding_distances = []
            for front in fronts:
                if len(front) > 0:
                    distances = self._crowding_distance_assignment(objectives, front)
                    all_crowding_distances.append(distances)
                else:
                    all_crowding_distances.append([])
            
            # Generate offspring
            offspring_population = []
            selected_indices = self._tournament_selection(population, objectives, fronts, all_crowding_distances)
            
            for i in range(0, self.pop_size, 2):
                parent1_idx = selected_indices[i]
                parent2_idx = selected_indices[min(i + 1, self.pop_size - 1)]
                
                parent1 = population[parent1_idx]
                parent2 = population[parent2_idx]
                
                # Crossover and mutation
                child1, child2 = self._crossover(parent1, parent2)
                child1 = self._mutate(child1)
                child2 = self._mutate(child2)
                
                offspring_population.extend([child1, child2])
            
            # Ensure we have exactly pop_size offspring
            offspring_population = offspring_population[:self.pop_size]
            offspring_population = np.array(offspring_population)
            
            # Evaluate offspring
            offspring_objectives = self._evaluate_population(offspring_population)
            
            # Combine parent and offspring populations
            combined_population = np.vstack((population, offspring_population))
            combined_objectives = np.vstack((objectives, offspring_objectives))
            
            # Store all solutions for analysis
            all_solutions.extend(combined_objectives.tolist())
            
            # Environmental selection
            population, objectives = self._environmental_selection(
                combined_population, combined_objectives
            )
            
            # Track convergence metrics
            current_pareto = pareto_front(objectives)
            self.convergence_history.append({
                'generation': generation,
                'pareto_size': len(current_pareto),
                'hypervolume': self._calculate_hypervolume(current_pareto),
                'evaluations': self.evaluation_count
            })
            
            # Optional: Print progress
            if generation % 50 == 0 or generation == self.num_generations - 1:
                print(f"Generation {generation}: {len(current_pareto)} Pareto solutions, "
                      f"HV: {self.convergence_history[-1]['hypervolume']:.6f}")
        
        self.execution_time = time.time() - start_time
        self.all_evaluated_solutions = np.array(all_solutions)
        
        final_objectives = self._evaluate_population(population)
        return pareto_front(final_objectives)
    
    def _environmental_selection(self, population, objectives):
        """Environmental selection to maintain population size"""
        fronts = self._fast_non_dominated_sort(objectives)
        
        new_population = []
        new_objectives = []
        
        for front_idx, front in enumerate(fronts):
            if len(new_population) + len(front) <= self.pop_size:
                # Add entire front
                for idx in front:
                    new_population.append(population[idx])
                    new_objectives.append(objectives[idx])
            else:
                # Partial front inclusion based on crowding distance
                remaining_slots = self.pop_size - len(new_population)
                distances = self._crowding_distance_assignment(objectives, front)
                
                # Sort by crowding distance (descending)
                sorted_indices = np.argsort(distances)[::-1]
                
                for i in range(remaining_slots):
                    idx = front[sorted_indices[i]]
                    new_population.append(population[idx])
                    new_objectives.append(objectives[idx])
                break
        
        return np.array(new_population), np.array(new_objectives)
    
    def _calculate_hypervolume(self, pareto_front_points, reference_point=None):
        """Calculate hypervolume indicator"""
        if len(pareto_front_points) == 0:
            return 0.0
        
        if reference_point is None:
            # Use worst point + 10% margin as reference
            reference_point = np.max(pareto_front_points, axis=0) * 1.1
        
        # Simple 2D hypervolume calculation
        if pareto_front_points.shape[1] == 2:
            sorted_points = pareto_front_points[pareto_front_points[:, 0].argsort()]
            hv = 0.0
            
            for i, point in enumerate(sorted_points):
                if i == 0:
                    width = reference_point[0] - point[0]
                    height = reference_point[1] - point[1]
                else:
                    width = reference_point[0] - point[0]
                    height = reference_point[1] - point[1]
                    # Subtract overlapping area
                    prev_point = sorted_points[i-1]
                    if point[1] < prev_point[1]:
                        height = prev_point[1] - point[1]
                
                if width > 0 and height > 0:
                    hv += width * height
            
            return hv
        
        return 0.0  # For higher dimensions, return 0 for now
    
    def get_performance_metrics(self):
        """Get comprehensive performance metrics"""
        return {
            'total_evaluations': self.evaluation_count,
            'execution_time': self.execution_time,
            'convergence_history': self.convergence_history,
            'cpp_used': self.use_cpp
        }


class MOEAD:
    def __init__(self, problem, pop_size=100, num_generations=250):
        self.problem = problem
        self.pop_size = pop_size
        self.num_generations = num_generations
        self.weights = self._generate_weights()
        self.neighborhood = self._define_neighborhood()

    def _generate_weights(self):
        # Simple uniform weights for demonstration
        return np.linspace(0, 1, self.pop_size)[:, None]

    def _define_neighborhood(self):
        # Simple neighborhood for demonstration
        return {i: list(range(self.pop_size)) for i in range(self.pop_size)}

    def _tchebycheff(self, objectives, ideal_point, weight_vector):
        return np.max(weight_vector * np.abs(objectives - ideal_point), axis=1)

    def optimize(self):
        # Simplified MOEA/D for demonstration
        population = np.random.rand(self.pop_size, self.problem.n_dim)

        # Initialize ideal point by evaluating initial population
        initial_objectives = np.array([self.problem.evaluate(ind) for ind in population])
        ideal_point = np.min(initial_objectives, axis=0)

        for gen in range(self.num_generations):
            for i in range(self.pop_size):
                # Select two parents from neighborhood
                p1_idx, p2_idx = np.random.choice(self.neighborhood[i], 2, replace=False)
                parent1, parent2 = population[p1_idx], population[p2_idx]

                # Simplified crossover and mutation
                child = (parent1 + parent2) / 2 + np.random.normal(0, 0.1, self.problem.n_dim)
                child = np.clip(child, self.problem.xl, self.problem.xu)

                child_obj = self.problem.evaluate(child)
                ideal_point = np.minimum(ideal_point, child_obj)

                # Update neighboring solutions
                for neighbor_idx in self.neighborhood[i]:
                    neighbor_obj = self.problem.evaluate(population[neighbor_idx])
                    if self._tchebycheff(child_obj[None, :], ideal_point, self.weights[neighbor_idx]) < \
                       self._tchebycheff(neighbor_obj[None, :], ideal_point, self.weights[neighbor_idx]):
                        population[neighbor_idx] = child

        final_objectives = np.array([self.problem.evaluate(ind) for ind in population])
        return pareto_front(final_objectives)




