
import numpy as np
from ..problems.test_functions import ZDT1
from ..metrics.performance_metrics import pareto_front

class NSGAII:
    def __init__(self, problem, pop_size=100, num_generations=250, crossover_rate=0.9, mutation_rate=0.01):
        self.problem = problem
        self.pop_size = pop_size
        self.num_generations = num_generations
        self.crossover_rate = crossover_rate
        self.mutation_rate = mutation_rate

    def _initialize_population(self):
        return np.random.rand(self.pop_size, self.problem.n_dim)

    def _evaluate_population(self, population):
        return np.array([self.problem.evaluate(ind) for ind in population])

    def _fast_non_dominated_sort(self, objectives):
        # Simplified non-dominated sort for demonstration
        # In a real NSGA-II, this would be more sophisticated
        fronts = []
        remaining_indices = set(range(objectives.shape[0]))

        while remaining_indices:
            current_front = []
            dominated_by_count = np.zeros(objectives.shape[0])
            dominates_set = [set() for _ in range(objectives.shape[0])]

            for i in list(remaining_indices):
                for j in list(remaining_indices):
                    if i == j: continue
                    if np.all(objectives[i] <= objectives[j]) and np.any(objectives[i] < objectives[j]):
                        dominates_set[i].add(j)
                    elif np.all(objectives[j] <= objectives[i]) and np.any(objectives[j] < objectives[i]):
                        dominated_by_count[i] += 1

                if dominated_by_count[i] == 0:
                    current_front.append(i)

            fronts.append(current_front)
            for idx in current_front:
                remaining_indices.remove(idx)

        return fronts

    def _crowding_distance_assignment(self, objectives, front_indices):
        # Simplified crowding distance for demonstration
        # In a real NSGA-II, this would be more sophisticated
        distances = np.zeros(objectives.shape[0])
        return distances

    def _selection(self, population, objectives):
        # Simplified selection for demonstration
        # In a real NSGA-II, this would be more sophisticated
        return population

    def _crossover(self, parent1, parent2):
        # Simplified crossover for demonstration
        # In a real NSGA-II, this would be more sophisticated
        child1 = parent1.copy()
        child2 = parent2.copy()
        return child1, child2

    def _mutate(self, individual):
        # Simplified mutation for demonstration
        # In a real NSGA-II, this would be more sophisticated
        return individual

    def optimize(self):
        population = self._initialize_population()
        for gen in range(self.num_generations):
            offspring_population = []
            for i in range(0, self.pop_size, 2):
                parent1, parent2 = population[i], population[i+1] # Simplified selection
                child1, child2 = self._crossover(parent1, parent2)
                offspring_population.append(self._mutate(child1))
                offspring_population.append(self._mutate(child2))

            combined_population = np.vstack((population, np.array(offspring_population)))
            combined_objectives = self._evaluate_population(combined_population)

            fronts = self._fast_non_dominated_sort(combined_objectives)
            new_population = []
            new_population_objectives = []

            for front in fronts:
                if len(new_population) + len(front) <= self.pop_size:
                    new_population.extend([combined_population[i] for i in front])
                    new_population_objectives.extend([combined_objectives[i] for i in front])
                else:
                    # Need to select based on crowding distance
                    remaining_count = self.pop_size - len(new_population)
                    # Simplified: just take the first 'remaining_count' individuals from the front
                    new_population.extend([combined_population[i] for i in front[:remaining_count]])
                    new_population_objectives.extend([combined_objectives[i] for i in front[:remaining_count]])
                    break
            population = np.array(new_population)

        final_objectives = self._evaluate_population(population)
        return pareto_front(final_objectives)


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




