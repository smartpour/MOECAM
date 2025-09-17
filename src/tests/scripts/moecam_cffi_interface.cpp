
#include "moecam_cffi_interface.h"
#include <vector>
#include <random>
#include <algorithm>
#include <cstring>

static std::vector<moecam_config_t> configs;
static std::random_device rd;
static std::mt19937 gen(rd());

int moecam_init(moecam_config_t* config) {
    configs.push_back(*config);
    return configs.size() - 1; // Return config ID
}

int moecam_run_ecam(int config_id, double* best_solution, double* best_objectives,
                   double* all_solutions, double* all_objectives, int* num_evaluations) {
    if (config_id >= configs.size()) return -1;

    moecam_config_t& config = configs[config_id];

    // Simple ECAM-like algorithm: random sampling with basic evolution
    int pop_size = 20;
    int total_evals = 0;
    std::vector<std::vector<double>> population(pop_size, std::vector<double>(config.n_dim));
    std::vector<std::vector<double>> objectives(pop_size, std::vector<double>(config.n_obj));

    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Initialize population
    for (int i = 0; i < pop_size; i++) {
        for (int j = 0; j < config.n_dim; j++) {
            population[i][j] = config.lower_bounds[j] +
                              dis(gen) * (config.upper_bounds[j] - config.lower_bounds[j]);
        }

        // Evaluate
        config.objective_func(population[i].data(), config.n_dim, objectives[i].data(), config.n_obj);

        // Store in all_solutions and all_objectives
        memcpy(&all_solutions[total_evals * config.n_dim], population[i].data(),
               config.n_dim * sizeof(double));
        memcpy(&all_objectives[total_evals * config.n_obj], objectives[i].data(),
               config.n_obj * sizeof(double));
        total_evals++;

        if (total_evals >= config.max_iterations) break;
    }

    // Evolution loop
    while (total_evals < config.max_iterations) {
        // Simple mutation and evaluation
        int parent_idx = total_evals % pop_size;
        std::vector<double> offspring = population[parent_idx];

        // Mutate
        for (int j = 0; j < config.n_dim; j++) {
            if (dis(gen) < 0.1) { // 10% mutation rate
                offspring[j] = config.lower_bounds[j] +
                              dis(gen) * (config.upper_bounds[j] - config.lower_bounds[j]);
            }
        }

        // Evaluate offspring
        std::vector<double> offspring_obj(config.n_obj);
        config.objective_func(offspring.data(), config.n_dim, offspring_obj.data(), config.n_obj);

        // Store
        memcpy(&all_solutions[total_evals * config.n_dim], offspring.data(),
               config.n_dim * sizeof(double));
        memcpy(&all_objectives[total_evals * config.n_obj], offspring_obj.data(),
               config.n_obj * sizeof(double));
        total_evals++;
    }

    // Find best solution (minimize first objective for simplicity)
    int best_idx = 0;
    for (int i = 1; i < total_evals; i++) {
        if (all_objectives[i * config.n_obj] < all_objectives[best_idx * config.n_obj]) {
            best_idx = i;
        }
    }

    memcpy(best_solution, &all_solutions[best_idx * config.n_dim], config.n_dim * sizeof(double));
    memcpy(best_objectives, &all_objectives[best_idx * config.n_obj], config.n_obj * sizeof(double));
    *num_evaluations = total_evals;

    return 0;
}

int moecam_run_random_start(int config_id, double* best_solution, double* best_objectives,
                           double* all_solutions, double* all_objectives, int* num_evaluations) {
    // Simplified implementation - just random sampling
    if (config_id >= configs.size()) return -1;

    moecam_config_t& config = configs[config_id];
    std::uniform_real_distribution<> dis(0.0, 1.0);

    double best_val = 1e100;

    for (int i = 0; i < config.max_iterations; i++) {
        std::vector<double> solution(config.n_dim);
        std::vector<double> obj(config.n_obj);

        // Generate random solution
        for (int j = 0; j < config.n_dim; j++) {
            solution[j] = config.lower_bounds[j] +
                         dis(gen) * (config.upper_bounds[j] - config.lower_bounds[j]);
        }

        // Evaluate
        config.objective_func(solution.data(), config.n_dim, obj.data(), config.n_obj);

        // Store
        memcpy(&all_solutions[i * config.n_dim], solution.data(), config.n_dim * sizeof(double));
        memcpy(&all_objectives[i * config.n_obj], obj.data(), config.n_obj * sizeof(double));

        // Update best
        if (obj[0] < best_val) {
            best_val = obj[0];
            memcpy(best_solution, solution.data(), config.n_dim * sizeof(double));
            memcpy(best_objectives, obj.data(), config.n_obj * sizeof(double));
        }
    }

    *num_evaluations = config.max_iterations;
    return 0;
}

int moecam_run_direct(int config_id, double* best_solution, double* best_objective, int* num_evaluations) {
    // Simplified DIRECT-like algorithm
    if (config_id >= configs.size()) return -1;

    moecam_config_t& config = configs[config_id];
    std::uniform_real_distribution<> dis(0.0, 1.0);

    double best_val = 1e100;
    *num_evaluations = 0;

    // Simple grid search with refinement
    int grid_size = 10;
    for (int iter = 0; iter < 3 && *num_evaluations < config.max_iterations; iter++) {
        for (int i = 0; i < grid_size && *num_evaluations < config.max_iterations; i++) {
            std::vector<double> solution(config.n_dim);

            // Generate grid point with some randomness
            for (int j = 0; j < config.n_dim; j++) {
                double t = (double)i / (grid_size - 1) + 0.1 * dis(gen);
                solution[j] = config.lower_bounds[j] +
                             t * (config.upper_bounds[j] - config.lower_bounds[j]);
            }

            // Evaluate (single objective for DIRECT)
            double obj;
            config.objective_func(solution.data(), config.n_dim, &obj, 1);
            (*num_evaluations)++;

            // Update best
            if (obj < best_val) {
                best_val = obj;
                memcpy(best_solution, solution.data(), config.n_dim * sizeof(double));
                *best_objective = obj;
            }
        }
        grid_size *= 2; // Refine grid
    }

    return 0;
}

void moecam_cleanup(int config_id) {
    // Simple cleanup - in real implementation would free allocated memory
}
