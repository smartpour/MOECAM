/*
 * MOECAM Algorithms C++ Wrapper Library
 * ====================================
 *
 * C-compatible wrapper for MOECAM optimization algorithms.
 * Provides interface for ECAM, Random Start, and DIRECT algorithms.
 */

#include <map>
#include <vector>
#include <memory>
#include <iostream>
#include <algorithm>
#include <random>
#include <functional>
#include <cmath>

extern "C" {

// Callback type for objective function evaluation
typedef void (*objective_function_t)(double* x, int n_dim, double* f, int n_obj);

// MOECAM algorithm configuration
typedef struct {
    int algorithm_id;
    int n_dim;
    int n_obj;
    int max_iterations;
    double* lower_bounds;
    double* upper_bounds;
    objective_function_t objective_func;
} moecam_config_t;

// Internal configuration storage
struct MOECAMConfig {
    int algorithm_id;
    int n_dim;
    int n_obj;
    int max_iterations;
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    objective_function_t objective_func;

    // Result storage
    std::vector<std::vector<double>> all_solutions;
    std::vector<std::vector<double>> all_objectives;
    int function_evaluations;

    MOECAMConfig() : function_evaluations(0) {}
};

// Global configuration storage
static std::map<int, std::unique_ptr<MOECAMConfig>> g_configs;
static int g_next_id = 1;

// Random number generator
static std::random_device rd;
static std::mt19937 rng(rd());

// Helper function to evaluate objective
void evaluate_objective(MOECAMConfig* config, const std::vector<double>& x, std::vector<double>& f) {
    config->objective_func(const_cast<double*>(x.data()), config->n_dim, f.data(), config->n_obj);
    config->function_evaluations++;

    // Store for output
    config->all_solutions.push_back(x);
    config->all_objectives.push_back(f);
}

// ECAM Algorithm Implementation
class ECAMAlgorithm {
private:
    MOECAMConfig* config;
    std::vector<std::vector<double>> population;
    std::vector<std::vector<double>> objectives;

public:
    ECAMAlgorithm(MOECAMConfig* cfg) : config(cfg) {}

    void initialize_population() {
        population.clear();
        objectives.clear();

        // Use small initial population that fits within iteration limit
        int init_pop_size = std::min(5, config->max_iterations);
        if (init_pop_size < 1) init_pop_size = 1;

        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (int i = 0; i < init_pop_size && config->function_evaluations < config->max_iterations; i++) {
            std::vector<double> individual(config->n_dim);
            for (int j = 0; j < config->n_dim; j++) {
                double range = config->upper_bounds[j] - config->lower_bounds[j];
                individual[j] = config->lower_bounds[j] + dist(rng) * range;
            }

            std::vector<double> obj(config->n_obj);
            evaluate_objective(config, individual, obj);

            population.push_back(individual);
            objectives.push_back(obj);
        }
    }

    void evolve_population() {
        if (population.empty() || config->function_evaluations >= config->max_iterations) return;

        int pop_size = population.size();
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        // Generate only a few offspring to stay within iteration limit
        int max_offspring = config->max_iterations - config->function_evaluations;
        int num_offspring = std::min(2, max_offspring);

        for (int i = 0; i < num_offspring && config->function_evaluations < config->max_iterations; i++) {
            // Select two parents
            int p1 = rng() % pop_size;
            int p2 = rng() % pop_size;

            // Simple crossover
            std::vector<double> offspring(config->n_dim);
            double alpha = dist(rng);

            for (int j = 0; j < config->n_dim; j++) {
                offspring[j] = alpha * population[p1][j] + (1.0 - alpha) * population[p2][j];

                // Ensure bounds
                offspring[j] = std::max(config->lower_bounds[j],
                                      std::min(config->upper_bounds[j], offspring[j]));
            }

            // Evaluate offspring
            std::vector<double> obj(config->n_obj);
            evaluate_objective(config, offspring, obj);

            // Replace worst individual if offspring is better
            int worst_idx = 0;
            double worst_score = 0;
            for (int k = 0; k < pop_size; k++) {
                double score = 0;
                for (int obj_idx = 0; obj_idx < config->n_obj; obj_idx++) {
                    score += objectives[k][obj_idx];
                }
                if (k == 0 || score > worst_score) {
                    worst_score = score;
                    worst_idx = k;
                }
            }

            // Check if offspring dominates worst
            double offspring_score = 0;
            for (int obj_idx = 0; obj_idx < config->n_obj; obj_idx++) {
                offspring_score += obj[obj_idx];
            }

            if (offspring_score < worst_score) {
                population[worst_idx] = offspring;
                objectives[worst_idx] = obj;
            }
        }
    }

    void get_best(std::vector<double>& best_x, std::vector<double>& best_f) {
        if (population.empty()) {
            best_x.resize(config->n_dim, 0.0);
            best_f.resize(config->n_obj, 1e6);
            return;
        }

        // Find best individual (lowest sum of objectives)
        int best_idx = 0;
        double best_score = 0;

        for (int i = 0; i < population.size(); i++) {
            double score = 0;
            for (int obj = 0; obj < config->n_obj; obj++) {
                score += objectives[i][obj];
            }

            if (i == 0 || score < best_score) {
                best_score = score;
                best_idx = i;
            }
        }

        best_x = population[best_idx];
        best_f = objectives[best_idx];
    }
};// Random Start Algorithm Implementation
class RandomStartAlgorithm {
private:
    MOECAMConfig* config;
    std::vector<double> best_solution;
    std::vector<double> best_objectives;

public:
    RandomStartAlgorithm(MOECAMConfig* cfg) : config(cfg) {}

    void run() {
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        bool first = true;

        while (config->function_evaluations < config->max_iterations) {
            std::vector<double> x(config->n_dim);
            for (int j = 0; j < config->n_dim; j++) {
                double range = config->upper_bounds[j] - config->lower_bounds[j];
                x[j] = config->lower_bounds[j] + dist(rng) * range;
            }

            std::vector<double> f(config->n_obj);
            evaluate_objective(config, x, f);

            if (first || dominates(f, best_objectives)) {
                best_solution = x;
                best_objectives = f;
                first = false;
            }
        }
    }

    bool dominates(const std::vector<double>& a, const std::vector<double>& b) {
        bool at_least_one_better = false;
        for (int i = 0; i < config->n_obj; i++) {
            if (a[i] > b[i]) return false;
            if (a[i] < b[i]) at_least_one_better = true;
        }
        return at_least_one_better;
    }

    void get_best(std::vector<double>& best_x, std::vector<double>& best_f) {
        best_x = best_solution;
        best_f = best_objectives;
    }
};

// DIRECT Algorithm Implementation (simplified)
class DIRECTAlgorithm {
private:
    MOECAMConfig* config;
    std::vector<double> best_solution;
    double best_objective;

public:
    DIRECTAlgorithm(MOECAMConfig* cfg) : config(cfg), best_objective(1e30) {}

    void run() {
        // Simplified DIRECT: systematic grid search with refinement
        int levels = std::min(5, (int)std::sqrt(config->max_iterations));

        for (int level = 0; level < levels; level++) {
            int points_per_dim = (int)std::pow(2, level) + 1;
            search_level(points_per_dim);

            if (config->function_evaluations >= config->max_iterations) break;
        }
    }

    void search_level(int points_per_dim) {
        std::vector<int> indices(config->n_dim, 0);

        do {
            if (config->function_evaluations >= config->max_iterations) break;

            std::vector<double> x(config->n_dim);
            for (int j = 0; j < config->n_dim; j++) {
                double t = (double)indices[j] / (points_per_dim - 1);
                double range = config->upper_bounds[j] - config->lower_bounds[j];
                x[j] = config->lower_bounds[j] + t * range;
            }

            std::vector<double> f(1);  // Single objective
            evaluate_objective(config, x, f);

            if (f[0] < best_objective) {
                best_objective = f[0];
                best_solution = x;
            }

        } while (next_combination(indices, points_per_dim));
    }

    bool next_combination(std::vector<int>& indices, int max_val) {
        for (int i = 0; i < config->n_dim; i++) {
            if (indices[i] < max_val - 1) {
                indices[i]++;
                return true;
            }
            indices[i] = 0;
        }
        return false;
    }

    void get_best(std::vector<double>& best_x, double& best_f) {
        best_x = best_solution;
        best_f = best_objective;
    }
};

// C Interface Implementation
int moecam_init(moecam_config_t* config) {
    try {
        auto internal_config = std::make_unique<MOECAMConfig>();

        internal_config->algorithm_id = config->algorithm_id;
        internal_config->n_dim = config->n_dim;
        internal_config->n_obj = config->n_obj;
        internal_config->max_iterations = config->max_iterations;
        internal_config->objective_func = config->objective_func;

        // Copy bounds
        internal_config->lower_bounds.resize(config->n_dim);
        internal_config->upper_bounds.resize(config->n_dim);
        for (int i = 0; i < config->n_dim; i++) {
            internal_config->lower_bounds[i] = config->lower_bounds[i];
            internal_config->upper_bounds[i] = config->upper_bounds[i];
        }

        int id = g_next_id++;
        g_configs[id] = std::move(internal_config);

        return id;
    } catch (...) {
        return -1;
    }
}

int moecam_run_ecam(int config_id, double* best_solution, double* best_objectives,
                   double* all_solutions, double* all_objectives, int* num_evaluations) {
    try {
        auto it = g_configs.find(config_id);
        if (it == g_configs.end()) return -1;

        MOECAMConfig* config = it->second.get();
        config->all_solutions.clear();
        config->all_objectives.clear();
        config->function_evaluations = 0;

        ECAMAlgorithm ecam(config);
        ecam.initialize_population();

        // Simple evolution loop with strict iteration control
        while (config->function_evaluations < config->max_iterations) {
            int evals_before = config->function_evaluations;
            ecam.evolve_population();

            // Safety check: if no progress, break to avoid infinite loop
            if (config->function_evaluations == evals_before) {
                break;
            }
        }

        std::vector<double> best_x, best_f;
        ecam.get_best(best_x, best_f);

        // Copy results safely
        for (int i = 0; i < config->n_dim && i < best_x.size(); i++) {
            best_solution[i] = best_x[i];
        }
        for (int i = 0; i < config->n_obj && i < best_f.size(); i++) {
            best_objectives[i] = best_f[i];
        }

        // Copy all evaluations safely
        int actual_evals = std::min((int)config->all_solutions.size(), config->max_iterations);
        for (int i = 0; i < actual_evals; i++) {
            for (int j = 0; j < config->n_dim && j < config->all_solutions[i].size(); j++) {
                all_solutions[i * config->n_dim + j] = config->all_solutions[i][j];
            }
            for (int j = 0; j < config->n_obj && j < config->all_objectives[i].size(); j++) {
                all_objectives[i * config->n_obj + j] = config->all_objectives[i][j];
            }
        }

        *num_evaluations = actual_evals;
        return 0;

    } catch (...) {
        return -1;
    }
}int moecam_run_random_start(int config_id, double* best_solution, double* best_objectives,
                           double* all_solutions, double* all_objectives, int* num_evaluations) {
    try {
        auto it = g_configs.find(config_id);
        if (it == g_configs.end()) return -1;

        MOECAMConfig* config = it->second.get();
        config->all_solutions.clear();
        config->all_objectives.clear();
        config->function_evaluations = 0;

        RandomStartAlgorithm rs(config);
        rs.run();

        std::vector<double> best_x, best_f;
        rs.get_best(best_x, best_f);

        // Copy results
        for (int i = 0; i < config->n_dim; i++) {
            best_solution[i] = best_x[i];
        }
        for (int i = 0; i < config->n_obj; i++) {
            best_objectives[i] = best_f[i];
        }

        // Copy all evaluations
        for (int i = 0; i < config->all_solutions.size(); i++) {
            for (int j = 0; j < config->n_dim; j++) {
                all_solutions[i * config->n_dim + j] = config->all_solutions[i][j];
            }
            for (int j = 0; j < config->n_obj; j++) {
                all_objectives[i * config->n_obj + j] = config->all_objectives[i][j];
            }
        }

        *num_evaluations = config->function_evaluations;
        return 0;

    } catch (...) {
        return -1;
    }
}

int moecam_run_direct(int config_id, double* best_solution, double* best_objective, int* num_evaluations) {
    try {
        auto it = g_configs.find(config_id);
        if (it == g_configs.end()) return -1;

        MOECAMConfig* config = it->second.get();
        config->all_solutions.clear();
        config->all_objectives.clear();
        config->function_evaluations = 0;

        DIRECTAlgorithm direct(config);
        direct.run();

        std::vector<double> best_x;
        double best_f;
        direct.get_best(best_x, best_f);

        // Copy results
        for (int i = 0; i < config->n_dim; i++) {
            best_solution[i] = best_x[i];
        }
        *best_objective = best_f;
        *num_evaluations = config->function_evaluations;

        return 0;

    } catch (...) {
        return -1;
    }
}

void moecam_cleanup(int config_id) {
    auto it = g_configs.find(config_id);
    if (it != g_configs.end()) {
        g_configs.erase(it);
    }
}

} // extern "C"
