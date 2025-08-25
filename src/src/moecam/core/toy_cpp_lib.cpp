
#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <cstring>

// C++ class for Multi-Objective Optimization
class MOOptimizer {
private:
    int n_dim;
    int n_obj;
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    std::mt19937 rng;
    std::uniform_real_distribution<double> dist;
    
public:
    MOOptimizer(int dimensions, int objectives, double* lb, double* ub) 
        : n_dim(dimensions), n_obj(objectives), rng(42), dist(0.0, 1.0) {
        lower_bounds.assign(lb, lb + dimensions);
        upper_bounds.assign(ub, ub + dimensions);
    }
    
    void generate_random_solution(double* solution) {
        for (int i = 0; i < n_dim; ++i) {
            double rand_val = dist(rng);
            solution[i] = lower_bounds[i] + rand_val * (upper_bounds[i] - lower_bounds[i]);
        }
    }
    
    void crossover(const double* parent1, const double* parent2, double* child1, double* child2) {
        // Simple uniform crossover
        for (int i = 0; i < n_dim; ++i) {
            if (dist(rng) < 0.5) {
                child1[i] = parent1[i];
                child2[i] = parent2[i];
            } else {
                child1[i] = parent2[i];
                child2[i] = parent1[i];
            }
        }
    }
    
    void mutate(double* individual, double mutation_rate) {
        for (int i = 0; i < n_dim; ++i) {
            if (dist(rng) < mutation_rate) {
                double rand_val = dist(rng);
                individual[i] = lower_bounds[i] + rand_val * (upper_bounds[i] - lower_bounds[i]);
            }
        }
    }
    
    int get_dimensions() const { return n_dim; }
    int get_objectives() const { return n_obj; }
};

// Global storage for optimizer instances
static std::vector<MOOptimizer*> optimizers;

extern "C" {
    // Legacy functions
    double add(double a, double b) {
        return a + b;
    }

    void print_message(const char* message) {
        std::cout << "C++ says: " << message << std::endl;
    }
    
    // New optimizer interface
    int create_optimizer(int n_dim, int n_obj, double* lower_bounds, double* upper_bounds) {
        MOOptimizer* opt = new MOOptimizer(n_dim, n_obj, lower_bounds, upper_bounds);
        optimizers.push_back(opt);
        return optimizers.size() - 1;  // Return handle/index
    }
    
    void destroy_optimizer(int handle) {
        if (handle >= 0 && handle < optimizers.size() && optimizers[handle]) {
            delete optimizers[handle];
            optimizers[handle] = nullptr;
        }
    }
    
    void generate_random_solution(int handle, double* solution) {
        if (handle >= 0 && handle < optimizers.size() && optimizers[handle]) {
            optimizers[handle]->generate_random_solution(solution);
        }
    }
    
    void crossover_solutions(int handle, const double* parent1, const double* parent2, 
                           double* child1, double* child2) {
        if (handle >= 0 && handle < optimizers.size() && optimizers[handle]) {
            optimizers[handle]->crossover(parent1, parent2, child1, child2);
        }
    }
    
    void mutate_solution(int handle, double* individual, double mutation_rate) {
        if (handle >= 0 && handle < optimizers.size() && optimizers[handle]) {
            optimizers[handle]->mutate(individual, mutation_rate);
        }
    }
    
    int get_optimizer_dimensions(int handle) {
        if (handle >= 0 && handle < optimizers.size() && optimizers[handle]) {
            return optimizers[handle]->get_dimensions();
        }
        return 0;
    }
    
    int get_optimizer_objectives(int handle) {
        if (handle >= 0 && handle < optimizers.size() && optimizers[handle]) {
            return optimizers[handle]->get_objectives();
        }
        return 0;
    }
}


