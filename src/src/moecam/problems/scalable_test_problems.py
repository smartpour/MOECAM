
import numpy as np

class ScalableProblem:
    def __init__(self, n_dim, n_obj, problem_type="sphere"):
        self.n_dim = n_dim
        self.n_obj = n_obj
        self.problem_type = problem_type
        self.xl = np.zeros(n_dim)
        self.xu = np.ones(n_dim)

    def evaluate(self, x):
        if self.problem_type == "sphere":
            # Simple sphere function for demonstration
            f1 = np.sum((x - 0.5)**2)
            f2 = np.sum((x + 0.5)**2)
            return np.array([f1, f2])
        elif self.problem_type == "linear":
            # Simple linear function
            f1 = np.sum(x)
            f2 = np.sum(1 - x)
            return np.array([f1, f2])
        else:
            raise ValueError("Unknown problem type")


