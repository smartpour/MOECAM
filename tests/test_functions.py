
import numpy as np

class ZDT1:
    def __init__(self, n_dim=30):
        self.n_dim = n_dim
        self.n_obj = 2
        self.xl = np.zeros(n_dim)
        self.xu = np.ones(n_dim)

    def evaluate(self, x):
        f1 = x[0]
        g = 1 + 9 * np.sum(x[1:]) / (self.n_dim - 1)
        h = 1 - np.sqrt(f1 / g) - (f1 / g) * np.sin(10 * np.pi * f1)
        f2 = g * h
        return np.array([f1, f2])

class ZDT2:
    def __init__(self, n_dim=30):
        self.n_dim = n_dim
        self.n_obj = 2
        self.xl = np.zeros(n_dim)
        self.xu = np.ones(n_dim)

    def evaluate(self, x):
        f1 = x[0]
        g = 1 + 9 * np.sum(x[1:]) / (self.n_dim - 1)
        h = 1 - (f1 / g)**2
        f2 = g * h
        return np.array([f1, f2])

class ZDT3:
    def __init__(self, n_dim=30):
        self.n_dim = n_dim
        self.n_obj = 2
        self.xl = np.zeros(n_dim)
        self.xu = np.ones(n_dim)

    def evaluate(self, x):
        f1 = x[0]
        g = 1 + 9 * np.sum(x[1:]) / (self.n_dim - 1)
        h = 1 - np.sqrt(f1 / g) - (f1 / g) * np.sin(10 * np.pi * f1)
        f2 = g * h
        return np.array([f1, f2])

class DTLZ1:
    def __init__(self, n_dim=7, n_obj=3):
        self.n_dim = n_dim
        self.n_obj = n_obj
        self.xl = np.zeros(n_dim)
        self.xu = np.ones(n_dim)

    def evaluate(self, x):
        g = 100 * (self.n_dim - self.n_obj + 1 + np.sum((x[self.n_obj-1:] - 0.5)**2 - np.cos(20 * np.pi * (x[self.n_obj-1:] - 0.5))))
        f = np.zeros(self.n_obj)
        for i in range(self.n_obj):
            f[i] = 0.5 * (1 + g)
            for j in range(self.n_obj - 1 - i):
                f[i] *= x[j]
            if i < self.n_obj - 1:
                f[i] *= (1 - x[self.n_obj - 1 - i])
        return f

class DTLZ2:
    def __init__(self, n_dim=12, n_obj=3):
        self.n_dim = n_dim
        self.n_obj = n_obj
        self.xl = np.zeros(n_dim)
        self.xu = np.ones(n_dim)

    def evaluate(self, x):
        g = np.sum((x[self.n_obj-1:] - 0.5)**2)
        f = np.zeros(self.n_obj)
        for i in range(self.n_obj):
            f[i] = (1 + g)
            for j in range(self.n_obj - 1 - i):
                f[i] *= np.cos(x[j] * np.pi / 2)
            if i < self.n_obj - 1:
                f[i] *= np.sin(x[self.n_obj - 1 - i] * np.pi / 2)
        return f


