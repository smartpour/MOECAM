
import time
import numpy as np

def pareto_front(points):
    is_pareto = np.ones(points.shape[0], dtype=bool)
    for i, p1 in enumerate(points):
        for j, p2 in enumerate(points):
            if i != j and np.all(p2 <= p1) and np.any(p2 < p1):
                is_pareto[i] = False
                break
    return points[is_pareto]

def hypervolume(pareto_front_points, reference_point):
    if not pareto_front_points.size:
        return 0.0

    sorted_points = pareto_front_points[pareto_front_points[:, 0].argsort()]

    hv = 0.0
    if sorted_points.shape[1] == 2:
        current_max_y = reference_point[1]
        for i in range(sorted_points.shape[0] - 1, -1, -1):
            p = sorted_points[i]
            if p[1] < current_max_y:
                width = reference_point[0] - p[0]
                height = current_max_y - p[1]
                hv += width * height
                current_max_y = p[1]
    else:
        print("Hypervolume calculation for >2 dimensions is complex and not fully implemented in this example.")
        return -1.0

    return hv

class EvaluationCounter:
    def __init__(self):
        self.count = 0
        self.start_time = None
        self.end_time = None

    def start(self):
        self.count = 0
        self.start_time = time.time()

    def increment(self):
        self.count += 1

    def stop(self):
        self.end_time = time.time()

    def get_count(self):
        return self.count

    def get_time(self):
        if self.start_time is None or self.end_time is None:
            return 0.0
        return self.end_time - self.start_time

# Example usage in an algorithm:
# counter = EvaluationCounter()
# counter.start()
# for _ in range(num_evaluations):
#     problem.evaluate(self.x)
#     counter.increment()
# counter.stop()
# print(f"Evaluations: {counter.get_count()}, Time: {counter.get_time():.2f}s")


