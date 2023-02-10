import numpy as np
import pandas as pd
from numpy.random import choice
from scipy.stats import beta, norm

class DirichletProcessSample():
    def __init__(self, base_measure, alpha):
        self.base_measure = base_measure
        self.alpha = alpha

        self.cache = []
        self.weights = []
        self.total_stick_used = 0

    def __call__(self):
        remaining = 1.0 - self.total_stick_used
        i = DirichletProcessSample.roll_die(self.weights + [remaining])
        if i is not None and i < len(self.weights):
            return self.cache[i]
        else:
            stick_piece = beta(1, self.alpha).rvs() * remaining
            self.total_stick_used += stick_piece
            self.weights.append(stick_piece)
            new_value = self.base_measure()
            self.cache.append(new_value)
            return new_value

    @staticmethod
    def roll_die(weights):
        if weights:
            return choice(range(len(weights)), p=weights)
        else:
            return None

if __name__ == '__main__':

    base_measure = lambda: norm().rvs()
    n_samples = 10000
    samples = {}
    for alpha in [1, 10, 100, 1000]:
        dirichlet_norm = DirichletProcessSample(base_measure=base_measure, alpha=alpha)
        samples["Alpha: %s" % alpha] = [dirichlet_norm() for _ in range(n_samples)]

    _ = pd.DataFrame(samples).hist()



