from functools import lru_cache
from itertools import combinations
import random


class Population:
    populations = ('EAS', 'EUR', 'AFR', 'AMR', 'SAS')
    subpopulations = ()

    def __init__(self, populations=None):
        self.all = populations or self.populations

    def __len__(self):
        return len(self.all)

    @property
    @lru_cache(1)
    def pair(self):
        """Sample pair of populations"""
        return random.sample(self.all, 2)

    @property
    @lru_cache(1)
    def pairs(self):
        return combinations(self.all, 2)

    @lru_cache(3)
    def combinations(self, n=2):
        return combinations(self.all, n)


