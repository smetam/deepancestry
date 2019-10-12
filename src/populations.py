from functools import lru_cache
from itertools import combinations
import random
import pathlib


BASE_PATH = pathlib.Path(__file__).parents[1]


class Population:
    superpopulations = ('EAS', 'EUR', 'AFR', 'AMR', 'SAS')
    populations = [
        'GBR', 'FIN', 'CHS', 'PUR', 'CDX', 'CLM', 'IBS', 'PEL', 'PJL', 'KHV', 'ACB', 'GWD', 'ESN',
        'BEB', 'MSL', 'STU', 'ITU', 'CEU', 'YRI', 'CHB', 'JPT', 'LWK', 'ASW', 'MXL', 'TSI', 'GIH'
    ]
    population_map = {
        'EUR': ['GBR', 'FIN', 'IBS', 'CEU', 'TSI'],
        'EAS': ['CHS', 'CDX', 'KHV', 'CHB', 'JPT'],
        'AMR': ['PUR', 'CLM', 'PEL', 'MXL'],
        'SAS': ['PJL', 'BEB', 'STU', 'ITU', 'GIH'],
        'AFR': ['ACB', 'GWD', 'ESN', 'MSL', 'YRI', 'LWK', 'ASW']
    }

    def __init__(self, populations=None):
        self.all = populations or self.superpopulations

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


