import pandas as pd
import numpy as np
import random

from pathlib import Path
from itertools import combinations
from collections import defaultdict
from functools import lru_cache


def calc_stats(part2pop, parts):
    """counts mean amount of wins for a real population"""
    for part in parts:
        populations = sorted(part2pop[part])
        d = {pop: [] for pop in populations}
        d['part'] = []
        file_path = f'../features_part{part}.csv'
        df = pd.read_csv(file_path, index_col='SAMPLE')
        d['part'].append(part)
        win_counts = pd.Series(np.zeros(len(df)), index=df.index.values)

        for index, row in df.iterrows():
            n = 0
            pop = row['REAL']
            for i, value in row.items():
                if i.startswith(pop) and value >= 0.5:
                    n += 1
                elif i.endswith(pop) and value <= 0.5:
                    n += 1
            win_counts[index] = n / 4

        df['WINS'] = win_counts
        for pop in populations:
            pop_df = df[df['REAL'] == pop][['WINS']]
            d[pop].append(pop_df['WINS'].mean())

    res = pd.DataFrame(d)
    return res


class BaseAnalyser:

    name = 'BaseAnalyser'

    def __init__(self, populations, part):
        self.populations = sorted(populations)
        self.part = part
        self.comb = list(combinations(self.populations, 2))

    @property
    def _path(self):
        return f'../features_part{self.part}.csv'

    def _score(self, score):
        """
        Return a tuple of scores for population1 and population2
        if in mutual comparison population1 got value score
        """
        if score > 0.5:
            return 1, 0
        else:
            return 0, 1

    def winner_for_sample(self, row):
        scores = {pop: 0 for pop in self.populations}
        for a, b in self.comb:
            column = f'{a}_{b}'
            s1, s2 = self._score(row[column])
            scores[a] += s1
            scores[b] += s2
        winner = sorted(self.populations, key=lambda x: scores[x])[len(self.populations) - 1]
        return winner

    @lru_cache(1)
    def calc(self):
        df = pd.read_csv(self._path, index_col='SAMPLE')
        s = pd.Series(['' for _ in range(len(df))], index=df.index.values)

        for index, row in df.iterrows():
            s[index] = self.winner_for_sample(row)

        df['WINNER'] = s
        df['MATCH'] = df['REAL'] == df['WINNER']
        return df[['REAL', 'WINNER', 'MATCH']]

    def n_correctly_classified(self):
        df = self.calc()
        return len(df[df['MATCH']])

    def correct_classifications(self):
        d = {}
        for pop in self.populations:
            df = self.calc()
            d[pop] = len(df[(df['REAL'] == pop) & (df['WINNER'] == pop)]) / (1 + len(df[df['WINNER'] == pop]))
        return d

    def n_pop(self):
        df = self.calc()
        d = {pop: len(df[df['WINNER'] == pop]) for pop in self.populations}
        return d


class WinsCountAnalyser(BaseAnalyser):
    name = 'WinsCountAnalyser'


class WinsCountAnalyserWithThreshold(BaseAnalyser):
    name = 'WinsCountAnalyserWithThreshold'
    threshold = 0.3

    def _score(self, score):
        if score > 1 - self.threshold:
            return 1, 0
        elif score < self.threshold:
            return 0, 1
        else:
            return 0.5, 0.5


class ScoreAnalyserWithThreshold(BaseAnalyser):
    name = 'ScoreAnalyserWithThreshold'
    threshold = 0.2

    def _score(self, score):
        if score > 1 - self.threshold:
            return score, 0
        elif score < self.threshold:
            return 0, score
        else:
            return 0.5, 0.5


class ScoreAnalyser(BaseAnalyser):
    name = 'ScoreAnalyser'

    def _score(self, score):
        return score, 1 - score


class BinaryTreeAnalyser(BaseAnalyser):

    name = 'BinaryTreeAnalyser'

    def vs(self, pair, row):
        pair = sorted(pair)
        a, b = pair
        if a == b:
            return a
        else:
            col = "_".join(pair)
            return a if row[col] > 0.5 else b

    def tree_compare(self, pops, row):
        if len(pops) == 2:
            return self.vs(pops, row)

        start = pops[:-2]
        prelast = [pops[-2]]
        last = [pops[-1]]
        return self.vs((self.tree_compare(start + prelast, row), self.tree_compare(start + last, row)), row)

    def winner_for_sample(self, row):

        winner = self.tree_compare(self.populations, row)
        return winner


class SymmetricalBinaryTreeAnalyser(BinaryTreeAnalyser):

    name = 'SymmetricalBinaryTreeAnalyser'

    def tree_compare(self, pops, row):
        if len(pops) == 1:
            return pops[0]

        elif len(pops) == 2:
            return self.vs(pops, row)

        odd = pops[::2]
        even = pops[1::2]
        return self.vs((self.tree_compare(odd, row), self.tree_compare(even, row)), row)


res = defaultdict(list)
parts = (102, 207, 312, 413)

part_to_pops = {
    102: ('ESN', 'GBR', 'JPT', 'LWK', 'PEL'),
    207: ('CEU', 'ESN', 'GBR', 'MSL', 'PEL'),
    312: ('FIN', 'GWD', 'LWK', 'MXL', 'PEL'),
    413: ('GBR', 'KHV', 'MSL', 'PEL', 'YRI'),
}

stats = calc_stats(part_to_pops, parts)
print(stats)


for p in parts:
    # pops = part_to_pops[p]
    pops = sorted([
        'GBR', 'FIN', 'CHS', 'PUR', 'CDX', 'CLM', 'IBS', 'PEL', 'PJL', 'KHV', 'ACB', 'GWD', 'ESN',
        'BEB', 'MSL', 'STU', 'ITU', 'CEU', 'YRI', 'CHB', 'JPT', 'LWK', 'ASW', 'MXL', 'TSI', 'GIH'
    ])
    print('ON PART: ', p, 'with populations: ', pops)

    for Analyser in (WinsCountAnalyserWithThreshold, WinsCountAnalyser,
                     ScoreAnalyserWithThreshold, ScoreAnalyser,
                     SymmetricalBinaryTreeAnalyser):  # BinaryTreeAnalyser,

        analyser = Analyser(pops, p)
        t = analyser.calc()
        output_dir = Path(f'predictions/part{p}/')
        output_dir.mkdir(parents=True, exist_ok=True)

        output_file = output_dir / f'{analyser.name}.csv'
        t.to_csv(str(output_file))
        print(analyser.correct_classifications())
        res[analyser.name].append(analyser.n_correctly_classified())

print(res)
for key, val in res.items():
    print(key, sum(val) / 495 / len(parts))
df = pd.DataFrame(res)
df.to_csv('analyser_comparison.csv')






