import sys
import time
import allel
import pathlib
import logging
import argparse
import tempfile
import subprocess
import pandas as pd
import numpy as np

from scipy.special import xlogy
from itertools import combinations
from functools import lru_cache, wraps
from queue import PriorityQueue
from operator import itemgetter
from random import choices
from math import log, factorial as fac


BASE_PATH = pathlib.Path(__file__).absolute().parents[1]

logger = logging.getLogger(__name__)
logger.setLevel('INFO')
logger.addHandler(logging.StreamHandler(sys.stdout))


POPULATIONS = [
    'GBR', 'FIN', 'CHS', 'PUR', 'CDX', 'CLM', 'IBS', 'PEL', 'PJL', 'KHV', 'ACB', 'GWD', 'ESN',
    'BEB', 'MSL', 'STU', 'ITU', 'CEU', 'YRI', 'CHB', 'JPT', 'LWK', 'ASW', 'MXL', 'TSI', 'GIH'
]


def run_time_logging(func):
    @wraps(func)
    def inner(*args, **kwargs):
        t = time.monotonic()
        result = func(*args, **kwargs)
        logger.info(f'Func {func.__name__} done in {time.monotonic() - t} seconds.')
        return result
    return inner


def xlogx(x):
    return xlogy(x, x)


def i4a(p_left, p_right):
    """Informativeness for assignment for 2 populations (Rosenberg et al.)"""
    p_avg = (p_left + p_right) / 2
    informativeness = -xlogx(p_avg) - xlogx(1 - p_avg) \
          + (xlogx(p_left) + xlogx(p_right) + xlogx(1 - p_left) + xlogx(1 - p_right)) / 2
    return informativeness


class BaseAnalyser:
    def __init__(self):
        self.base_path = BASE_PATH
        self.populations = sorted(POPULATIONS)
        self.pairs = list(combinations(self.populations, 2))

    @staticmethod
    def _check_local(path):
        return pathlib.Path(path).exists()

    def _path(self, file, absolute=True):
        return str(self.base_path / file) if absolute else file

    def _samples_path(self, pop='ALL'):
        """Path to samples from population"""
        return str(self.base_path / 'sample_lists' / f'{pop}.samples.list')

    def _prepare_samples(self):
        for pop in self.populations:
            if not self._check_local(self._samples_path(pop)):
                logger.info(f'Retrieving samples for {pop} to {self._samples_path(pop)}')
                command = ' '.join([
                    'grep', pop, self._samples_path('ALL'), '|', 'cut', '-f1', '>', self._samples_path(pop)
                ])
                subprocess.call(command, shell=True, stdout=subprocess.PIPE)

        logger.info(f'Samples for {self.populations} are ready.')

    def samples_df(self, pop='ALL'):
        header = 0 if pop == 'ALL' else None
        index = 'sample' if pop == 'ALL' else 0
        df = pd.read_csv(self._samples_path(pop), sep='\t', header=header, index_col=index)
        return df[df.iloc[:, 0].isin(self.populations)]

    def samples(self, pop='ALL'):
        df = self.samples_df(pop)
        return df.index.values

    def random_samples(self, k=20):
        return choices(self.samples(), k=k)


class PartialInfoCalc(BaseAnalyser):

    def __init__(self, chromosome, recombination, n_markers=15, directory=None):
        super().__init__()
        self.chr = chromosome
        self.n = n_markers
        self._init_dir(directory, recombination)
        self._prepare_samples()

    def _init_dir(self, directory, recombination):
        self.directory = directory or f'chromosome{self.chr}_data'
        self.data_path = self.base_path / self.directory
        self._chr_path = str(self.data_path / f'chr{self.chr}_full.vcf.gz')
        self._rec_path = str(self.data_path / recombination)
        self._groups_path = str(self.data_path / 'groups.tsv')

        freq_path = self.data_path / 'frequencies'
        marker_path = self.data_path / 'informativeness'
        genotype_path = self.data_path / 'genotypes'

        for p in (freq_path, marker_path, genotype_path):
            p.mkdir(parents=True, exist_ok=True)

    def _prepare_samples(self):
        self._freq_query = ' %AF_'.join(['%POS'] + self.populations) + '\n'
        self._freq_names = ['POS'] + [f'AF_{pop}' for pop in self.populations]

        df = pd.read_csv(self._samples_path('ALL'), sep='\t', index_col='sample')
        df = df.apply(lambda row: ','.join((row['pop'], row['super_pop'])), axis=1)
        df.to_csv(self._groups_path, sep='\t', header=False)

    def _freq_path(self, part):
        return str(self.data_path / 'frequencies' / f'ALL.freq.part{part}')

    def _inf_path(self, part):
        return str(self.data_path / 'informativeness' / f'ALL.inf.part{part}.csv')

    def _markers_path(self, part):
        return str(self.data_path / 'informativeness' / f'markers.part{part}.csv')

    def _sample_genotype_path(self, sample, part):
        return str(self.data_path / 'genotypes' / f'{sample}.part{part}.vcf.gz')

    @property
    @lru_cache(1)
    def segments(self):
        logger.info(f'Getting segments for chromosome{self.chr} from {self._rec_path}')
        segment_df = pd.read_csv(self._rec_path, sep=' ')
        return segment_df['POS'].values

    def _query(self, part):
        start = self.segments[part]
        end = self.segments[part + 1]
        return f'chr{self.chr}:{start}-{end}'

    @run_time_logging
    def _calculate_frequency(self, part):
        """Calculate frequencies for biallelic SNPs on the specified part of the genome."""
        cmd = f'tabix -h {self._chr_path} {self._query(part)} ' \
            f'| bcftools +fill-tags -Ou -- -S {self._groups_path} -t AF' \
            f'| bcftools query -f "{self._freq_query}" > {self._freq_path(part)}'
        subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)

    def frequency(self, part):
        if not self._check_local(self._freq_path(part)):
            self._calculate_frequency(part)
        return pd.read_csv(self._freq_path(part), sep=' ', names=self._freq_names, index_col=0)

    @run_time_logging
    def _calculate_informativeness(self, part):
        freq_df = self.frequency(part)
        informativeness = {}
        for pop_left, pop_right in self.pairs:
            col_name = '_'.join((pop_left, pop_right))
            informativeness[col_name] = i4a(freq_df[f'AF_{pop_left}'].values,
                                            freq_df[f'AF_{pop_right}'].values)
        df = pd.DataFrame(informativeness, index=freq_df.index)
        df.to_csv(self._inf_path(part))

    def informativeness(self, part):
        if not self._check_local(self._inf_path(part)):
            self._calculate_informativeness(part)
        return pd.read_csv(self._inf_path(part), index_col=0)

    @run_time_logging
    def _calculate_markers(self, part):
        inf_df = self.informativeness(part)
        with open(self._markers_path(part), 'w') as f:
            for pop_left, pop_right in self.pairs:
                col_name = '_'.join((pop_left, pop_right))
                best_markers = inf_df[col_name].nlargest(n=self.n)
                idx_inf_pairs = zip(best_markers.index.values, best_markers.values)
                s = col_name + ',' + ','.join([str(item) for sublist in idx_inf_pairs for item in sublist])
                f.write(s + '\n')

    def markers(self, part):
        if not self._check_local(self._markers_path(part)):
            self._calculate_markers(part)
        return pd.read_csv(self._markers_path(part), header=None, index_col=0)

    # def sample_genotype_array(self, sample, part):
    #     """Get a genotype array for the specified individual"""
    #     file = self._sample_genotype_path(sample, part)
    #
    #     if not self._check_local(file):
    #         cmd = f'bcftools view -s {sample} -v snps -m2 -M2 -Oz -o {file} {self._chr_path} {self._query(part)}'
    #         subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)
    #
    #     gt = allel.read_vcf(file, fields=['GT', 'POS'])
    #     return pd.Series(allel.GenotypeArray(gt['calldata/GT'])[:, 0], index=gt['variants/POS'])
    #
    # def closest_populations(self, sample, part, n=10):
    #     """Return an ordered list of Populations from best to worst"""
    #
    #     if n >= len(self.populations):
    #         logger.info(f'Closest populations for {sample}: {self.populations}')
    #         return self.populations
    #
    #     freq = self.frequency(part)
    #     freq[sample] = self.sample_genotype_array(sample, part)
    #
    #     likeliness = {pop: 0 for pop in self.populations}
    #     for pos, row in freq.iterrows():
    #         left, right = row[sample]
    #         for pop in self.populations:
    #             likeliness[pop] += row[pop] if left == 0 else 1 - row[pop]
    #             likeliness[pop] += row[pop] if right == 0 else 1 - row[pop]
    #
    #     best = [k for k in sorted(likeliness, key=likeliness.get, reverse=True)]
    #     logger.info(f'Closest populations for {sample}: {best[:n]}')
    #     return best[:n]
    #
    # def real_population(self, sample):
    #     """Return a real population for the specified individual"""
    #     df = self.samples_df()
    #     return df.loc[sample, 'pop']
    #
    # def generate_features_as_csv(self, part):
    #     t = time.monotonic()
    #     logger.info(f'Generating features. Start time {t}')
    #     markers_df = self.markers_df(part)
    #     freq = self.frequency(part)
    #     for sample in self.samples():
    #         logger.info(f'Running feature calculator for {sample}')
    #         genotype = self.sample_genotype_array(sample, part)
    #         populations = sorted(self.closest_populations(sample, part))
    #         real = self.real_population(sample)
    #         scores = []
    #         for first, second in combinations(populations, 2):
    #             first_score, second_score = 0, 0
    #             line = markers_df.xs((first, second)).values
    #             for marker, inf in zip(line[1::2], line[::2]):
    #                 for allele in genotype[marker]:
    #                     first_freq, second_freq = freq.loc[int(marker), first], freq.loc[int(marker), second]
    #                     if allele == 0:
    #                         if first_freq > second_freq:
    #                             first_score += inf
    #                         else:
    #                             second_score += inf
    #                     else:
    #                         if first_freq > second_freq:
    #                             second_score += inf
    #                         else:
    #                             first_score += inf
    #             scores.append(f"{first_score / (first_score + second_score):.5f}")
    #         logger.info(f'Current execution time: {time.monotonic() - t}')
    #         with open(f'features_{self._query(part)}.csv', 'a') as f:
    #             f.write(f'{sample} {real} {" ".join(populations)} {" ".join(scores)}\n')


if __name__ == '__main__':
    pic = PartialInfoCalc(chromosome=18, recombination='recombination_spots_18.tsv',
                          n_markers=15, directory='chromosome18_data')

    print(pic.markers(0))
