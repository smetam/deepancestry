import sys
import time
import pathlib
import logging
import argparse
import subprocess
import pandas as pd
import numpy as np

from itertools import combinations
from functools import lru_cache, wraps
from random import choices
from scipy.special import xlogy


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


def lazy_property(f):
    return property(lru_cache(1)(f))


def xlogx(x):
    return xlogy(x, x)


def i4a(p_left, p_right):
    """Informativeness for assignment for 2 populations (Rosenberg et al.)"""
    p_avg = (p_left + p_right) / 2
    informativeness = -xlogx(p_avg) - xlogx(1 - p_avg) \
        + (xlogx(p_left) + xlogx(p_right) + xlogx(1 - p_left) + xlogx(1 - p_right)) / 2
    return informativeness


def calc_scores(pl, pr, allele1, allele2):
    score = 0
    for allele in (allele1, allele2):
        if allele == 0:
            score += (pr - pl) / (2 - pl - pr)
        if allele == 1:
            score += (pl - pr) / (pl + pr)


def calc_score_array(pl, pr, inf, allele1, allele2):
    p_sum = pl + pr
    p_diff = pl - pr
    scores = np.vstack((-p_diff / (2 - p_sum), p_diff / p_sum)).T
    score = scores[np.arange(pl.size), allele1] + scores[np.arange(pl.size), allele2]
    return np.sum(score * inf)


class PartialInfoCalc:

    def __init__(self, chromosome, recombination_file, groups_file, n_markers=15, directory=None):
        self.chr = chromosome
        self.n = n_markers
        self.base_path = BASE_PATH

        self._init_dir(directory, recombination_file, groups_file)
        self._prepare_samples()

    def _init_dir(self, directory, recombination_file, groups_file):
        self.directory = directory or f'chromosome{self.chr}_data'
        self.data_path = self.base_path / self.directory
        self._chr_path = str(self.data_path / f'chr{self.chr}_full.vcf.gz')
        self._rec_path = str(self.data_path / recombination_file)
        self._groups_path = str(self.data_path / groups_file)

        for folder in ('frequencies', 'markers', 'genotypes', 'informativeness', 'features'):
            (self.data_path / folder).mkdir(parents=True, exist_ok=True)

    def _prepare_samples(self):
        self.populations = sorted(self._samples_df['GROUP'].unique())
        self.pairs = list(combinations(self.populations, 2))

        self._freq_query = ' %AF_'.join(['%POS'] + self.populations) + '\n'
        self._freq_names = ['POS'] + self.populations

    @lazy_property
    def _samples_df(self):
        return pd.read_csv(self._groups_path, sep='\t', header=None, index_col=0, names=['SAMPLE', 'GROUP'])

    def samples(self, population_list=None):
        if not population_list:
            return self._samples_df.index.values
        return self._samples_df[self._samples_df['GROUP'].isin(population_list)].index.values

    def true_group(self, sample):
        return self._samples_df.at[sample, 'GROUP']

    @staticmethod
    def _check_local(path):
        return pathlib.Path(path).exists()

    def _freq_path(self, part):
        return str(self.data_path / 'frequencies' / f'ALL.freq.part{part}')

    def _inf_path(self, part):
        return str(self.data_path / 'informativeness' / f'ALL.inf.part{part}.csv')

    def _markers_path(self, part):
        return str(self.data_path / 'markers' / f'ALL.markers.part{part}.csv')

    def _features_path(self, part, pop_list=None):
        prefix = "_".join(sorted(pop_list)) if pop_list else 'ALL'
        return str(self.data_path / 'features' / f'{prefix}.features.part{part}')

    def _sample_genotype_path(self, sample, part):
        return str(self.data_path / 'genotypes' / f'{sample}.genotype.part{part}')

    @lazy_property
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
        """Calculate frequencies for biallelic SNPs for the specified part of the genome."""
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
        """Calculate informativeness for assignment for each snp position."""
        freq_df = self.frequency(part)
        informativeness = {}
        for pop_left, pop_right in self.pairs:
            col_name = '_'.join((pop_left, pop_right))
            informativeness[col_name] = i4a(freq_df[pop_left].values, freq_df[pop_right].values)
        df = pd.DataFrame(informativeness, index=freq_df.index)
        df.to_csv(self._inf_path(part))

    def informativeness(self, part):
        if not self._check_local(self._inf_path(part)):
            self._calculate_informativeness(part)
        return pd.read_csv(self._inf_path(part), index_col=0)

    @run_time_logging
    def _calculate_markers(self, part):
        """Find best markers for each pair of populations."""
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

    def _calculate_genotype_array(self, sample, part):
        cmd = f'bcftools view -s {sample} -Ou {self._chr_path} {self._query(part)} ' \
            f'| bcftools query -f "%POS|[%GT]\n" > {self._sample_genotype_path(sample, part)}'
        subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)

    def sample_genotype_array(self, sample, part):
        if not self._check_local(self._sample_genotype_path(sample, part)):
            self._calculate_genotype_array(sample, part)
        return pd.read_csv(self._sample_genotype_path(sample, part), header=None, sep='|', index_col=0)

    # TODO: rewrite file
    def _write_scores(self, sample, part, scores, pop_list):
        real = self.true_group(sample)
        with open(self._features_path(part, pop_list), 'a') as f:
            f.write(f'{sample} {real} {" ".join(scores)}\n')

    @run_time_logging
    def _calculate_features(self, part, pop_list=None):
        markers_df = self.markers(part)
        freq_df = self.frequency(part)
        t = time.monotonic()
        pairs = list(combinations(sorted(pop_list), 2)) if pop_list else self.pairs
        for sample in self.samples(pop_list):
            logger.info(f'Calculating scores for {sample}, {time.monotonic()-t} since start.')
            genotype_df = self.sample_genotype_array(sample, part)
            scores = []
            for pop_left, pop_right in pairs:
                markers_list = markers_df.loc[f'{pop_left}_{pop_right}'].values
                idx = markers_list[::2]
                inf = markers_list[1::2]
                allele1 = genotype_df.loc[idx][1].values
                allele2 = genotype_df.loc[idx][2].values
                p_left, p_right = freq_df[pop_left].loc[idx], freq_df[pop_right].loc[idx]
                score = calc_score_array(p_left, p_right, inf, allele1, allele2)
                scores.append(f"{score:.5f}")

            self._write_scores(sample, part, scores, pop_list)

    def features(self, part):
        if not self._check_local(self._features_path(part)):
            self._calculate_features(part)
        return pd.read_csv(self._features_path(part), index_col=0)

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


if __name__ == '__main__':
    pic = PartialInfoCalc(chromosome=18, recombination_file='recombination_spots_18.tsv',
                          groups_file='single_groups.tsv',
                          n_markers=15, directory='chromosome18_data')

    print(pic.features(0))
