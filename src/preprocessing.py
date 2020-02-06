import vcf
import sys
import time
import json
import allel
import pathlib
import logging
import argparse
import tempfile
import subprocess
import pandas as pd

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


def x_logx(x):
    return 0 if x == 0 else x * log(x)


def i4a(p1, p2):
    """Calculate informativeness for assignment for 2 populations"""
    inf = x_logx(p1) + x_logx(1 - p1) + x_logx(p2) + x_logx(1 - p2)
    inf /= 2
    inf -= x_logx((p1 + p2) / 2) + x_logx((2 - p1 - p2) / 2)
    return inf


def make_header(pairs):
    return 'SAMPLE,REAL,' + ','.join(list(map(lambda x: "_".join(x), pairs)))


class BaseAnalyser:
    def __init__(self):
        self.base_path = BASE_PATH
        self.populations = sorted(POPULATIONS)
        self._prepare_samples()

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

    def load(self):
        # TODO: load requested region of chromosome
        pass


class PartialInfoCalc(BaseAnalyser):

    def __init__(self, chromosome=1, n_markers=15, part_size=1):
        super().__init__()
        self.chr = chromosome
        self.n = n_markers
        self.part_size = part_size
        self._init_dir()

    def _init_dir(self):
        freq_path = self.base_path / 'data' / f'part_size{self.part_size}' / 'frequencies'
        marker_path = self.base_path / 'data' / f'part_size{self.part_size}' / 'informativeness'

        for p in (freq_path, marker_path):
            p.mkdir(parents=True, exist_ok=True)

    @property
    @lru_cache(1)
    def segments(self):
        segment_path = str(BASE_PATH / 'data' / 'plink' / f'plink.chr{self.chr}.GRCh37.map')
        logger.info(f'Getting segments for chromosome{self.chr} from {segment_path}')
        segment_df = pd.read_csv(segment_path, header=None, sep='\t')
        segments = list(segment_df[3])
        return segments

    @property
    def _chr_path(self):
        """Path to vcf with chr"""
        chr_filename = f'ALL.chr{self.chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
        return str(self.base_path / 'data' / 'chromosomes' / chr_filename)

    def _freq_path(self, part):
        return str(self.base_path / 'data' / f'part_size{self.part_size}' / 'frequencies' / f'ALL.frq.part{part}')

    def _query(self, part):
        start = self.segments[self.part_size * part]
        end = self.segments[self.part_size * (part + 1)]
        return f'{self.chr}:{start}-{end}'

    def _calculate_freq(self, part: int, with_ref=False):
        """Calculate frequencies for biallelic SNPs on the specified part of the genome."""
        start_time = time.monotonic()
        frequency_df = pd.DataFrame(columns=self.populations)
        with tempfile.TemporaryDirectory(dir=BASE_PATH) as tmp:
            for pop in self.populations:
                file = str(pathlib.Path(tmp) / f'{pop}.chr{self.chr}')
                cmd1 = f'bcftools view -S {self._samples_path(pop)} -v snps -m2 -M2 ' \
                    f'-Oz -o {file}.vcf.gz {self._chr_path} {self._query(part)}'
                subprocess.call(cmd1, shell=True, stdout=subprocess.PIPE)

                cmd2 = f'vcftools --gzvcf {file}.vcf.gz --freq2 --stdout | cut -f "2 5" > {file}.part{part}'
                subprocess.call(cmd2, shell=True, stdout=subprocess.PIPE)

                df = pd.read_csv(f'{file}.part{part}', index_col='POS', sep='\t')
                frequency_df[pop] = df.iloc[:, 0]

            if with_ref:
                ref = allel.read_vcf(f'{file}.vcf.gz', fields=['POS', 'REF'])
                df = pd.DataFrame({'REF': ref['variants/REF']}, index=ref['variants/POS'])
                frequency_df = frequency_df.join(df)

            frequency_df.to_csv(self._freq_path(part))

        logger.info(f'Job for {self._query(part)} done in {time.monotonic() - start_time} seconds.')

    def frequency(self, part):
        if not self._check_local(self._freq_path(part)):
            self._calculate_freq(part)
        return pd.read_csv(self._freq_path(part), index_col=0)

    def _sample_genotype_path(self, sample, part):
        return str(self.base_path / 'data' / f'part_size{self.part_size}' / 'samples' / f'{sample}.part{part}.vcf.gz')

    def sample_genotype_array(self, sample, part):
        """Get a genotype array for the specified individual"""
        file = self._sample_genotype_path(sample, part)

        if not self._check_local(file):
            cmd = f'bcftools view -s {sample} -v snps -m2 -M2 -Oz -o {file} {self._chr_path} {self._query(part)}'
            subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)

        gt = allel.read_vcf(file, fields=['GT', 'POS'])
        return pd.Series(allel.GenotypeArray(gt['calldata/GT'])[:, 0], index=gt['variants/POS'])

    def closest_populations(self, sample, part, n=10):
        """Return an ordered list of Populations from best to worst"""

        if n >= len(self.populations):
            logger.info(f'Closest populations for {sample}: {self.populations}')
            return self.populations

        freq = self.frequency(part)
        freq[sample] = self.sample_genotype_array(sample, part)

        likeliness = {pop: 0 for pop in self.populations}
        for pos, row in freq.iterrows():
            left, right = row[sample]
            for pop in self.populations:
                likeliness[pop] += row[pop] if left == 0 else 1 - row[pop]
                likeliness[pop] += row[pop] if right == 0 else 1 - row[pop]

        best = [k for k in sorted(likeliness, key=likeliness.get, reverse=True)]
        logger.info(f'Closest populations for {sample}: {best[:n]}')
        return best[:n]

    def real_population(self, sample):
        """Return a real population for the specified individual"""
        df = self.samples_df()
        return df.loc[sample, 'pop']

    def calculate_informativeness(self, part):
        start_time = time.monotonic()
        pairs = list(combinations(self.populations, 2))
        with open(str(self.base_path / 'data' / 'informativeness' / f'inf.part{part}.csv'), 'w') as f:
            f.write('POS' + '\t'.join([f'{a}_{b}' for a, b in pairs]) + '\n')
            for pos, row in self.frequency(part).iterrows():
                inf = [i4a(row[first], row[second]) for first, second in pairs]
                f.write('\t'.join([str(pos)] + inf) + '\n')
                logger.info(f'Row #{pos} done, current duration: {time.monotonic() - start_time}')

    def _markers_path(self, part):
        return str(self.base_path / 'data' / f'part_size{self.part_size}'
                   / 'informativeness' / f'markers.part{part}.csv')

    def get_best_markers(self, p1, p2):
        """Find the best markers based on probabilities of REF allele for the pair of populations"""
        q = PriorityQueue()
        for a, b, pos in zip(p1, p2, p1.index.values):
            if a != b:
                if q.qsize() < self.n:
                    q.put((i4a(a, b), pos))
                else:
                    last = q.get()
                    q.put(max((i4a(a, b), pos), last, key=itemgetter(0)))
        return reversed([q.get() for _ in range(q.qsize())])

    def _calculate_markers(self, part):
        """
        Write best markers with informativeness to a file.
        The output file has info about all pairs of populations.
        """
        freq = self.frequency(part)
        pairs = list(combinations(self.populations, 2))
        with open(self._markers_path(part), 'w') as f:
            for first, second in pairs:
                f.write(' '.join((first, second, '')))
                f.write(' '.join(f'{i} {m}' for i, m in self.get_best_markers(freq[first], freq[second])) + '\n')

                logger.info(f'{(first, second)} done.')

    def markers_df(self, part):
        if not self._check_local(self._markers_path(part)):
            self._calculate_markers(part)
        header = sum([[f'INF{i}', f'MARKER{i}'] for i in range(self.n)], [])
        return pd.read_csv(self._markers_path(part), index_col=(0, 1), names=header, sep=' ')

    def generate_features_as_json(self, part):
        t = time.monotonic()
        markers_df = self.markers_df(part)
        freq = self.frequency(part)
        result = {}
        for sample in self.samples():
            logger.info(f'Running feature calculator for {sample}')
            genotype = self.sample_genotype_array(sample, part)
            populations = sorted(self.closest_populations(sample, part))
            score = {'real': self.real_population(sample), 'closest': populations}
            for first, second in combinations(populations, 2):
                first_score, second_score = 0, 0
                line = markers_df.xs((first, second)).values
                for marker, inf in zip(line[1::2], line[::2]):
                    for allele in genotype[marker]:
                        first_freq, second_freq = freq.loc[int(marker), first], freq.loc[int(marker), second]
                        if allele == 0:
                            if first_freq > second_freq:
                                first_score += inf
                            else:
                                second_score += inf
                        else:
                            if first_freq > second_freq:
                                second_score += inf
                            else:
                                first_score += inf
                score[f'{first}_{second}'] = (first_score, second_score)
            logger.info(f'Current execution time: {time.monotonic() - t}')
            result[sample] = score

        json.encoder.FLOAT_REPR = lambda o: format(o, '.4f')
        with open(f'features_{self._query(part)}.json', 'w') as f:
            json.dump(result, f, indent=4)

    def generate_features_as_csv(self, part):
        t = time.monotonic()
        logger.info(f'Generating features. Start time {t}')
        markers_df = self.markers_df(part)
        freq = self.frequency(part)
        for sample in self.samples():
            logger.info(f'Running feature calculator for {sample}')
            genotype = self.sample_genotype_array(sample, part)
            populations = sorted(self.closest_populations(sample, part))
            real = self.real_population(sample)
            scores = []
            for first, second in combinations(populations, 2):
                first_score, second_score = 0, 0
                line = markers_df.xs((first, second)).values
                for marker, inf in zip(line[1::2], line[::2]):
                    for allele in genotype[marker]:
                        first_freq, second_freq = freq.loc[int(marker), first], freq.loc[int(marker), second]
                        if allele == 0:
                            if first_freq > second_freq:
                                first_score += inf
                            else:
                                second_score += inf
                        else:
                            if first_freq > second_freq:
                                second_score += inf
                            else:
                                first_score += inf
                scores.append(f"{first_score / (first_score + second_score):.5f}")
            logger.info(f'Current execution time: {time.monotonic() - t}')
            with open(f'features_{self._query(part)}.csv', 'a') as f:
                f.write(f'{sample} {real} {" ".join(populations)} {" ".join(scores)}\n')


class ContinentalFeatureCalc(PartialInfoCalc):

    def __init__(self, chromosome=1, n_markers=7, populations=None, part_size=1):
        super().__init__(chromosome=chromosome, n_markers=n_markers, part_size=part_size)
        self.subset = sorted(populations or ['CEU', 'CDX', 'CLM', 'GIH', 'MSL'])

    def subset_samples_df(self, pop='ALL'):
        header = 0 if pop == 'ALL' else None
        index = 'sample' if pop == 'ALL' else 0
        df = pd.read_csv(self._samples_path(pop), sep='\t', header=header, index_col=index)
        return df[df.iloc[:, 0].isin(self.subset)]

    def subset_samples(self, pop='ALL'):
        df = self.subset_samples_df(pop)
        return df.index.values

    def generate_features_as_csv(self, part):
        t = time.monotonic()
        logger.info(f'Generating features for populations: {self.populations}.\n'
                    f'Start time {t}')
        markers_df = self.markers_df(part)
        pairs = list(combinations(self.subset, 2))
        freq = self.frequency(part)
        out_file = f'features_part{part}.csv'
        with open(out_file, 'w') as f:
            f.write(f'{make_header(pairs)}\n')
            for sample in self.subset_samples():
                logger.info(f'Running feature calculator for {sample}')
                genotype = self.sample_genotype_array(sample, part)
                real = self.real_population(sample)
                scores = []
                for first, second in pairs:
                    first_score, second_score = 0, 0
                    line = markers_df.xs((first, second)).values
                    for marker, inf in zip(line[1::2], line[::2]):
                        for allele in genotype[marker]:
                            first_freq, second_freq = freq.loc[int(marker), first], freq.loc[int(marker), second]
                            if allele == 0:
                                if first_freq > second_freq:
                                    first_score += inf
                                else:
                                    second_score += inf
                            else:
                                if first_freq > second_freq:
                                    second_score += inf
                                else:
                                    first_score += inf
                    scores.append(f"{first_score / (first_score + second_score):.5f}")

                logger.info(f'Current execution time: {time.monotonic() - t}')

                f.write(f'{sample},{real},{",".join(scores)}\n')

        logger.info(f'Done! Results are avaliable at {out_file}\n')


class LargeRegionFeatureCalc(ContinentalFeatureCalc):

    def _init_dir(self):
        freq_path = self.base_path / 'data' / 'ld50' / 'frequencies'
        marker_path = self.base_path / 'data' / f'ld50' / 'informativeness'
        genotype_path = self.base_path / 'data' / f'ld50' / 'genotypes'

        for p in (freq_path, marker_path, genotype_path):
            p.mkdir(parents=True, exist_ok=True)

    def _freq_path(self, part):
        return str(self.base_path / 'data' / 'ld50' / 'frequencies' / f'ALL.frq.part{part}')

    def _markers_path(self, part):
        return str(self.base_path / 'data' / 'ld50' / 'informativeness' / f'markers.part{part}.csv')

    def _sample_genotype_path(self, sample, part):
        return str(self.base_path / 'data' / 'ld50' / 'genotypes' / f'{sample}.part{part}.vcf.gz')

    @property
    @lru_cache(1)
    def segments(self):
        segment_path = str(BASE_PATH / 'data' / 'plink' / 'genetic_map_chr1_003.txt')
        logger.info(f'Getting segments for chromosome{self.chr} from {segment_path}')
        segment_df = pd.read_csv(segment_path)
        segments = [int(x) for x in segment_df['position']]
        return segments

    def find_separable_populations(self, part, num=5):
        names = ['POP1', 'POP2'] + " ".join([f'SCORE{i} MARKER{i}' for i in range(self.n)]).split()
        df = pd.read_csv(self._markers_path(part), header=None, index_col=[0, 1], names=names, sep=' ')
        df['TOTAL_SCORE'] = sum(df[f'SCORE{i}'] for i in range(self.n))

        df = df.sort_values(by=['TOTAL_SCORE'], ascending=False)
        df.to_csv(f'pop_scores_part{part}', index=False)

        logger.info(f'Choosing best 5 populations from {fac(26) / fac(num) / fac(26 - num)} for part {part}')
        best = []
        best_score = 0
        for pops in combinations(self.populations, num):
            score = 0
            for pop1, pop2 in combinations(pops, 2):
                score += df.xs([pop1, pop2])['TOTAL_SCORE']
            if score > best_score:
                best = pops
                best_score = score
        return best

    def set_subset(self, subset):
        self.subset = subset


def main(args):
    pic = LargeRegionFeatureCalc(
        chromosome=args.chromosome,
        populations=sorted(POPULATIONS),
        n_markers=args.n
    )
    for part in [102, 207, 312, 413]:
        pops = pic.find_separable_populations(part, 3)
        print(part)
        print(pops)
        pic.set_subset(pops)
        pic.generate_features_as_csv(part)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate best markers and their informativeness for specified '
                                                 'segment of the chromosome')

    parser.add_argument('population codes', metavar='pop', type=str, nargs='*', default=('GBR', 'FIN'),
                        help='Populations for the informativeness calculation')
    parser.add_argument('-n', action='store', type=int, default=10, help='number of markers')
    parser.add_argument('-c', '--chromosome', action='store', type=int, default=1, help='Chromosome number')
    parser.add_argument('-s', '--start', action='store', type=int, default=1, help='Starting position')
    parser.add_argument('-e', '--end', action='store', type=int, default=1000000, help='Ending position')

    args = parser.parse_args()

    main(args)
