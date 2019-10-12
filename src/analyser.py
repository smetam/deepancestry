import vcf
import sys
import pathlib
import logging
import argparse
import subprocess
import pandas as pd

from queue import PriorityQueue
from operator import itemgetter
from math import log


BASE_PATH = pathlib.Path(__file__).absolute().parents[1]

logger = logging.getLogger(__name__)
logger.setLevel('INFO')
logger.addHandler(logging.StreamHandler(sys.stdout))


def x_logx(x):
    return 0 if x == 0 else x * log(x)


class InformativenessCalculator:

    def __init__(self, chromosome, start, end, populations=None):
        self.populations = populations
        self._configure()
        self._set_segment(chromosome, start, end)

    def _set_segment(self, chromosome, start, end):
        """Set chromosome and segment"""
        self._alleles = None
        self._frequencies = None
        self._informativeness = None
        self.chromosome = chromosome
        self.start = start
        self.end = end

    def _configure(self, local_dir=None, file_pattern=None, base_url=None):
        """Configure paths and patterns"""
        self.base_url = base_url or 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
        self.pattern = file_pattern or 'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
        self.local_dir = local_dir or BASE_PATH / 'data'

    def _path(self, file, absolute=True):
        return str(self.local_dir / file) if absolute else file

    def _samples_path(self, pop):
        """Path to samples from population"""
        return str(BASE_PATH / 'samples' / f'{pop}.samples.list')

    def _segment_path(self):
        """Path to current """
        return self._path(f'chr{self.chromosome}.{self.start}_{self.end}.vcf.gz')

    def _pop_segment_path(self, pop):
        return self._path(f'{pop}.chr{self.chromosome}.{self.start}-{self.end}.vcf.gz')

    @staticmethod
    def _check_local(path):
        return pathlib.Path(path).exists()

    def fetch(self):
        file = f'{self.base_url}{self.pattern.format(self.chromosome)}'
        query = f'{self.chromosome}:{self.start}-{self.end}'
        local_filename = self._segment_path()

        logger.info(f'Fetching from {file}')
        load_command = ' '.join(['tabix', '-h', file, query, '|', 'bgzip', '-c', '>', local_filename])
        subprocess.call(load_command, shell=True, stdout=subprocess.PIPE)

        logger.info(f'Indexing {local_filename}')
        index_command = ' '.join(['tabix', '-p', 'vcf', local_filename])
        subprocess.call(index_command, shell=True, stdout=subprocess.PIPE)

    def samples(self, pop):
        logger.info(f'Retrieving samples for {pop} to {self._samples_path(pop)}')
        command = ' '.join(['grep', pop, self._samples_path('ALL'), '|', 'cut', '-f1', '>', self._samples_path(pop)])
        subprocess.call(command, shell=True, stdout=subprocess.PIPE)

    def subset(self, pop):
        local_filename = self._segment_path()

        logger.info(f'Preparing subset for population {pop} at {self._pop_segment_path(pop)}')
        command = ' '.join(['vcf-subset', '-c', self._samples_path(pop), local_filename, '|', 'fill-an-ac',
                           '|', 'bgzip', '-c', '>',  self._pop_segment_path(pop)])
        subprocess.call(command, shell=True)

    def calculate_multiallelic_frequency(self, pop):
        reader = vcf.Reader(filename=self._pop_segment_path(pop))
        af = pd.Series()
        for record in reader:
            allele_num = record.INFO['AN']
            allele_counts = [allele_num - sum(record.INFO['AC'])] + record.INFO['AC']
            af.loc[record.POS] = [ac / allele_num for ac in allele_counts]
        return af

    def _calculate_alleles(self):
        reader = vcf.Reader(filename=self._segment_path())
        alleles = pd.DataFrame(columns=['REF', 'ALT'])
        logger.info('Filtering biallelic and snp records')
        for record in reader:
            if len(record.ALT) == 1 and record.is_snp:
                alleles.loc[record.POS] = record.REF, record.ALT[0]
        self._alleles = alleles

    @property
    def alleles(self):
        if self._alleles is None:
            self._calculate_alleles()
        return self._alleles

    def _calculate_frequency(self, pop):
        logger.info(f'Calculating frequencies for {pop}')
        reader = vcf.Reader(filename=self._pop_segment_path(pop))
        ref, alt = f'{pop}_REF', f'{pop}_ALT'
        frequency_df = pd.DataFrame(columns=[ref, alt])
        for record in reader:
            if len(record.ALT) == 1 and record.is_snp:
                allele_num = record.INFO['AN']
                frequency_df.loc[record.POS] = \
                    (allele_num - record.INFO['AC'][0]) / allele_num, record.INFO['AC'][0] / allele_num
        return frequency_df

    def _calculate_frequencies(self):
        if not self._check_local(self._segment_path()):
            logger.info(f'Could not find local segment at path: {self._segment_path()}')
            self.fetch()
        population_frequencies = self.alleles
        for pop in self.populations:
            if not self._check_local(self._pop_segment_path(pop)):
                logger.info(f'Could not find local segment for {pop} at path: {self._pop_segment_path(pop)}')
                self.samples(pop)
                self.subset(pop)
            population_frequencies = population_frequencies.join(self._calculate_frequency(pop))
        self._frequencies = population_frequencies

    def _i4a(self, af):
        """Calculate informativeness for assignment"""
        i4a = 0
        for allele in ('REF', 'ALT'):
            for pop in self.populations:
                freq = af[f'{pop}_{allele}']
                i4a += x_logx(freq)
        i4a /= len(self.populations)
        average_ref_frequency = sum([af[f'{pop}_REF'] for pop in self.populations]) / len(self.populations)
        average_alt_frequency = sum([af[f'{pop}_ALT'] for pop in self.populations]) / len(self.populations)
        i4a -= x_logx(average_alt_frequency) + x_logx(average_ref_frequency)
        return i4a

    def _calculate_informativeness(self):
        logger.info('Calculating informativeness')
        inf = pd.Series()
        for pos, row in self.frequencies.iterrows():
            inf.loc[pos] = self._i4a(row)
        self._informativeness = inf

    @property
    def frequencies(self):
        if not self._frequencies:
            self._calculate_frequencies()
        return self._frequencies

    @property
    def informativeness(self):
        if self._informativeness is None:
            self._calculate_informativeness()
        return self._informativeness

    def best_markers(self, n=10):
        q = PriorityQueue()
        for pos, informativeness in self.informativeness.iteritems():
            if q.qsize() < n:
                q.put((informativeness, pos))
            else:
                last = q.get()
                q.put(max((informativeness, pos), last, key=itemgetter(0)))
        return [q.get() for _ in range(q.qsize())]


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate best markers and their informativeness for specified '
                                                 'segment of the chromosome')

    parser.add_argument('population codes', metavar='pop', type=str, nargs='*', default=('GBR', 'FIN'),
                        help='Populations for the informativeness calculation')
    parser.add_argument('-n', action='store', type=int, default=10, help='number of markers')
    parser.add_argument('-c', '--chromosome', action='store', type=int, default=1, help='Chromosome number')
    parser.add_argument('-s', '--start', action='store', type=int, default=1, help='Starting position')
    parser.add_argument('-e', '--end', action='store', type=int, default=100000, help='Ending position')

    args = parser.parse_args()

    fc = InformativenessCalculator(args.chromosome, args.start, args.end, populations=('GBR', 'FIN'))
    markers = fc.best_markers(args.n)

    print('Best markers for the segment:')
    for marker in markers:
        print(f'marker at POS {marker[1]} with informativeness {marker[0]}')
