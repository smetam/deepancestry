import pathlib
import vcf
import sys
import logging
import subprocess
import pandas as pd

from queue import PriorityQueue
from operator import itemgetter
from math import log

from src.populations import Population


BASE_PATH = pathlib.Path(__file__).parents[1]

logger = logging.getLogger(__name__)
logger.setLevel('INFO')
logger.addHandler(logging.StreamHandler(sys.stdout))


class InformativenessCalculator:

    def __init__(self, chromosome, start, end, populations=None):
        self.populations = populations or Population().pair
        self._configure()
        self._set_segment(chromosome, start, end)

    def _set_segment(self, chromosome, start, end):
        """Set chromosome and segment"""
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
        return self._path(f'{pop}.samples.list')

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

    def index(self, filename):
        pass

    def samples(self, pop):
        logger.info(f'Retrieving samples for {pop} to {self._samples_path(pop)}')
        command = ' '.join(['grep', pop, self._path('samples.tsv'), '|', 'cut', '-f1', '>', self._samples_path(pop)])
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

    def calculate_frequency(self, pop):
        logger.info(f'Calculating frequencies for {pop}')
        reader = vcf.Reader(filename=self._pop_segment_path(pop))
        ref, alt = f'{pop}_REF', f'{pop}_ALT'
        frequency_df = pd.DataFrame(columns=[ref, alt])
        for record in reader:
            if len(record.ALT) == 1 and record.is_snp:
                allele_num = record.INFO['AN']
                frequency_df.loc[record.POS] = (allele_num - record.INFO['AC'][0]) / allele_num, record.INFO['AC'][0] / allele_num
        return frequency_df

    def alleles(self):
        reader = vcf.Reader(filename=self._segment_path())
        alleles = pd.DataFrame(columns=['REF', 'ALT'])
        logger.info('Filtering biallelic and snp records')
        for record in reader:
            if len(record.ALT) == 1 and record.is_snp:
                alleles.loc[record.POS] = record.REF, record.ALT[0]
        return alleles

    def calculate_frequencies(self):
        if not self._check_local(self._segment_path()):
            logger.info(f'Could not find local segment at path: {self._segment_path()}')
            self.fetch()
        population_frequencies = self.alleles()
        for pop in self.populations:
            if not self._check_local(self._pop_segment_path(pop)):
                logger.info(f'Could not find local segment for {pop} at path: {self._pop_segment_path(pop)}')
                self.samples(pop)
                self.subset(pop)
            population_frequencies = population_frequencies.join(self.calculate_frequency(pop))
        return population_frequencies


def i4a(record, frequency_df, populations=None):
    """Calculate informativeness for assignment"""
    alleles = [record.REF] + record.ALT
    populations = populations or list(frequency_df)
    s = 0
    N = len(alleles)
    K = len(populations)
    # population allele frequencies
    paf = {
        pop: [1 - sum(record.INFO[pop + '_AF'])] + record.INFO[pop + '_AF'] for pop in populations
    }
    # average allele frequencies
    aaf = [sum(paf[pop][i] for pop in populations) / K for i in range(N)]
    for j in range(N):
        p = aaf[j]
        s -= p * log(p) if p != 0 else 0
        for pop in populations:
            q = paf[pop][j]
            s += q / K * log(q) if q != 0 else 0
    return s


def get_best_markers(records, populations, capacity=10):
    """Get a list of most informative markers in tuples (informativeness, record)"""
    q = PriorityQueue()
    for record in records:
        informativeness = i4a(record, populations)
        if q.qsize() < capacity:
            q.put((informativeness, record))
        else:
            last = q.get()
            q.put(max((informativeness, record), last, key=itemgetter(0)))
    return [q.get() for _ in range(q.qsize())]


if __name__ == '__main__':
    fc = InformativenessCalculator(1, 1, 100000, populations=('GBR', 'FIN'))
    df = fc.calculate_frequencies()
    print(df)
