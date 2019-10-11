import pathlib
import vcf
import subprocess

from queue import PriorityQueue
from operator import itemgetter
from math import log

from src.populations import Population


BASE_PATH = pathlib.Path(__file__).parents[1]


class FrequencyCalculator:

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

    def _check_local(self):
        # TODO: implement this
        return False

    def fetch(self):
        if self._check_local():
            return
        file = f'{self.base_url}{self.pattern.format(self.chromosome)}'
        query = f'{self.chromosome}:{self.start}-{self.end}'
        local_filename = self._segment_path()

        load_command = ' '.join(['tabix', '-h', file, query, '|', 'bgzip', '-c', '>', local_filename])
        subprocess.call(load_command, shell=True, stdout=subprocess.PIPE)

        index_command = ' '.join(['tabix', '-p', 'vcf', local_filename])
        subprocess.call(index_command, shell=True, stdout=subprocess.PIPE)

    def index(self, filename):
        pass

    def samples(self, pop):
        command = ' '.join(['grep', pop, self._path('samples.tsv'), '|', 'cut', '-f1', '>', self._samples_path(pop)])
        subprocess.call(command, shell=True, stdout=subprocess.PIPE)

    def subset(self, pop):
        local_filename = self._segment_path()
        command = ' '.join(['vcf-subset', '-c', self._samples_path(pop), local_filename, '|', 'fill-an-ac',
                           '|', 'bgzip', '-c', '>',  self._pop_segment_path(pop)])
        subprocess.call(command, shell=True)

    def calculate_frequency(self, pop):
        reader = vcf.Reader(filename=self._pop_segment_path(pop))
        dict = {}
        for record in reader:
            allele_num = record.INFO['AN']
            allele_counts = [allele_num - sum(record.INFO['AC'])] + record.INFO['AC']
            dict[record.POS] = [ac / allele_num for ac in allele_counts]
        return dict

    def execute(self):
        self.fetch()
        population_frequencies = {}
        for pop in self.populations:
            self.samples(pop)
            self.subset(pop)
            population_frequencies[pop] = self.calculate_frequency(pop)
        print(population_frequencies)


def i4a(record, populations=None):
    """Calculate informativeness for assignment"""
    alleles = [record.REF] + record.ALT
    populations = populations or Population.superpopulations
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


