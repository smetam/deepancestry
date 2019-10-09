from src.populations import Population
from math import log
from queue import PriorityQueue
from operator import itemgetter


def i4a(record, populations=None):
    """Calculate informativeness for assignment"""
    alleles = [record.REF] + record.ALT
    populations = populations or Population.all
    s = 0
    N = len(alleles)
    K = len(populations)
    # population_allele_frequencies
    paf = {
        pop: [1 - sum(record.INFO[pop + '_AF'])] + record.INFO[pop + '_AF'] for pop in populations
    }
    mean_frequencies = [sum(paf[pop][i] for pop in populations) / K for i in range(N)]
    for j in range(N):
        p = mean_frequencies[j]
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

