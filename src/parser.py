import pathlib
import vcf

from src.analyser import get_best_markers
from src.populations import Population

path = pathlib.Path(__file__).parents[1] / 'data' / 'ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
reader = vcf.Reader(filename=str(path))


pop = Population()
print(pop.pair)

markers = get_best_markers(reader.fetch('1', 10000, 100000), pop.pair)
for marker in markers:
    inf, rec = marker
    print(inf, rec.POS)

