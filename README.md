# deepancestry
Deep learning approach to determination of local ancestry

## Purpose of the project
Evaluate efficiency of neural networks and 

## Main ideas:
I. Split genome into regions according to recombination hot spots map.

II. For each region of the genome: 

1. Take biallelic SNPs
2. Count frequencies for alleles in each population
3. For each pair of populations find n = 15 most informative marker using infocalc algorithm
4. Take a persons genotype array and find 10 closest populations according to Hamming distance
5. For each pair of populations find out which one is more likely out of the two
6. Sort these populations by looking at the comparison of pairs

III. Compare efficiency of different methods

## Dependencies:
* Python 3.7
* vcftools
* samtools/tabix/bcftools
 
 Python dependencies can be installed with pip:
 
 `
 pip install -r requirements.txt
 `
 
## Usage

Calculate markers informativeness by specifying chromosome and populations

`
python3 src/preprocessing -c 1 FIN ESN TSI
`
Example of results

| SAMPLE  | REAL | CLM_ESN | CLM_KHV | CLM_PJL | CLM_TSI | ESN_KHV | ESN_PJL | ESN_TSI | KHV_PJL | KHV_TSI | PJL_TSI |
|---------|------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
| HG01112 | CLM  | 0.10682 | 0.54587 | 0.60590 | 0.74482 | 0.92378 | 0.85410 | 0.89658 | 0.27562 | 0.63720 | 1.00000 |
| HG01113 | CLM  | 0.65587 | 0.93624 | 0.72565 | 0.90181 | 0.92378 | 0.48376 | 0.65378 | 0.05737 | 0.38286 | 1.00000 |
| HG01119 | CLM  | 0.40114 | 0.00000 | 0.49395 | 0.15712 | 0.00000 | 0.36776 | 0.36066 | 0.94487 | 0.88404 | 0.41088 |
| HG01121 | CLM  | 0.43143 | 0.77294 | 0.66578 | 0.82665 | 0.92378 | 0.61627 | 0.73108 | 0.13781 | 0.51003 | 1.00000 |
| HG01122 | CLM  | 0.10682 | 0.54587 | 0.60590 | 0.74482 | 0.92378 | 0.85410 | 0.89658 | 0.27562 | 0.63720 | 1.00000 |

## Links 
1. 1000 genomes project
https://www.internationalgenome.org

2. Infocalc 
https://rosenberglab.stanford.edu/infocalc.html
