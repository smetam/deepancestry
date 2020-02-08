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
5. Sort these populations by looking at the informative markers from step 3

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