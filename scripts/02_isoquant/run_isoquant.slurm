#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J isoquant
#SBATCH -o ./logs/ref_all.%J.out
#SBATCH -e ./logs/ref_all.%J.err
#SBATCH --time=168:30:00
#SBATCH --mem=100G


#run the application:
# based on the https://github.com/ablab/IsoQuant
conda activate epigenetic

isoquant.py --reference ../raw_data/aiptasia_genome_final.fa \
--genedb ../raw_data/aiptasia_genome.dups_removed.gtf \
--fastq  ../raw_data/A1/fastq/combind_all.fastq ../raw_data/A2/fastq/combind_all.fastq ../raw_data/A3/fastq/combind_all.fastq ../raw_data/A4/fastq/combind_all.fastq ../raw_data/H1/fastq/combind_all.fastq ../raw_data/H2/fastq/combind_all.fastq ../raw_data/H3/fastq/combind_all.fastq ../raw_data/H4/fastq/combind_all.fastq \
--labels A1 A2 A3 A4 H1 H2 H3 H4 \
--threads 36 \
--data_type nanopore -o ./isoquant/outputs/final_all_ref --sqanti_output --check_canonical --count_exons 
