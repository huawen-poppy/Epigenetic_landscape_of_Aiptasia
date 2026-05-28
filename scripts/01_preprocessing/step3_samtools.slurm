#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J step3_samtools
#SBATCH -o ./logs/step3.%J.out
#SBATCH -e ./logs/step3.%J.err
#SBATCH --time=03:30:00
#SBATCH --mem=64G

#run the application:
conda activate epigenetic
module load samtools 


# You get a bam (machine readable sam) file after mapping, and it contains information about mapped and unmapped reads.
# https://www.biostars.org/p/56246/
# this step is to index the mapping file
samtools view -Sb $1 | samtools sort -o $2 &>> $3

samtools index $4 &>> $5
