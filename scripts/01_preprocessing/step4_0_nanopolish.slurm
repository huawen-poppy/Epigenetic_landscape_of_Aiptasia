#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J step3_samtools
#SBATCH -o ./logs/step3.%J.out
#SBATCH -e ./logs/step3.%J.err
#SBATCH --time=03:30:00
#SBATCH --mem=64G
#
##run the application:
conda activate epigenetic
# nanopolish index -d <PATH/TO/FAST5_DIR> <PATH/TO/FASTQ_FILE>
# for this step, the fastq file need to be gunzip 
# we need to provide the sequencing summary files

nanopolish index -d $1 -s $2 $3
 
