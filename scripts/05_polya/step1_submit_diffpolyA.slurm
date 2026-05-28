#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J polya_diff
#SBATCH -o ./pipeline-polya-diff/logs/polya.%J.out
#SBATCH -e ./pipeline-polya-diff/logs/polya.%J.err
#SBATCH --time=168:30:00
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1

#run the application:
conda activate epigenetic
module load snakemake/7.32.3

snakemake -j 36 all
