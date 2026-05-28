#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J step1
#SBATCH -o ./logs/step1.%J.out
#SBATCH -e ./logs/step1.%J.err
#SBATCH --time=48:30:00
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks-per-node=1

conda activate epigenetic
#multi_to_single_fast5 -i /ibex/scratch/projects/c2101/aip_epitrans_analysis/raw_data_pass/A1/guppy -s /ibex/scratch/projects/c2101/aip_epitrans_analysis/raw_data_pass/A1/guppy_single -t 40 --recursive

multi_to_single_fast5 -i $1 -s $2 -t 40 --recursive
