#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J analysis_06
#SBATCH -o ./logs/ana_step6.%J.out
#SBATCH -e ./logs/ana_step6.%J.err
#SBATCH --time=18:30:00
#SBATCH --mem=200G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --ntasks-per-node=1

conda activate scanpy
module load gcc   
# this step is to analyze the differential modification transcripts position
python analyze_differential_modification_position.py

