#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J analysis_01
#SBATCH -o ./logs/ana_step1.%J.out
#SBATCH -e ./logs/ana_step1.%J.err
#SBATCH --time=48:30:00
#SBATCH --mem=200G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --ntasks-per-node=1

conda activate epigenetic
python summarize_tandemmod_predictions.py
