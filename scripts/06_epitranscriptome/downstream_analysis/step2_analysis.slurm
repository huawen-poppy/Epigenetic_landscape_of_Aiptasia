#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J analysis_02
#SBATCH -o ./logs/ana_step2.%J.out
#SBATCH -e ./logs/ana_step2.%J.err
#SBATCH --time=18:30:00
#SBATCH --mem=200G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --ntasks-per-node=1

conda activate epigenetic
module load gcc
#python summarize_and_improve_tandemmod_results.py
# this step is to analysis the differential modification rate between the two conditions for all the modification types.
python tandemmod_site_level_differential_analysis.py
python summarize_and_improve_tandemmod_results.py
