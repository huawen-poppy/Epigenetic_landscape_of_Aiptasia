#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J step0_prepare
#SBATCH -o ./logs/step0_prepare.%J.out
#SBATCH -e ./logs/step0_prepare.%J.err
#SBATCH --time=24:30:00
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1

conda activate TandemMod
module load guppy
guppy_basecaller -i $1 -s $2 --recursive --fast5_out -c rna_r9.4.1_70bps_hac.cfg -x cuda:all --num_callers 4 --gpu_runners_per_device 4 --cpu_threads_per_caller 4
