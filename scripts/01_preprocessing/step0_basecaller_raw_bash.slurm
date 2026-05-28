#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J step0_prepare
#SBATCH -o ./logs/step0_prepare.%J.out
#SBATCH -e ./logs/step0_prepare.%J.err
#SBATCH --time=24:30:00
#SBATCH --mem=64G
#SBATCH --gres=gpu:4
#SBATCH --constraint=[gpu]
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-gpu=6


#run the application:
conda activate epigenetic
module load guppy/3.4.1
# --flowcell FLO-MIN106 --kit SQK-RNA002
guppy_basecaller -i $1 -s $2 -c rna_r9.4.1_70bps_hac.cfg --device auto --cpu_threads_per_caller 24 --gpu_runners_per_device 4 --num_callers 4 -q 0 -r
