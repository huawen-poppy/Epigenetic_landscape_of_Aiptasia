#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J step5
#SBATCH -o ./logs/step5.%J.out
#SBATCH -e ./logs/step5.%J.err
#SBATCH --time=12:30:00
#SBATCH --mem=700G
#SBATCH --gres=gpu:1
#SBATCH --constraint=v100

conda activate epigenetic 
module load gcc
python scripts/TandemMod.py --run_mode predict \
    --pretrained_model $1 \
    --feature_file $2 \
    --predict_result $3

