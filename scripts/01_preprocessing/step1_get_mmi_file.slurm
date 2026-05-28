#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J step1_mmi
#SBATCH -o ./logs/step1_mmi.%J.out
#SBATCH -e ./logs/step1_mmi.%J.err
#SBATCH --time=10:30:00
#SBATCH --mem=64G
#SBATCH --gres=gpu:4
#SBATCH --constraint=[gpu]
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-gpu=6


#run the application:
conda activate epigenetic

#the isoquant_extend_transcript.fa is the isoquant assembly transcript fasta file
minimap2 -ax map-ont -t 8 -uf -k14 -d ../raw_data/aiptasia_mrna_new_isoquant_transcript.fasta.mmi ../raw_data/isoquant_extend_transcript.fa  
