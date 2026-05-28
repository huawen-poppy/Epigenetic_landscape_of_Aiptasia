#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J step2_minimap2
#SBATCH -o ./logs/step2_minimap2.%J.out
#SBATCH -e ./logs/step2_minimap2.%J.err
#SBATCH --time=05:30:00
#SBATCH --mem=64G
#SBATCH --gres=gpu:4
#SBATCH --constraint=[gpu]
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-gpu=6


#run the application:
conda activate epigenetic

#minimap2 -ax map-ont splice -t 8 -uf -k14 [reference].fasta.mmi [sample].fastq > [sample].sam
#need to cat *.fastq > combind_all.fastq
#then gzip the combind_all.fastq file

minimap2 -ax map-ont -uf -t 8 --secondary=no ../raw_data/aiptasia_mrna_new_isoquant_transcript.fasta.mmi $1 > $2 2>> $3
