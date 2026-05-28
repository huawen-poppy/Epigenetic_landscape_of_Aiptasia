#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J step4_nanopolish
#SBATCH -o ./logs/step4_nanopolish.%J.out
#SBATCH -e ./logs/step4_nanopolish.%J.err
#SBATCH --time=68:30:00
#SBATCH --mem=100G
#SBATCH --gpus=1

#run the application:
conda activate epigenetic


nanopolish eventalign --reads $1 \
--bam $2 \
--genome ../raw_data/isoquant_extend_transcript.fa \
--signal-index \
--scale-events \
--progress \
--verbose \
--summary $3 \
--threads 36 > $4
