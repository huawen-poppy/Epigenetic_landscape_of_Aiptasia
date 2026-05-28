#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J nanopolish_polyA
#SBATCH -o ./logs/nanopolish.%J.out
#SBATCH -e ./logs/nanopolish.%J.err
#SBATCH --time=48:30:00
#SBATCH --mem=100G


#run the application:
conda activate epigenetic

#Usage: nanopolish polya [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa
# after that, we need to filter the output file to select only the 'PASS' tag in the last column qc_tag
# grep 'PASS' file > polya_filter_pass_only.tsv


nanopolish polya --reads $1 \
--bam $2 \
--genome /ibex/project/c2101/aip_epitrans_analysis/rnabloom/busco_test/isoquant_input/isoquant_extend_transcript.fa \
--threads 36 > $3
