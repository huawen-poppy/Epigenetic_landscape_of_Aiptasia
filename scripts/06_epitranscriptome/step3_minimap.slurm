#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J step3
#SBATCH -o ./logs/step3.%J.out
#SBATCH -e ./logs/step3.%J.err
#SBATCH --time=24:30:00
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --ntasks-per-node=1

conda activate epigenetic 
#cat demo/IVET/IVET_m6A_guppy/pass/*.fastq >demo/IVET/IVET_m6A.fastq
#minimap2 -ax map-ont demo/IVET_reference.fa demo/IVET/IVET_m6A.fastq >demo/IVET/IVET_m6A.sam

#minimap2 -d ../raw_data_pass/isoquant_extend_transcript.mmi \
#    /ibex/project/c2101/aip_epitrans_analysis/rnabloom/busco_test/isoquant_input/isoquant_extend_transcript.fa -t 12

#minimap2 -ax map-ont ../raw_data_pass/isoquant_extend_transcript.mmi /ibex/scratch/projects/c2101/aip_epitrans_analysis/raw_data_pass/A1/A1.fastq > /ibex/scratch/projects/c2101/aip_epitrans_analysis/raw_data_pass/A1/A1.sam -t 12

minimap2 -ax map-ont ../raw_data_pass/isoquant_extend_transcript.mmi $1 > $2 -t 12
