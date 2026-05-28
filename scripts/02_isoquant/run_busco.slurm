#!/bin/bash
#SBATCH --job-name=prokka
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=200G

# BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs
# https://academic.oup.com/bioinformatics/article/31/19/3210/211866

module load busco

#set -vex
conda activate epigenetic

inFile='./isoquant_input/isoquant_extend_transcript.fa'
outPref='./out_isoquant_extend_data'

#run_BUSCO.py -i $inFile -o $outPref -l $lib -c 5 -sp E_coli_K12 -m prot > ecoli_80x_busco.log 2>&1
busco -i $inFile -o $outPref -m transcriptome -l metazoa_odb10 -c 36
