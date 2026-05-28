#!/bin/bash

targets=('A1' 'A2' 'A3' 'H1' 'H2' 'H3')
pre='/ibex/scratch/projects/c2101/aip_epitrans_analysis/raw_data/'
input_path='/fastq/combind_all.fastq'
output_path='/bamtx/aligned_new_isoquant_transcriptome.bam'
output_path2='/nanopolish/polya_results_new_isoquant_transcriptome.tsv'

for target in ${targets[@]};do
    input_a=$pre$target$input_path
    input_b=$pre$target$output_path
    input_c=$pre$target$output_path2
    echo $input_a
    echo $input_b
    echo $input_c
    sbatch step0_nanopolish_polya.slurm $input_a $input_b $input_c 
done
