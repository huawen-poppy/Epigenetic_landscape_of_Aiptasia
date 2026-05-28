#!/bin/bash

targets=('A1' 'A2' 'A3' 'A4' 'H1' 'H2' 'H3' 'H4')
pre='../raw_data/'
input_path='/fastq/combind_all.fastq'
output_path='/bamtx/aligned_new_isoquant_transcriptome.bam'
output_path2='/nanopolish/summary_new_isoquant_transcriptome.txt'
output_path3='/nanopolish/eventalign_new_isoquant_transcriptome.txt'

for target in ${targets[@]};do
    input_a=$pre$target$input_path
    input_b=$pre$target$output_path
    input_c=$pre$target$output_path2
    input_d=$pre$target$output_path3
    echo $input_a
    echo $input_b
    echo $input_c
    echo $input_d
    sbatch step4_1_nanopolish.slurm $input_a $input_b $input_c $input_d
done
