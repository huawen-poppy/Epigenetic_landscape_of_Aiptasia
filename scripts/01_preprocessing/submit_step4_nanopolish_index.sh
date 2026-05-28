#!/bin/bash

targets=('A1' 'A2' 'A3' 'A4' 'H1' 'H2' 'H3' 'H4')
input_path='/fast5'
pre2='../raw_data/'
output_path='/fastq/sequencing_summary.txt'
output_path2='/fastq/combind_all.fastq'

for target in ${targets[@]};do
    input_a=$pre$target$input_path
    input_b=$pre2$target$output_path
    input_c=$pre2$target$output_path2
    echo $input_a
    echo $input_b
    echo $input_c
    sbatch step4_0_nanopolish.slurm $input_a $input_b $input_c 
done
