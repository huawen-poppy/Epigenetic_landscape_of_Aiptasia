#!/bin/bash

targets=('A1' 'A2' 'A3' 'A4' 'H1' 'H2' 'H3' 'H4')
pre='../raw_data/'
input_path='/fast5'
output_path='/fastq'

for target in ${targets[@]};do
    input_a=$pre$target$input_path
    input_b=$pre$target$output_path
    echo $input_a
    echo $input_b
    sbatch step0_basecaller_raw_bash.slurm $input_a $input_b
done
