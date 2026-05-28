#!/bin/bash

targets=('A1' 'A3' 'A2' 'H1' 'H3' 'H2')
pre='/ibex/scratch/projects/c2101/aip_epitrans_analysis/raw_data_pass/'
input_path='/guppy'
output_path='/guppy_single'

for target in ${targets[@]};do
    input_a=$pre$target$input_path
    input_b=$pre$target$output_path
    echo $input_a
    echo $input_b
    sbatch step1_multi_fast5_to_single.slurm  $input_a $input_b
done 
