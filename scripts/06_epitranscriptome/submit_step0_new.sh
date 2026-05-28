#!/bin/bash

targets=('A2' 'A3' 'H1' 'H2' 'H3')
pre='/ibex/scratch/projects/c2101/aip_epitrans_analysis/raw_data_pass/'
input_path='/fast5'
output_path='/guppy'

for target in ${targets[@]};do
    input_a=$pre$target$input_path
    input_b=$pre$target$output_path
    echo $input_a
    echo $input_b
    sbatch step0_guppy.slurm  $input_a $input_b
done 
