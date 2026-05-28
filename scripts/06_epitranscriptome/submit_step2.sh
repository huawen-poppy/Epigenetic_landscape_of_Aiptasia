#!/bin/bash

targets=('A2' 'A3' 'H2' 'H1' 'H3' 'A1')
pre='/ibex/scratch/projects/c2101/aip_epitrans_analysis/raw_data_pass/'
input_path='/guppy'
output_path='/guppy_single'

for target in ${targets[@]};do
    input_a=$pre$target$input_path
    input_b=$pre$target$output_path
    echo $input_b
    sbatch step2_resquiggling.slurm $input_b
done 
