#!/bin/bash

targets=('A2' 'A3' 'H2' 'H1' 'H3' 'A1')
pre='/ibex/scratch/projects/c2101/aip_epitrans_analysis/raw_data_pass/'
input_path='/guppy'
output_path='/guppy_single'

for target in ${targets[@]};do
    input_a=$pre$target"/"$target".fastq"
    input_b=$pre$target"/"$target".sam"
    echo $input_a
    echo $input_b
    sbatch step3_minimap.slurm $input_a $input_b
done 
