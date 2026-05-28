#!/bin/bash

targets=('A1' 'A2' 'A3' 'A4' 'H1' 'H2' 'H3' 'H4')
pre='../raw_data/'
input_path='/bamtx/aligned_new_isoquant_transcriptome.sam'
output_path='/bamtx/aligned_new_isoquant_transcriptome.bam'
output_path2='/bamtx/aligned_new_isoquant_transcriptome.bam.log'

for target in ${targets[@]};do
    input_a=$pre$target$input_path
    input_b=$pre$target$output_path
    input_c=$pre$target$output_path2
    echo $input_a
    echo $input_b
    echo $input_c
    sbatch step3_samtools.slurm $input_a $input_b $input_c $input_b $input_c
done
