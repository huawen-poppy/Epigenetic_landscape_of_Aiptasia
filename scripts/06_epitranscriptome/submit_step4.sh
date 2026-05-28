#!/bin/bash

targets=('A1' 'A2' 'A3' 'H2' 'H1''A1')
pre='/ibex/scratch/projects/c2101/aip_epitrans_analysis/raw_data_pass/'
input_path='/guppy'
output_path='/guppy_single'

motifs=("DRACH" "NNANN" "")
for target in ${targets[@]};do
    input_a=$pre$target"/guppy_single"
    input_b=$pre$target"/"$target".sam"
    input_c=$pre$target"/"$target".signal.tsv"
    input_d=$pre$target"/"$target"_DRACH.feature.tsv"
    input_e=$pre$target"/"$target"_NNANN.feature.tsv"
    input_f=$pre$target"/"$target"_NNCNN.feature.tsv"
    input_g=$pre$target"/"$target"_NNGNN.feature.tsv"
    input_h=$pre$target"/"$target"_NNTNN.feature.tsv"
    echo $input_a
    echo $input_b
    sbatch step4_signal_extract.slurm $input_a $input_b $input_c $input_d $input_e $input_f $input_g $input_h
done 

#python scripts/extract_signal_from_fast5.py -p 36 --fast5 ../raw_data_pass/A1/guppy_single/ --reference /ibex/project/c2101/aip_epitrans_analysis/rnabloom/busco_test/isoquant_input/isoquant_extend_transcript.fa --sam ../raw_data_pass/A1/A1.sam --output ../raw_data_pass/A1/A1.signal.tsv --clip 10
#python scripts/extract_feature_from_signal.py  --signal_file ../raw_data_pass/A1/A1.signal.tsv --clip 10 --output ../raw_data_pass/A1/A1.feature.tsv
