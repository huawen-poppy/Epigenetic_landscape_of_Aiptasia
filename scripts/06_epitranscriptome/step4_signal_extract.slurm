#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J step4
#SBATCH -o ./logs/step4.%J.out
#SBATCH -e ./logs/step4.%J.err
#SBATCH --time=24:30:00
#SBATCH --mem=300G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=1

conda activate epigenetic 
#python scripts/extract_signal_from_fast5.py -p 36 --fast5 ../raw_data_pass/A1/guppy_single/ --reference /ibex/project/c2101/aip_epitrans_analysis/rnabloom/busco_test/isoquant_input/isoquant_extend_transcript.fa --sam ../raw_data_pass/A1/A1.sam --output ../raw_data_pass/A1/A1.signal.tsv --clip 10
#python scripts/extract_feature_from_signal.py  --signal_file ../raw_data_pass/A1/A1.signal.tsv --clip 10 --output ../raw_data_pass/A1/A1.feature.tsv 

python scripts/extract_signal_from_fast5.py -p 36 --fast5 $1 --reference /ibex/project/c2101/aip_epitrans_analysis/rnabloom/busco_test/isoquant_input/isoquant_extend_transcript.fa --sam $2 --output $3 --clip 10
python scripts/extract_feature_from_signal.py  --signal_file $3 --clip 10 --motif DRACH --output $4
python scripts/extract_feature_from_signal.py  --signal_file $3 --clip 10 --motif NNANN --output $5
python scripts/extract_feature_from_signal.py  --signal_file $3 --clip 10 --motif NNCNN --output $6
python scripts/extract_feature_from_signal.py  --signal_file $3 --clip 10 --motif NNGNN --output $7
python scripts/extract_feature_from_signal.py  --signal_file $3 --clip 10 --motif NNTNN --output $8


# the NN*NN process takes very long time, so do the below:
# (1) chunk the signal files into pieces 
# (2) submit the nn*nn job for each piece
# (3) cat the results together

