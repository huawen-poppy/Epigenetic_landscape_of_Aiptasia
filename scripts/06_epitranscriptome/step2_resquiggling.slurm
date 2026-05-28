#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J step0_prepare
#SBATCH -o ./logs/step2.%J.out
#SBATCH -e ./logs/step2.%J.err
#SBATCH --time=24:30:00
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks-per-node=1

conda activate tombo_env 

tombo resquiggle --rna --overwrite --basecall-group Basecall_1D_000 --basecall-subgroups BaseCalled_template --processes 40 $1 /ibex/project/c2101/aip_epitrans_analysis/rnabloom/busco_test/isoquant_input/isoquant_extend_transcript.fa --ignore-read-locks
