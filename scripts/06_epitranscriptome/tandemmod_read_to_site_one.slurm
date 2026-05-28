#!/bin/bash
#SBATCH --job-name=TM_site
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

INPUT_FILE="$1"
OUTPUT_FILE="$2"

echo "======================================"
echo "Input read-level prediction: ${INPUT_FILE}"
echo "Output site-level prediction: ${OUTPUT_FILE}"
echo "Started at: $(date)"
echo "======================================"

source ~/.bashrc
conda activate epigenetic
module load gcc

TANDEMMOD_DIR="/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod"

if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: input file not found: ${INPUT_FILE}"
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT_FILE")"

cd "$TANDEMMOD_DIR"

python scripts/read_level_prediction_to_site_level_prediction_fixed.py \
    --input "$INPUT_FILE" \
    --output "$OUTPUT_FILE"

echo "Finished at: $(date)"
ls -lh "$OUTPUT_FILE"
