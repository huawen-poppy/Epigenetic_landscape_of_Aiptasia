#!/bin/bash

mkdir -p logs

TANDEMMOD_DIR="/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod"

PRED_DIR="${TANDEMMOD_DIR}/TandemMod_prediction"
SITE_OUT_DIR="${TANDEMMOD_DIR}/TandemMod_site_level_prediction"

mkdir -p "$SITE_OUT_DIR"

SAMPLES=("A1" "A2" "A3" "A4" "H1" "H2" "H3" "H4")

MODS=("m6A" "m1A" "A_I" "m5C" "hm5C" "m7G" "G_I" "psU")

FEATURE_TYPES=("DRACH" "NNANN" "NNANN" "NNCNN" "NNCNN" "NNGNN" "NNGNN" "NNTNN")

N_MODS=${#MODS[@]}

for SAMPLE in "${SAMPLES[@]}"; do
    for ((i=0; i<N_MODS; i++)); do

        MOD="${MODS[$i]}"
        FEATURE_TYPE="${FEATURE_TYPES[$i]}"

        INPUT_FILE="${PRED_DIR}/${SAMPLE}_${MOD}_${FEATURE_TYPE}_prediction.tsv"
        OUTPUT_FILE="${SITE_OUT_DIR}/${SAMPLE}_${MOD}_${FEATURE_TYPE}_site_level_prediction.tsv"

        echo "--------------------------------------"
        echo "Sample: ${SAMPLE}"
        echo "Modification: ${MOD}"
        echo "Feature type: ${FEATURE_TYPE}"
        echo "Input: ${INPUT_FILE}"
        echo "Output: ${OUTPUT_FILE}"

        if [ ! -f "$INPUT_FILE" ]; then
            echo "WARNING: input file missing, skipping:"
            echo "$INPUT_FILE"
            continue
        fi

        if [ -s "$OUTPUT_FILE" ]; then
            echo "Output already exists and is non-empty, skipping:"
            echo "$OUTPUT_FILE"
            continue
        fi

        sbatch \
            --job-name="site_${SAMPLE}_${MOD}" \
            tandemmod_read_to_site_one.slurm \
            "$INPUT_FILE" \
            "$OUTPUT_FILE"

        sleep 1

    done
done
