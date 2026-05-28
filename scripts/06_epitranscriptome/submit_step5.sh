#!/bin/bash

targets=('A1' 'A2' 'A3' 'H3' 'H2' 'H1')
pre='/ibex/scratch/projects/c2101/aip_epitrans_analysis/raw_data_pass/'

modifications=("m6A" "m1A" "m5C" "hm5C" "m7G" "psU" "A_I_classifier" "G_I_classifier")

declare -A MODEL_MAP
MODEL_MAP["m6A"]="m6A_train_on_rice_cDNA.pkl"
MODEL_MAP["m1A"]="m1A_train_on_rice_cDNA.pkl"
MODEL_MAP["m5C"]="m5C_train_on_rice_cDNA.pkl"
MODEL_MAP["hm5C"]="hm5C_transfered_from_m5C.pkl"
MODEL_MAP["m7G"]="m7G_transfered_from_m5C.pkl"
MODEL_MAP["psU"]="psU_transfered_from_m5C.pkl"
MODEL_MAP["A_I_classifier"]="A_I_classifier_transfered_from_m5C.pkl"
MODEL_MAP["G_I_classifier"]="G_I_classifier_transfered_from_m5C.pkl"

# Modification -> feature motif suffix
declare -A MOTIF_MAP
MOTIF_MAP["m6A"]="DRACH"
MOTIF_MAP["m1A"]="NNANN"
MOTIF_MAP["m5C"]="NNCNN"
MOTIF_MAP["hm5C"]="NNCNN"
MOTIF_MAP["m7G"]="NNGNN"
MOTIF_MAP["psU"]="NNTNN"
MOTIF_MAP["A_I_classifier"]="NNANN"
MOTIF_MAP["G_I_classifier"]="NNGNN"

PRE="/ibex/scratch/projects/c2101/aip_epitrans_analysis/raw_data_pass"
TANDEMMOD_DIR="/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod"
MODEL_DIR="${TANDEMMOD_DIR}/models"
SLURM_SCRIPT="${TANDEMMOD_DIR}/step5_prediction.slurm"

for target in "${targets[@]}"; do
    SAMPLE_DIR="${PRE}/${target}"
    OUT_DIR="${TANDEMMOD_DIR}/TandemMod_prediction"
    mkdir -p "$OUT_DIR"

    for mod in "${!MODEL_MAP[@]}"; do
        motif="${MOTIF_MAP[$mod]}"
        model="${MODEL_DIR}/${MODEL_MAP[$mod]}"
        feature_file="${SAMPLE_DIR}/${target}_${motif}.feature.tsv"
        predict_result="${OUT_DIR}/${target}_${mod}.prediction.tsv"

        if [[ ! -f "$model" ]]; then
            echo "WARNING: model not found, skipping: $model"
            continue
        fi

        if [[ ! -f "$feature_file" ]]; then
            echo "WARNING: feature file not found, skipping: $feature_file"
            continue
        fi

        echo "Submitting prediction:"
        echo "  sample:  $target"
        echo "  mod:     $mod"
        echo "  motif:   $motif"
        echo "  model:   $model"
        echo "  feature: $feature_file"
        echo "  output:  $predict_result"

        #sbatch "$SLURM_SCRIPT" "$model" "$feature_file" "$predict_result"
    done
done

