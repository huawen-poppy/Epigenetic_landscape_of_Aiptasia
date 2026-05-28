#!/bin/bash --login
#SBATCH -N 1
#SBATCH --constraint=v100
#SBATCH --gpus=1
#SBATCH -J step5_predict
#SBATCH -o ./logs/step5_predict.%A_%a.out
#SBATCH -e ./logs/step5_predict.%A_%a.err
#SBATCH --time=36:00:00
#SBATCH --mem=700G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=0-8%8

set -eo pipefail
set -x

source ~/miniconda3/etc/profile.d/conda.sh
conda activate epigenetic
module load gcc
TANDEMMOD_DIR="/ibex/scratch/projects/c2101/aip_epitrans_analysis/TandemMod"
BASE_DIR="/ibex/scratch/projects/c2101/aip_epitrans_analysis/raw_data_pass"
OUT_DIR="${TANDEMMOD_DIR}/TandemMod_prediction"

mkdir -p "$OUT_DIR"

SAMPLES=("A1" "A2" "A3" "A4" "H1" "H2" "H3" "H4")
MODS=("m6A" "m1A" "A_I" "m5C" "hm5C" "m7G" "G_I" "psU")

FEATURE_TYPES=("DRACH" "NNANN" "NNANN" "NNCNN" "NNCNN" "NNGNN" "NNGNN" "NNTNN")

MODELS=(
    "m6A_train_on_rice_cDNA.pkl"
    "m1A_train_on_rice_cDNA.pkl"
    "A_I_classifier_transfered_from_m5C.pkl"
    "m5C_train_on_rice_cDNA.pkl"
    "hm5C_transfered_from_m5C.pkl"
    "m7G_transfered_from_m5C.pkl"
    "G_I_classifier_transfered_from_m5C.pkl"
    "psU_transfered_from_m5C.pkl"
)


N_MODS=${#MODS[@]}

SAMPLE_INDEX=$((SLURM_ARRAY_TASK_ID / N_MODS))
MOD_INDEX=$((SLURM_ARRAY_TASK_ID % N_MODS))

SAMPLE="${SAMPLES[$SAMPLE_INDEX]}"
MOD="${MODS[$MOD_INDEX]}"
FEATURE_TYPE="${FEATURE_TYPES[$MOD_INDEX]}"
MODEL_NAME="${MODELS[$MOD_INDEX]}"

MODEL_FILE="${TANDEMMOD_DIR}/models/${MODEL_NAME}"

if [ "$FEATURE_TYPE" = "DRACH" ]; then
    FEATURE_FILE="${BASE_DIR}/${SAMPLE}/${SAMPLE}_DRACH.feature.tsv"
else
    FEATURE_FILE="${BASE_DIR}/${SAMPLE}/${SAMPLE}_${FEATURE_TYPE}_from_chunk.tsv"
fi

PREDICT_RESULT="${OUT_DIR}/${SAMPLE}_${MOD}_${FEATURE_TYPE}_prediction.tsv"

echo "======================================"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Sample: ${SAMPLE}"
echo "Modification: ${MOD}"
echo "Feature type: ${FEATURE_TYPE}"
echo "Feature file: ${FEATURE_FILE}"
echo "Model file: ${MODEL_FILE}"
echo "Prediction output: ${PREDICT_RESULT}"
echo "Started at: $(date)"
echo "======================================"

if [ ! -f "$FEATURE_FILE" ]; then
    echo "ERROR: feature file not found: ${FEATURE_FILE}"
    exit 1
fi

if [ ! -f "$MODEL_FILE" ]; then
    echo "ERROR: model file not found: ${MODEL_FILE}"
    exit 1
fi

cd "$TANDEMMOD_DIR"

python scripts/TandemMod.py \
    --run_mode predict \
    --pretrained_model "$MODEL_FILE" \
    --feature_file "$FEATURE_FILE" \
    --predict_result "$PREDICT_RESULT"

echo "Finished at: $(date)"
ls -lh "$PREDICT_RESULT"
