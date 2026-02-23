#!/bin/bash
# submit_pipeline.sh
# Submits the full EEG analysis pipeline as chained SLURM jobs.
# Stages 1-3 are submitted as subject-level SLURM arrays.
# Optionally submits an epoch->ERP branch after Stage 3.
# Edit RAW_INPUT, PLOTS, epoch logic, ERP settings, condition order, and contrasts below.

set -euo pipefail

# ============================================================
# CONFIGURE THESE PATHS
# ============================================================
RAW_INPUT="/home/devon7y/scratch/devon7y/mydata"
PLOTS="/home/devon7y/scratch/devon7y/plots"
SCRIPTS="/home/devon7y/scratch/devon7y/hpc_eeg_analysis"

# ============================================================
# ARRAY SETTINGS (Stages 1-3)
# Tune throttles based on cluster availability and storage load.
# ============================================================
ARRAY_THROTTLE_STAGE1="24"
ARRAY_THROTTLE_STAGE2="8"
ARRAY_THROTTLE_STAGE3="24"

# ============================================================
# EPOCH + DESIGN SETTINGS
# ============================================================
VOLTAGE_DIFF="20"
VOLTAGE_ABS="1000"
EPOCH_TRIGGERS_CSV="11,21,22"
EPOCH_GROUP_SPEC="SME:11;Test_Intact:21;Test_Recombined:22"
CONDITION_ORDER=""  # optional: comma-separated labels; empty = infer from trial_type

# ============================================================
# OPTIONAL ERP BRANCH SETTINGS (Stage 4A)
# ============================================================
# 1 = submit hpc_epoch_to_erp_plot after Stage 3
# 0 = skip ERP branch and run only the LIMO branch
RUN_EPOCH_ERP_BRANCH="0"
TRIAL_TYPES_CSV="Study_hits,Study_misses"
ERP_CHANNELS_CSV="21"
ERP_OUTPUT_DIR=""     # empty => <EPOCH>/erp_plots
ERP_FIGURE_TITLE=""   # empty => default title from hpc_epoch_to_erp_plot.m

# Contrast entries:
#   key|condition_label_1|condition_label_2|display_title
COMPARISONS=(
  "study_hits_vs_study_misses|Study_hits|Study_misses|Subsequent Memory Effect (Study Hits vs Misses)"
  "test_hits_vs_test_misses|Test_hits|Test_misses|Retrieval Success Effect (Test Hits vs Misses)"
  "test_hits_vs_correct_rejections|Test_hits|Correct_rejections|Old/New Effect (Test Hits vs Correct Rejections)"
)

# ============================================================
# DERIVED PATHS (no need to edit)
# ============================================================
PARENT_DIR="$(dirname "${RAW_INPUT}")"
INITIAL_SET="${PARENT_DIR}/initial_set"
BEHAVIORAL_SET="${INITIAL_SET}/behavioral_set"
INTERPOL="${BEHAVIORAL_SET}/interpol"
EPOCH="${INTERPOL}/epoch"
FIRST_LEVEL="${EPOCH}/limo_first_level"
if [[ -z "${ERP_OUTPUT_DIR}" ]]; then
  ERP_OUTPUT_DIR="${EPOCH}/erp_plots"
fi
MANIFEST_DIR="${PARENT_DIR}/pipeline_manifests"
mkdir -p "${MANIFEST_DIR}"

# Per-stage log folders (each under that stage's output location)
LOG_STAGE1_DIR="${INITIAL_SET}/logs/raw_to_set"
LOG_STAGE2_DIR="${INTERPOL}/logs/set_to_interpol"
LOG_STAGE3_DIR="${EPOCH}/logs/interpol_to_epoch"
LOG_STAGE4A_DIR="${ERP_OUTPUT_DIR}/logs/epoch_to_erp_plot"
LOG_STAGE4B_DIR="${FIRST_LEVEL}/logs/limo_first_level"
LOG_PLOTS_DIR="${PLOTS}/logs/channel_time_plots"

mkdir -p "${LOG_STAGE1_DIR}" "${LOG_STAGE2_DIR}" "${LOG_STAGE3_DIR}" "${LOG_STAGE4B_DIR}" "${LOG_PLOTS_DIR}"
if [[ "${RUN_EPOCH_ERP_BRANCH}" == "1" ]]; then
  mkdir -p "${LOG_STAGE4A_DIR}"
fi

if [[ "${RUN_EPOCH_ERP_BRANCH}" != "0" && "${RUN_EPOCH_ERP_BRANCH}" != "1" ]]; then
  echo "ERROR: RUN_EPOCH_ERP_BRANCH must be 0 or 1, got: ${RUN_EPOCH_ERP_BRANCH}" >&2
  exit 1
fi

if [[ ! -d "${RAW_INPUT}" ]]; then
  echo "ERROR: RAW_INPUT folder does not exist: ${RAW_INPUT}" >&2
  exit 1
fi

REQUIRED_SCRIPTS=(
  "hpc_raw_to_set.slurm"
  "hpc_set_to_interpol.slurm"
  "hpc_interpol_to_epoch.slurm"
  "hpc_limo_first_level.slurm"
  "hpc_limo_second_level.slurm"
  "hpc_limo_channel_time_plots.slurm"
)
if [[ "${RUN_EPOCH_ERP_BRANCH}" == "1" ]]; then
  REQUIRED_SCRIPTS+=("hpc_epoch_to_erp_plot.slurm")
fi

for REQUIRED_SCRIPT in "${REQUIRED_SCRIPTS[@]}"; do
  if [[ ! -f "${SCRIPTS}/${REQUIRED_SCRIPT}" ]]; then
    echo "ERROR: Missing SLURM script: ${SCRIPTS}/${REQUIRED_SCRIPT}" >&2
    exit 1
  fi
done

# ============================================================
# BUILD SUBJECT MANIFEST FROM RAW_INPUT/*.raw
# ============================================================
SUBJECTS_FILE="${MANIFEST_DIR}/subjects_$(date +%Y%m%d_%H%M%S).txt"
TMP_SUBJECTS="${SUBJECTS_FILE}.tmp"

shopt -s nullglob
RAW_FILES=( "${RAW_INPUT}"/*.raw )
shopt -u nullglob

if [[ ${#RAW_FILES[@]} -eq 0 ]]; then
  echo "ERROR: No .raw files found in RAW_INPUT: ${RAW_INPUT}" >&2
  exit 1
fi

: > "${TMP_SUBJECTS}"
for RAW_FILE in "${RAW_FILES[@]}"; do
  RAW_BASE="$(basename "${RAW_FILE}" .raw)"
  SID="$(echo "${RAW_BASE}" | grep -oE '[0-9]+' | head -n 1 || true)"
  if [[ -n "${SID}" ]]; then
    echo "${SID}" >> "${TMP_SUBJECTS}"
  else
    echo "WARNING: Could not extract subject ID from raw file, skipping: ${RAW_FILE}" >&2
  fi
done

if [[ ! -s "${TMP_SUBJECTS}" ]]; then
  echo "ERROR: No parseable subject IDs were found in ${RAW_INPUT}" >&2
  rm -f "${TMP_SUBJECTS}"
  exit 1
fi

sort -n -u "${TMP_SUBJECTS}" > "${SUBJECTS_FILE}"
rm -f "${TMP_SUBJECTS}"

NUM_SUBJECTS="$(wc -l < "${SUBJECTS_FILE}" | tr -d '[:space:]')"
if [[ "${NUM_SUBJECTS}" -lt 1 ]]; then
  echo "ERROR: Subject manifest is empty: ${SUBJECTS_FILE}" >&2
  exit 1
fi

echo "Subject manifest created: ${SUBJECTS_FILE}"
echo "Subjects detected:        ${NUM_SUBJECTS}"

# ============================================================
# STAGE 1: raw -> set + behavioral alignment (array)
# ============================================================
JOB1=$(sbatch --parsable \
  --array=1-${NUM_SUBJECTS}%${ARRAY_THROTTLE_STAGE1} \
  --output="${LOG_STAGE1_DIR}/%x_%A_%a.out" \
  --error="${LOG_STAGE1_DIR}/%x_%A_%a.err" \
  --export=ALL,INPUT_FOLDER="${RAW_INPUT}",SUBJECTS_FILE="${SUBJECTS_FILE}" \
  "${SCRIPTS}/hpc_raw_to_set.slurm")
echo "Submitted job1 array (raw_to_set):            ${JOB1}"

# ============================================================
# STAGE 2: behavioral set -> interpol (array)
# ============================================================
JOB2=$(sbatch --parsable \
  --dependency=afterok:${JOB1} \
  --array=1-${NUM_SUBJECTS}%${ARRAY_THROTTLE_STAGE2} \
  --output="${LOG_STAGE2_DIR}/%x_%A_%a.out" \
  --error="${LOG_STAGE2_DIR}/%x_%A_%a.err" \
  --export=ALL,INPUT_FOLDER="${BEHAVIORAL_SET}",SUBJECTS_FILE="${SUBJECTS_FILE}" \
  "${SCRIPTS}/hpc_set_to_interpol.slurm")
echo "Submitted job2 array (set_to_interpol):       ${JOB2}"

# ============================================================
# STAGE 3: interpol -> epoch (array)
# ============================================================
JOB3=$(sbatch --parsable \
  --dependency=afterok:${JOB2} \
  --array=1-${NUM_SUBJECTS}%${ARRAY_THROTTLE_STAGE3} \
  --output="${LOG_STAGE3_DIR}/%x_%A_%a.out" \
  --error="${LOG_STAGE3_DIR}/%x_%A_%a.err" \
  --export=ALL,INPUT_FOLDER="${INTERPOL}",SUBJECTS_FILE="${SUBJECTS_FILE}",VOLTAGE_DIFF="${VOLTAGE_DIFF}",VOLTAGE_ABS="${VOLTAGE_ABS}",EPOCH_TRIGGERS_CSV="${EPOCH_TRIGGERS_CSV}",EPOCH_GROUP_SPEC="${EPOCH_GROUP_SPEC}" \
  "${SCRIPTS}/hpc_interpol_to_epoch.slurm")
echo "Submitted job3 array (interpol_to_epoch):     ${JOB3}"

# ============================================================
# STAGE 4A (optional): epoch -> ERP grand average plot
# ============================================================
JOB4A=""
if [[ "${RUN_EPOCH_ERP_BRANCH}" == "1" ]]; then
  JOB4A=$(sbatch --parsable \
    --dependency=afterok:${JOB3} \
    --output="${LOG_STAGE4A_DIR}/%x_%j.out" \
    --error="${LOG_STAGE4A_DIR}/%x_%j.err" \
    --export=ALL,INPUT_FOLDER="${EPOCH}",TRIAL_TYPES_CSV="${TRIAL_TYPES_CSV}",CHANNELS_CSV="${ERP_CHANNELS_CSV}",OUTPUT_DIR="${ERP_OUTPUT_DIR}",FIGURE_TITLE="${ERP_FIGURE_TITLE}" \
    "${SCRIPTS}/hpc_epoch_to_erp_plot.slurm")
  echo "Submitted job4A (epoch_to_erp_plot):         ${JOB4A}"
fi

# ============================================================
# STAGE 4B: epoch -> limo first level
# ============================================================
JOB4B=$(sbatch --parsable \
  --dependency=afterok:${JOB3} \
  --output="${LOG_STAGE4B_DIR}/%x_%j.out" \
  --error="${LOG_STAGE4B_DIR}/%x_%j.err" \
  --export=ALL,INPUT_FOLDER="${EPOCH}",CONDITION_ORDER="${CONDITION_ORDER}" \
  "${SCRIPTS}/hpc_limo_first_level.slurm")
echo "Submitted job4B (limo_first_level):          ${JOB4B}"

# ============================================================
# STAGE 5: second-level contrasts (parallel after job4B)
# ============================================================
declare -A JOB5
for ITEM in "${COMPARISONS[@]}"; do
  IFS='|' read -r KEY C1 C2 TITLE <<< "${ITEM}"
  SECOND_LEVEL_DIR="${FIRST_LEVEL}/limo_second_level_${KEY}"
  LOG_STAGE5_DIR="${SECOND_LEVEL_DIR}/logs/limo_second_level"
  mkdir -p "${LOG_STAGE5_DIR}"

  JOB5["${KEY}"]=$(sbatch --parsable \
    --dependency=afterok:${JOB4B} \
    --job-name="hpc_second_level_${KEY}" \
    --output="${LOG_STAGE5_DIR}/%x_%j.out" \
    --error="${LOG_STAGE5_DIR}/%x_%j.err" \
    --export=ALL,INPUT_FOLDER="${FIRST_LEVEL}",C1_LABEL="${C1}",C2_LABEL="${C2}",CONTRAST_KEY="${KEY}" \
    "${SCRIPTS}/hpc_limo_second_level.slurm")

  echo "Submitted job5_${KEY} (second_level ${C1} vs ${C2}): ${JOB5[${KEY}]}"
done

# ============================================================
# STAGE 6: plots (one per contrast, each after its job5)
# ============================================================
for ITEM in "${COMPARISONS[@]}"; do
  IFS='|' read -r KEY C1 C2 TITLE <<< "${ITEM}"
  SECOND_LEVEL_DIR="${FIRST_LEVEL}/limo_second_level_${KEY}"
  LOG_STAGE6_DIR="${LOG_PLOTS_DIR}/${KEY}"
  mkdir -p "${LOG_STAGE6_DIR}"

  JOB6=$(sbatch --parsable \
    --dependency=afterok:${JOB5[${KEY}]} \
    --job-name="hpc_plots_${KEY}" \
    --output="${LOG_STAGE6_DIR}/%x_%j.out" \
    --error="${LOG_STAGE6_DIR}/%x_%j.err" \
    --export=ALL,INPUT_FOLDER="${SECOND_LEVEL_DIR}",OUTPUT_DIR="${PLOTS}" \
    "${SCRIPTS}/hpc_limo_channel_time_plots.slurm")

  echo "Submitted job6_${KEY} (plots):                  ${JOB6} (after ${JOB5[${KEY}]})"
done

echo ""
echo "All jobs submitted. Monitor with: squeue -u \$USER"
echo "Subject manifest: ${SUBJECTS_FILE}"
echo "ERP branch enabled: ${RUN_EPOCH_ERP_BRANCH}"
echo "Array throttles:"
echo "  stage1: ${ARRAY_THROTTLE_STAGE1}"
echo "  stage2: ${ARRAY_THROTTLE_STAGE2}"
echo "  stage3: ${ARRAY_THROTTLE_STAGE3}"
echo "Log folders:"
echo "  stage1: ${LOG_STAGE1_DIR}"
echo "  stage2: ${LOG_STAGE2_DIR}"
echo "  stage3: ${LOG_STAGE3_DIR}"
if [[ "${RUN_EPOCH_ERP_BRANCH}" == "1" ]]; then
  echo "  stage4A: ${LOG_STAGE4A_DIR}"
fi
echo "  stage4B: ${LOG_STAGE4B_DIR}"
echo "  stage5: ${FIRST_LEVEL}/limo_second_level_<key>/logs/limo_second_level"
echo "  stage6: ${LOG_PLOTS_DIR}"
echo "Derived folders:"
echo "  initial_set:    ${INITIAL_SET}"
echo "  behavioral_set: ${BEHAVIORAL_SET}"
echo "  interpol:       ${INTERPOL}"
echo "  epoch:          ${EPOCH}"
if [[ "${RUN_EPOCH_ERP_BRANCH}" == "1" ]]; then
  echo "  erp_plots:      ${ERP_OUTPUT_DIR}"
fi
echo "  first_level:    ${FIRST_LEVEL}"
