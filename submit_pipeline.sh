#!/bin/bash
# submit_pipeline.sh
# Submits the full EEG analysis pipeline as chained SLURM jobs.
# Edit RAW_INPUT, PLOTS, epoch logic, condition order, and contrasts below.

set -euo pipefail

# ============================================================
# CONFIGURE THESE PATHS
# ============================================================
RAW_INPUT="/home/devon7y/scratch/devon7y/mydata"
PLOTS="/home/devon7y/scratch/devon7y/plots"
SCRIPTS="/home/devon7y/scratch/devon7y/hpc_eeg_analysis"

# ============================================================
# EPOCH + DESIGN SETTINGS
# ============================================================
VOLTAGE_DIFF="20"
VOLTAGE_ABS="1000"
EPOCH_TRIGGERS_CSV="11,21,22"
EPOCH_GROUP_SPEC="SME:11;Test_Intact:21;Test_Recombined:22"
CONDITION_ORDER=""  # optional: comma-separated labels; empty = infer from trial_type

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

# ============================================================
# STAGE 1: raw -> set + behavioral alignment
# ============================================================
JOB1=$(sbatch --parsable \
  --export=ALL,INPUT_FOLDER="${RAW_INPUT}" \
  "${SCRIPTS}/hpc_raw_to_set.slurm")
echo "Submitted job1 (raw_to_set):                  ${JOB1}"

# ============================================================
# STAGE 2: behavioral set -> interpol
# ============================================================
JOB2=$(sbatch --parsable \
  --dependency=afterok:${JOB1} \
  --export=ALL,INPUT_FOLDER="${BEHAVIORAL_SET}" \
  "${SCRIPTS}/hpc_set_to_interpol.slurm")
echo "Submitted job2 (set_to_interpol):             ${JOB2}"

# ============================================================
# STAGE 3: interpol -> epoch
# ============================================================
JOB3=$(sbatch --parsable \
  --dependency=afterok:${JOB2} \
  --export=ALL,INPUT_FOLDER="${INTERPOL}",VOLTAGE_DIFF="${VOLTAGE_DIFF}",VOLTAGE_ABS="${VOLTAGE_ABS}",EPOCH_TRIGGERS_CSV="${EPOCH_TRIGGERS_CSV}",EPOCH_GROUP_SPEC="${EPOCH_GROUP_SPEC}" \
  "${SCRIPTS}/hpc_interpol_to_epoch.slurm")
echo "Submitted job3 (interpol_to_epoch):           ${JOB3}"

# ============================================================
# STAGE 4: epoch -> limo first level
# ============================================================
JOB4=$(sbatch --parsable \
  --dependency=afterok:${JOB3} \
  --export=ALL,INPUT_FOLDER="${EPOCH}",CONDITION_ORDER="${CONDITION_ORDER}" \
  "${SCRIPTS}/hpc_limo_first_level.slurm")
echo "Submitted job4 (limo_first_level):            ${JOB4}"

# ============================================================
# STAGE 5: second-level contrasts (parallel after job4)
# ============================================================
declare -A JOB5
for ITEM in "${COMPARISONS[@]}"; do
  IFS='|' read -r KEY C1 C2 TITLE <<< "${ITEM}"

  JOB5["${KEY}"]=$(sbatch --parsable \
    --dependency=afterok:${JOB4} \
    --job-name="hpc_second_level_${KEY}" \
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

  JOB6=$(sbatch --parsable \
    --dependency=afterok:${JOB5[${KEY}]} \
    --job-name="hpc_plots_${KEY}" \
    --export=ALL,INPUT_FOLDER="${SECOND_LEVEL_DIR}",OUTPUT_DIR="${PLOTS}" \
    "${SCRIPTS}/hpc_limo_channel_time_plots.slurm")

  echo "Submitted job6_${KEY} (plots):                  ${JOB6} (after ${JOB5[${KEY}]})"
done

echo ""
echo "All jobs submitted. Monitor with: squeue -u \$USER"
echo "Derived folders:"
echo "  initial_set:    ${INITIAL_SET}"
echo "  behavioral_set: ${BEHAVIORAL_SET}"
echo "  interpol:       ${INTERPOL}"
echo "  epoch:          ${EPOCH}"
echo "  first_level:    ${FIRST_LEVEL}"
