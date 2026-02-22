#!/bin/bash
# submit_pipeline.sh
# Submits the full EEG analysis pipeline as chained SLURM jobs.
# Edit RAW_INPUT, PLOTS, and comparisons below, then run:
#   bash submit_pipeline.sh

set -euo pipefail

# ============================================================
# CONFIGURE THESE PATHS
# ============================================================
RAW_INPUT="/home/devon7y/scratch/devon7y/mydata"
PLOTS="/home/devon7y/scratch/devon7y/plots"
SCRIPTS="/home/devon7y/scratch/devon7y/hpc_eeg_analysis"

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
echo "Submitted job1 (raw_to_set):          ${JOB1}"

# ============================================================
# STAGE 2: behavioral set -> interpol
# ============================================================
JOB2=$(sbatch --parsable \
  --dependency=afterok:${JOB1} \
  --export=ALL,INPUT_FOLDER="${BEHAVIORAL_SET}" \
  "${SCRIPTS}/hpc_set_to_interpol.slurm")
echo "Submitted job2 (set_to_interpol):      ${JOB2}"

# ============================================================
# STAGE 3: interpol -> epoch
# ============================================================
JOB3=$(sbatch --parsable \
  --dependency=afterok:${JOB2} \
  --export=ALL,INPUT_FOLDER="${INTERPOL}" \
  "${SCRIPTS}/hpc_interpol_to_epoch.slurm")
echo "Submitted job3 (interpol_to_epoch):    ${JOB3}"

# ============================================================
# STAGE 4: epoch -> limo first level
# ============================================================
JOB4=$(sbatch --parsable \
  --dependency=afterok:${JOB3} \
  --export=ALL,INPUT_FOLDER="${EPOCH}" \
  "${SCRIPTS}/hpc_limo_first_level.slurm")
echo "Submitted job4 (limo_first_level):     ${JOB4}"

# ============================================================
# STAGE 5: second level comparisons (run in parallel after job4)
# Add or remove comparisons here as needed.
# ============================================================
declare -A COMPARISONS
COMPARISONS["1_2"]="Subsequent Memory Effect (Study Hits vs Misses)"
COMPARISONS["3_4"]="Retrieval Success Effect (Test Hits vs Misses)"
COMPARISONS["3_5"]="Old/New Effect (Test Hits vs Correct Rejections)"

declare -A JOB5
for KEY in "${!COMPARISONS[@]}"; do
  P1=$(echo "$KEY" | cut -d_ -f1)
  P2=$(echo "$KEY" | cut -d_ -f2)
  JOB5[$KEY]=$(sbatch --parsable \
    --dependency=afterok:${JOB4} \
    --job-name="hpc_second_level_${KEY}" \
    --export=ALL,INPUT_FOLDER="${FIRST_LEVEL}",P1=${P1},P2=${P2} \
    "${SCRIPTS}/hpc_limo_second_level.slurm")
  echo "Submitted job5_${KEY} (second_level):   ${JOB5[$KEY]}"
done

# ============================================================
# STAGE 6: plots (one per comparison, each after its job5)
# ============================================================
for KEY in "${!COMPARISONS[@]}"; do
  TITLE="${COMPARISONS[$KEY]}"
  SECOND_LEVEL_DIR="${FIRST_LEVEL}/limo_second_level_${KEY}"
  JOB6=$(sbatch --parsable \
    --dependency=afterok:${JOB5[$KEY]} \
    --job-name="hpc_plots_${KEY}" \
    --export=ALL,INPUT_FOLDER="${SECOND_LEVEL_DIR}",TEST_TITLE="${TITLE}",OUTPUT_DIR="${PLOTS}" \
    "${SCRIPTS}/hpc_limo_channel_time_plots.slurm")
  echo "Submitted job6_${KEY} (plots):          ${JOB6} (after ${JOB5[$KEY]})"
done

echo ""
echo "All jobs submitted. Monitor with: squeue -u \$USER"
echo "Derived folders:"
echo "  initial_set:    ${INITIAL_SET}"
echo "  behavioral_set: ${BEHAVIORAL_SET}"
echo "  interpol:       ${INTERPOL}"
echo "  epoch:          ${EPOCH}"
echo "  first_level:    ${FIRST_LEVEL}"
