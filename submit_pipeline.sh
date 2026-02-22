#!/bin/bash
# submit_pipeline.sh
# Submits the full EEG analysis pipeline as chained SLURM jobs.
# Edit the BASE, PLOTS, and comparisons below, then run:
#   bash submit_pipeline.sh

# ============================================================
# CONFIGURE THESE PATHS
# ============================================================
BASE="/home/devon7y/scratch/devon7y/mydata"
PLOTS="/home/devon7y/scratch/devon7y/plots"
SCRIPTS="/home/devon7y/scratch/devon7y/hpc_eeg_analysis"

# ============================================================
# DERIVED PATHS (no need to edit)
# ============================================================
INTERPOL="${BASE}/interpol"
EPOCH="${INTERPOL}/epoch"
FIRST_LEVEL="${EPOCH}/limo_first_level"

# ============================================================
# STAGE 1: raw -> interpol
# ============================================================
JOB1=$(sbatch --parsable \
  --export=ALL,INPUT_FOLDER="${BASE}" \
  "${SCRIPTS}/hpc_raw_to_interpol.slurm")
echo "Submitted job1 (raw_to_interpol):   ${JOB1}"

# ============================================================
# STAGE 2: interpol -> epoch
# ============================================================
JOB2=$(sbatch --parsable \
  --dependency=afterok:${JOB1} \
  --export=ALL,INPUT_FOLDER="${INTERPOL}" \
  "${SCRIPTS}/hpc_interpol_to_epoch.slurm")
echo "Submitted job2 (interpol_to_epoch): ${JOB2}"

# ============================================================
# STAGE 3: epoch -> limo first level
# ============================================================
JOB3=$(sbatch --parsable \
  --dependency=afterok:${JOB2} \
  --export=ALL,INPUT_FOLDER="${EPOCH}" \
  "${SCRIPTS}/hpc_limo_first_level.slurm")
echo "Submitted job3 (limo_first_level):  ${JOB3}"

# ============================================================
# STAGE 4: second level comparisons (run in parallel after job3)
# Add or remove comparisons here as needed.
# ============================================================
declare -A COMPARISONS
COMPARISONS["1_2"]="Subsequent Memory Effect (Study Hits vs Misses)"
COMPARISONS["3_4"]="Retrieval Success Effect (Test Hits vs Misses)"
COMPARISONS["3_5"]="Old/New Effect (Test Hits vs Correct Rejections)"

declare -A JOB4
for KEY in "${!COMPARISONS[@]}"; do
  P1=$(echo $KEY | cut -d_ -f1)
  P2=$(echo $KEY | cut -d_ -f2)
  JOB4[$KEY]=$(sbatch --parsable \
    --dependency=afterok:${JOB3} \
    --job-name="hpc_second_level_${KEY}" \
    --export=ALL,INPUT_FOLDER="${FIRST_LEVEL}",P1=${P1},P2=${P2} \
    "${SCRIPTS}/hpc_limo_second_level.slurm")
  echo "Submitted job4_${KEY} (second_level): ${JOB4[$KEY]}"
done

# ============================================================
# STAGE 5: plots (one per comparison, each after its job4)
# ============================================================
for KEY in "${!COMPARISONS[@]}"; do
  TITLE="${COMPARISONS[$KEY]}"
  SECOND_LEVEL_DIR="${FIRST_LEVEL}/limo_second_level_${KEY}"
  sbatch \
    --dependency=afterok:${JOB4[$KEY]} \
    --job-name="hpc_plots_${KEY}" \
    --export=ALL,INPUT_FOLDER="${SECOND_LEVEL_DIR}",TEST_TITLE="${TITLE}",OUTPUT_DIR="${PLOTS}" \
    "${SCRIPTS}/hpc_limo_channel_time_plots.slurm"
  echo "Submitted job5_${KEY} (plots):         (after ${JOB4[$KEY]})"
done

echo ""
echo "All jobs submitted. Monitor with: squeue -u \$USER"
