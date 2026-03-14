#!/bin/bash
#SBATCH --account=def-jcaplan
#SBATCH --job-name=hpc_limo_first_level
#SBATCH --output=/home/devon7y/scratch/devon7y/logs/%x_%j.out
#SBATCH --error=/home/devon7y/scratch/devon7y/logs/%x_%j.err
# Default log path above is a fallback.
# Prefer passing --output/--error at sbatch time to stage-local logs folders.
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G

# ============================================================
# INPUT: the epoch/ folder produced by hpc_interpol_to_epoch
# Output will be written to: <parent_of_INPUT_FOLDER>/limo_first_level/
#
# CONDITION_ORDER is a comma-separated list of trial_type labels in the
# order that should map to LIMO parameter indices.
# Leave empty to infer order from trial_type values present in datasets.
# SUBJECT_ID is intentionally unsupported for this aggregate stage.
# Re-run Stage 4B on the full epoch folder after any subject-scoped
# Stage 1-3 updates.
# ============================================================
INPUT_FOLDER="${INPUT_FOLDER:-/home/devon7y/scratch/devon7y/epoch}"
CONDITION_ORDER="${CONDITION_ORDER:-}"
SUBJECT_ID="${SUBJECT_ID:-}"

if [[ -n "${SUBJECT_ID}" ]]; then
  echo "ERROR: SUBJECT_ID is not supported for hpc_limo_first_level.slurm. Re-run Stage 4B on the full epoch folder." >&2
  exit 1
fi

module load matlab/2023b

matlab -batch "\
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1'); \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
hpc_limo_first_level('${INPUT_FOLDER}', '${CONDITION_ORDER}');"
