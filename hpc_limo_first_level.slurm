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
# Output will be written to: INPUT_FOLDER/limo_first_level/
#
# CONDITION_ORDER is a comma-separated list of trial_type labels in the
# order that should map to LIMO parameter indices.
# Leave empty to infer order from trial_type values present in datasets.
# ============================================================
INPUT_FOLDER="${INPUT_FOLDER:-/home/devon7y/scratch/devon7y/initial_set/behavioral_set/interpol/epoch}"
CONDITION_ORDER="${CONDITION_ORDER:-}"

module load matlab/2023b

matlab -nodisplay -nosplash -r "\
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1'); \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
hpc_limo_first_level('${INPUT_FOLDER}', '${CONDITION_ORDER}'); \
exit;"
