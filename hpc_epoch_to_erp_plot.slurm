#!/bin/bash
#SBATCH --account=def-jcaplan
#SBATCH --job-name=hpc_epoch_to_erp_plot
#SBATCH --output=/home/devon7y/scratch/devon7y/logs/%x_%j.out
#SBATCH --error=/home/devon7y/scratch/devon7y/logs/%x_%j.err
# Default log path above is a fallback.
# Prefer passing --output/--error at sbatch time to stage-local logs folders.
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G

# ============================================================
# INPUT_FOLDER: epoch/ folder produced by hpc_interpol_to_epoch
# TRIAL_TYPES_CSV: trial_type labels to plot.
#                  Supports comma or semicolon separators.
#                  For sbatch --export, semicolon is safer because commas are
#                  also used as env separators.
# CHANNELS_CSV: optional comma-separated channel indices (default: 21)
# OUTPUT_DIR: output folder for ERP plot files (default: INPUT_FOLDER/erp_plots)
# FIGURE_TITLE: optional full custom figure title
# ============================================================
INPUT_FOLDER="${INPUT_FOLDER:-/home/devon7y/scratch/devon7y/initial_set/behavioral_set/interpol/epoch}"
TRIAL_TYPES_CSV="${TRIAL_TYPES_CSV:-Study_hits;Study_misses}"
CHANNELS_CSV="${CHANNELS_CSV:-21}"
OUTPUT_DIR="${OUTPUT_DIR:-${INPUT_FOLDER}/erp_plots}"
FIGURE_TITLE="${FIGURE_TITLE:-}"

module load matlab/2023b

matlab -nodisplay -nosplash -r "\
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1'); \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
hpc_epoch_to_erp_plot('${INPUT_FOLDER}', '${TRIAL_TYPES_CSV}', '${CHANNELS_CSV}', '${OUTPUT_DIR}', getenv('FIGURE_TITLE')); \
exit;"
