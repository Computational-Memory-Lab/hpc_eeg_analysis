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
# CHANNELS_CSV: required channel indices (comma/semicolon/space separated)
#               examples: "21", "21,22", or "21;22"
# OUTPUT_DIR: output folder for ERP plot files (default: <parent_of_INPUT_FOLDER>/erp_plots)
# SUBJECT_ID: optional numeric subject ID to limit ERP plotting to one
#             epoched dataset within INPUT_FOLDER
# FIGURE_TITLE: optional full custom figure title
# TIME_WINDOW_MS: optional time window(s) for stats/highlight
#                 examples: "300-500" or "300-500;600-800"
# SHOW_ERROR_BARS: optional boolean for within-subject SEM shading
#                  accepted: true/false, 1/0, yes/no, on/off (default: true)
# PLOT_DIMENSIONS: optional figure size/position [left bottom width height]
#                  example: "100 100 1200 700" (default: [100 100 1200 700])
# X_AXIS_RANGE_MS: optional x-axis range or duration in ms
#                  examples: "-100-1000", "-100,1000", or "1100" (duration)
# ============================================================
INPUT_FOLDER="${INPUT_FOLDER:-/home/devon7y/scratch/devon7y/epoch}"
TRIAL_TYPES_CSV="${TRIAL_TYPES_CSV:-Study_hits;Study_misses}"
CHANNELS_CSV="${CHANNELS_CSV:-}"
DEFAULT_ERP_OUTPUT_DIR="$(dirname "${INPUT_FOLDER}")/erp_plots"
OUTPUT_DIR="${OUTPUT_DIR:-${DEFAULT_ERP_OUTPUT_DIR}}"
FIGURE_TITLE="${FIGURE_TITLE:-}"
TIME_WINDOW_MS="${TIME_WINDOW_MS:-}"
SHOW_ERROR_BARS="${SHOW_ERROR_BARS:-true}"
PLOT_DIMENSIONS="${PLOT_DIMENSIONS:-100 100 1200 700}"
X_AXIS_RANGE_MS="${X_AXIS_RANGE_MS:-}"
SUBJECT_ID="${SUBJECT_ID:-}"

if [[ -z "${CHANNELS_CSV//[[:space:]]/}" ]]; then
  echo "ERROR: CHANNELS_CSV is required (example: 21 or 21,22)." >&2
  exit 1
fi

if [[ -n "${SUBJECT_ID}" ]] && [[ ! "${SUBJECT_ID}" =~ ^[0-9]+$ ]]; then
  echo "ERROR: SUBJECT_ID must be numeric, got: ${SUBJECT_ID}" >&2
  exit 1
fi

module load matlab/2023b

if [[ -n "${SUBJECT_ID}" ]]; then
  matlab -batch "\
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1'); \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
hpc_epoch_to_erp_plot('${INPUT_FOLDER}', '${TRIAL_TYPES_CSV}', '${CHANNELS_CSV}', '${OUTPUT_DIR}', getenv('FIGURE_TITLE'), getenv('TIME_WINDOW_MS'), getenv('SHOW_ERROR_BARS'), getenv('PLOT_DIMENSIONS'), getenv('X_AXIS_RANGE_MS'), ${SUBJECT_ID});"
else
  matlab -batch "\
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1'); \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
hpc_epoch_to_erp_plot('${INPUT_FOLDER}', '${TRIAL_TYPES_CSV}', '${CHANNELS_CSV}', '${OUTPUT_DIR}', getenv('FIGURE_TITLE'), getenv('TIME_WINDOW_MS'), getenv('SHOW_ERROR_BARS'), getenv('PLOT_DIMENSIONS'), getenv('X_AXIS_RANGE_MS'));"
fi
