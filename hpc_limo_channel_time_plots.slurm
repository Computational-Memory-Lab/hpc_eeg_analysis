#!/bin/bash
#SBATCH --account=def-jcaplan
#SBATCH --job-name=hpc_limo_channel_time_plots
#SBATCH --output=/home/devon7y/scratch/devon7y/logs/%x_%j.out
#SBATCH --error=/home/devon7y/scratch/devon7y/logs/%x_%j.err
# Default log path above is a fallback.
# Prefer passing --output/--error at sbatch time to stage-local logs folders.
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
# Observed MaxRSS for recent Step 6 runs is ~8.1G-11.4G.
# Keep ~2x headroom for larger layouts and export variability.
#SBATCH --mem=24G

# ============================================================
# INPUT_FOLDER: a limo_second_level_<contrast_key>/ folder
# OUTPUT_DIR:   folder where .png files will be saved
#
# Optional title overrides:
#   TITLE_BASE          Base contrast text used to build default titles
#   CHANNEL_TIME_TITLE  Full custom channel-time title
#   TOPOPLOT_TITLE      Full custom topoplot title
#   LR_TITLE            Full custom LR title
#   TOPOPLOT_LAYOUT     Topoplot layout type: 'grid' (default), 'zigzag', or 'line'
# ============================================================
INPUT_FOLDER="${INPUT_FOLDER:-/home/devon7y/scratch/devon7y/initial_set/behavioral_set/interpol/epoch/limo_first_level/limo_second_level_test_hits_vs_test_misses}"
OUTPUT_DIR="${OUTPUT_DIR:-/home/devon7y/scratch/devon7y/plots}"
TITLE_BASE="${TITLE_BASE:-}"
CHANNEL_TIME_TITLE="${CHANNEL_TIME_TITLE:-}"
TOPOPLOT_TITLE="${TOPOPLOT_TITLE:-}"
LR_TITLE="${LR_TITLE:-}"
TOPOPLOT_LAYOUT="${TOPOPLOT_LAYOUT:-grid}"

module load matlab/2023b

# Fir's MATLAB module sets JAVA_TOOL_OPTIONS=-Xmx2g by default, which is
# too small for large SVG scene serialization and can trigger Java GC/OOM.
export JAVA_TOOL_OPTIONS="-Xmx8g"

matlab -nodisplay -nosplash -r "\
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1'); \
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1/plugins/limo_tools-3.4/'); \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
hpc_limo_channel_time_plots('${INPUT_FOLDER}', '${OUTPUT_DIR}', getenv('TITLE_BASE'), getenv('CHANNEL_TIME_TITLE'), getenv('TOPOPLOT_TITLE'), getenv('LR_TITLE'), getenv('TOPOPLOT_LAYOUT')); \
exit;"
