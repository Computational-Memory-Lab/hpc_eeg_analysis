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
# INPUT_FOLDER auto-resolution:
#   If INPUT_FOLDER is not a second-level folder directly, this wrapper will
#   try common layouts such as:
#     <run>/limo_second_level_<contrast_key>
#     <run>/epoch/limo_first_level/limo_second_level_<contrast_key>
#
# Optional title overrides:
#   TITLE_BASE          Base contrast text used to build default titles
#   CHANNEL_TIME_TITLE  Full custom channel-time title
#   TOPOPLOT_TITLE      Full custom topoplot title
#   LR_TITLE            Full custom LR title
#   TOPOPLOT_LAYOUT     Topoplot layout type: 'grid' (default), 'zigzag', or 'line'
#   TOPOPLOT_STEP_MS    Topoplot interval in ms (default: 100)
# ============================================================
INPUT_FOLDER="${INPUT_FOLDER:-/home/devon7y/scratch/devon7y/limo_second_level_test_hits_vs_test_misses}"
OUTPUT_DIR="${OUTPUT_DIR:-/home/devon7y/scratch/devon7y/plots}"
TITLE_BASE="${TITLE_BASE:-}"
CHANNEL_TIME_TITLE="${CHANNEL_TIME_TITLE:-}"
TOPOPLOT_TITLE="${TOPOPLOT_TITLE:-}"
LR_TITLE="${LR_TITLE:-}"
TOPOPLOT_LAYOUT="${TOPOPLOT_LAYOUT:-grid}"
TOPOPLOT_STEP_MS="${TOPOPLOT_STEP_MS:-100}"

has_second_level_artifacts() {
  local dir="$1"
  [[ -d "$dir" ]] || return 1
  [[ -f "$dir/LIMO.mat" ]] || return 1
  compgen -G "$dir/paired_samples_ttest_parameter_*.mat" > /dev/null || return 1
  return 0
}

resolve_second_level_dir() {
  local raw="$1"
  local base parent
  base="$(basename "$raw")"
  parent="$(dirname "$raw")"
  local candidates=(
    "$raw"
    "$raw/$base"
    "$parent/$base"
    "$parent/epoch/limo_first_level/$base"
    "$parent/limo_first_level/$base"
  )
  local c
  for c in "${candidates[@]}"; do
    if has_second_level_artifacts "$c"; then
      echo "$c"
      return 0
    fi
  done
  return 1
}

if ! RESOLVED_INPUT_FOLDER="$(resolve_second_level_dir "$INPUT_FOLDER")"; then
  echo "ERROR: Could not resolve a valid second-level folder from INPUT_FOLDER=${INPUT_FOLDER}" >&2
  echo "Expected LIMO.mat and paired_samples_ttest_parameter_*.mat in the same folder." >&2
  echo "Tried common layouts such as <run>/limo_second_level_* and <run>/epoch/limo_first_level/limo_second_level_*." >&2
  exit 1
fi
INPUT_FOLDER="$RESOLVED_INPUT_FOLDER"
echo "Using second-level folder: ${INPUT_FOLDER}"

module load matlab/2023b

# Fir's MATLAB module sets JAVA_TOOL_OPTIONS=-Xmx2g by default, which is
# too small for large SVG scene serialization and can trigger Java GC/OOM.
export JAVA_TOOL_OPTIONS="-Xmx8g"

matlab -nodisplay -nosplash -r "\
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1'); \
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1/plugins/limo_tools-3.4/'); \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
hpc_limo_channel_time_plots('${INPUT_FOLDER}', '${OUTPUT_DIR}', getenv('TITLE_BASE'), getenv('CHANNEL_TIME_TITLE'), getenv('TOPOPLOT_TITLE'), getenv('LR_TITLE'), getenv('TOPOPLOT_LAYOUT'), getenv('TOPOPLOT_STEP_MS')); \
exit;"
