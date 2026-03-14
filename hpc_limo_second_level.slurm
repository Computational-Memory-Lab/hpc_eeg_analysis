#!/bin/bash
#SBATCH --account=def-jcaplan
#SBATCH --job-name=hpc_limo_second_level
#SBATCH --output=/home/devon7y/scratch/devon7y/logs/%x_%j.out
#SBATCH --error=/home/devon7y/scratch/devon7y/logs/%x_%j.err
# Default log path above is a fallback.
# Prefer passing --output/--error at sbatch time to stage-local logs folders.
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G

# ============================================================
# INPUT: the limo_first_level/ folder produced by hpc_limo_first_level
#
# You can drive contrast selection either by labels or by numeric indices:
# 1) Preferred label mode:
#    C1_LABEL, C2_LABEL, CONTRAST_KEY
# 2) Backward-compatible numeric mode:
#    P1, P2, CONTRAST_KEY (optional)
#
# Output will be written to:
#   <parent_of_INPUT_FOLDER>/limo_second_level_<CONTRAST_KEY>/
#
# SUBJECT_ID is intentionally unsupported for this aggregate stage.
# Re-run second-level analysis on the full first-level folder after any
# subject-scoped upstream updates.
#
# INPUT_FOLDER auto-resolution:
#   If INPUT_FOLDER is not a first-level folder directly, this wrapper will
#   try common layouts such as:
#     <run>/limo_first_level
#     <run>/epoch/limo_first_level
# ============================================================
INPUT_FOLDER="${INPUT_FOLDER:-/home/devon7y/scratch/devon7y/limo_first_level}"
P1="${P1:-3}"
P2="${P2:-4}"
C1_LABEL="${C1_LABEL:-}"
C2_LABEL="${C2_LABEL:-}"
CONTRAST_KEY="${CONTRAST_KEY:-}"
SUBJECT_ID="${SUBJECT_ID:-}"

if [[ -n "${SUBJECT_ID}" ]]; then
  echo "ERROR: SUBJECT_ID is not supported for hpc_limo_second_level.slurm. Re-run second-level analysis on the full first-level folder." >&2
  exit 1
fi

has_first_level_artifacts() {
  local dir="$1"
  [[ -d "$dir" ]] || return 1
  compgen -G "$dir/*.study" > /dev/null || return 1
  [[ -d "$dir/derivatives" ]] || return 1
  return 0
}

resolve_first_level_dir() {
  local raw="$1"
  local parent
  parent="$(dirname "$raw")"
  local candidates=(
    "$raw"
    "$raw/limo_first_level"
    "$raw/epoch/limo_first_level"
    "$parent/limo_first_level"
    "$parent/epoch/limo_first_level"
  )
  local c
  for c in "${candidates[@]}"; do
    if has_first_level_artifacts "$c"; then
      echo "$c"
      return 0
    fi
  done
  return 1
}

if ! RESOLVED_INPUT_FOLDER="$(resolve_first_level_dir "$INPUT_FOLDER")"; then
  echo "ERROR: Could not resolve a valid first-level folder from INPUT_FOLDER=${INPUT_FOLDER}" >&2
  echo "Expected a folder containing *.study and derivatives/." >&2
  echo "Tried common layouts such as <run>/limo_first_level and <run>/epoch/limo_first_level." >&2
  exit 1
fi
INPUT_FOLDER="$RESOLVED_INPUT_FOLDER"
echo "Using first-level folder: ${INPUT_FOLDER}"

if [[ -z "${CONTRAST_KEY}" ]]; then
  if [[ -n "${C1_LABEL}" && -n "${C2_LABEL}" ]]; then
    CONTRAST_KEY="$(echo "${C1_LABEL}_vs_${C2_LABEL}" | tr '[:upper:]' '[:lower:]' | sed -E 's/[^a-z0-9]+/_/g; s/_+/_/g; s/^_+|_+$//g')"
  else
    CONTRAST_KEY="${P1}_${P2}"
  fi
fi

module load matlab/2023b

if [[ -n "${C1_LABEL}" && -n "${C2_LABEL}" ]]; then
  matlab -batch "\
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1'); \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
hpc_limo_second_level('${INPUT_FOLDER}', {'${C1_LABEL}','${C2_LABEL}'}, '${CONTRAST_KEY}');"
else
  matlab -batch "\
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1'); \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
hpc_limo_second_level('${INPUT_FOLDER}', [${P1} ${P2}], '${CONTRAST_KEY}');"
fi
