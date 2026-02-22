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
#   INPUT_FOLDER/limo_second_level_<CONTRAST_KEY>/
# ============================================================
INPUT_FOLDER="${INPUT_FOLDER:-/home/devon7y/scratch/devon7y/initial_set/behavioral_set/interpol/epoch/limo_first_level}"
P1="${P1:-3}"
P2="${P2:-4}"
C1_LABEL="${C1_LABEL:-}"
C2_LABEL="${C2_LABEL:-}"
CONTRAST_KEY="${CONTRAST_KEY:-}"

if [[ -z "${CONTRAST_KEY}" ]]; then
  if [[ -n "${C1_LABEL}" && -n "${C2_LABEL}" ]]; then
    CONTRAST_KEY="$(echo "${C1_LABEL}_vs_${C2_LABEL}" | tr '[:upper:]' '[:lower:]' | sed -E 's/[^a-z0-9]+/_/g; s/_+/_/g; s/^_+|_+$//g')"
  else
    CONTRAST_KEY="${P1}_${P2}"
  fi
fi

module load matlab/2023b

if [[ -n "${C1_LABEL}" && -n "${C2_LABEL}" ]]; then
  matlab -nodisplay -nosplash -r "\
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1'); \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
hpc_limo_second_level('${INPUT_FOLDER}', {'${C1_LABEL}','${C2_LABEL}'}, '${CONTRAST_KEY}'); \
exit;"
else
  matlab -nodisplay -nosplash -r "\
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1'); \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
hpc_limo_second_level('${INPUT_FOLDER}', [${P1} ${P2}], '${CONTRAST_KEY}'); \
exit;"
fi
