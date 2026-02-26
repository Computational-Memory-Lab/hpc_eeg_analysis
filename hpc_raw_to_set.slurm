#!/bin/bash
#SBATCH --account=def-jcaplan
#SBATCH --job-name=hpc_raw_to_set
#SBATCH --output=/home/devon7y/scratch/devon7y/logs/%x_%j.out
#SBATCH --error=/home/devon7y/scratch/devon7y/logs/%x_%j.err
# Default log path above is a fallback.
# Prefer passing --output/--error at sbatch time to stage-local logs folders.
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

# ============================================================
# INPUT: folder containing subject .raw + session .log + .eeglog files
# Output folders are created under the parent of INPUT_FOLDER:
#   <parent_of_INPUT_FOLDER>/initial_set/
#   <parent_of_INPUT_FOLDER>/behavioral_set/
# Optional subject scoping (for SLURM arrays):
#   SUBJECT_ID=<numeric_id>
#   SUBJECTS_FILE=<path_with_one_subject_id_per_line> + --array
# ============================================================
INPUT_FOLDER="${INPUT_FOLDER:-/home/devon7y/scratch/devon7y/mydata}"
SUBJECT_ID="${SUBJECT_ID:-}"
SUBJECTS_FILE="${SUBJECTS_FILE:-}"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  if [[ -z "${SUBJECTS_FILE}" ]]; then
    echo "ERROR: SUBJECTS_FILE is required when using --array for hpc_raw_to_set.slurm" >&2
    exit 1
  fi
  if [[ ! -f "${SUBJECTS_FILE}" ]]; then
    echo "ERROR: SUBJECTS_FILE not found: ${SUBJECTS_FILE}" >&2
    exit 1
  fi
  SUBJECT_ID="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SUBJECTS_FILE}" | tr -d '[:space:]')"
fi

if [[ -n "${SUBJECT_ID}" ]] && [[ ! "${SUBJECT_ID}" =~ ^[0-9]+$ ]]; then
  echo "ERROR: SUBJECT_ID must be numeric, got: ${SUBJECT_ID}" >&2
  exit 1
fi

if [[ -n "${SUBJECT_ID}" ]]; then
  MATLAB_CALL="hpc_raw_to_set('${INPUT_FOLDER}', ${SUBJECT_ID});"
  echo "Running Stage 1 for subject ${SUBJECT_ID}"
else
  MATLAB_CALL="hpc_raw_to_set('${INPUT_FOLDER}');"
  echo "Running Stage 1 for all subjects in ${INPUT_FOLDER}"
fi

module load matlab/2023b
module load StdEnv/2023

matlab -nodisplay -nosplash -r "\
restoredefaultpath; \
addpath('/home/devon7y/scratch/devon7y/eeglab2025.0.0'); \
eeglab nogui; \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
${MATLAB_CALL} \
exit;"
