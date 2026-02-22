#!/bin/bash
#SBATCH --account=def-jcaplan
#SBATCH --job-name=hpc_interpol_to_epoch
#SBATCH --output=/home/devon7y/scratch/devon7y/logs/%x_%j.out
#SBATCH --error=/home/devon7y/scratch/devon7y/logs/%x_%j.err
# Default log path above is a fallback.
# Prefer passing --output/--error at sbatch time to stage-local logs folders.
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G

# ============================================================
# INPUT: the interpol/ folder produced by hpc_set_to_interpol
# Output will be written to: INPUT_FOLDER/epoch/
#
# Optional controls (all overrideable via sbatch --export):
#   VOLTAGE_DIFF, VOLTAGE_ABS
#   EPOCH_TRIGGERS_CSV      e.g., "11,21,22"
#   EPOCH_GROUP_SPEC        e.g., "SME:11;Test_Intact:21;Test_Recombined:22"
# Optional subject scoping (for SLURM arrays):
#   SUBJECT_ID=<numeric_id>
#   SUBJECTS_FILE=<path_with_one_subject_id_per_line> + --array
# ============================================================
INPUT_FOLDER="${INPUT_FOLDER:-/home/devon7y/scratch/devon7y/initial_set/behavioral_set/interpol}"
VOLTAGE_DIFF="${VOLTAGE_DIFF:-20}"
VOLTAGE_ABS="${VOLTAGE_ABS:-1000}"
EPOCH_TRIGGERS_CSV="${EPOCH_TRIGGERS_CSV:-11,21,22}"
EPOCH_GROUP_SPEC="${EPOCH_GROUP_SPEC:-SME:11;Test_Intact:21;Test_Recombined:22}"
SUBJECT_ID="${SUBJECT_ID:-}"
SUBJECTS_FILE="${SUBJECTS_FILE:-}"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  if [[ -z "${SUBJECTS_FILE}" ]]; then
    echo "ERROR: SUBJECTS_FILE is required when using --array for hpc_interpol_to_epoch.slurm" >&2
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
  MATLAB_CALL="hpc_interpol_to_epoch('${INPUT_FOLDER}', ${VOLTAGE_DIFF}, ${VOLTAGE_ABS}, '${EPOCH_TRIGGERS_CSV}', '${EPOCH_GROUP_SPEC}', ${SUBJECT_ID});"
  echo "Running Stage 3 for subject ${SUBJECT_ID}"
else
  MATLAB_CALL="hpc_interpol_to_epoch('${INPUT_FOLDER}', ${VOLTAGE_DIFF}, ${VOLTAGE_ABS}, '${EPOCH_TRIGGERS_CSV}', '${EPOCH_GROUP_SPEC}');"
  echo "Running Stage 3 for all subjects in ${INPUT_FOLDER}"
fi

module load matlab/2023b

matlab -nodisplay -nosplash -r "\
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1'); \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
${MATLAB_CALL} \
exit;"
