#!/bin/bash
#SBATCH --account=def-jcaplan
#SBATCH --job-name=hpc_interpol_to_epoch
#SBATCH --output=/home/devon7y/scratch/devon7y/logs/%x_%j.out
#SBATCH --error=/home/devon7y/scratch/devon7y/logs/%x_%j.err
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
# ============================================================
INPUT_FOLDER="${INPUT_FOLDER:-/home/devon7y/scratch/devon7y/initial_set/behavioral_set/interpol}"
VOLTAGE_DIFF="${VOLTAGE_DIFF:-20}"
VOLTAGE_ABS="${VOLTAGE_ABS:-1000}"
EPOCH_TRIGGERS_CSV="${EPOCH_TRIGGERS_CSV:-11,21,22}"
EPOCH_GROUP_SPEC="${EPOCH_GROUP_SPEC:-SME:11;Test_Intact:21;Test_Recombined:22}"

module load matlab/2023b

matlab -nodisplay -nosplash -r "\
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1'); \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
hpc_interpol_to_epoch('${INPUT_FOLDER}', ${VOLTAGE_DIFF}, ${VOLTAGE_ABS}, '${EPOCH_TRIGGERS_CSV}', '${EPOCH_GROUP_SPEC}'); \
exit;"
