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
# INPUT: the interpol/ folder produced by hpc_raw_to_interpol
# Output will be written to: INPUT_FOLDER/epoch/
#
# Optional: set voltage thresholds (uncomment to override defaults)
#   Default voltage_diff_threshold = 20 uV
#   Default voltage_abs_threshold  = 1000 uV
# ============================================================
INPUT_FOLDER="/home/devon7y/scratch/devon7y/mydata/interpol"
# VOLTAGE_DIFF=20
# VOLTAGE_ABS=1000

module load matlab/2023b

matlab -nodisplay -nosplash -r "\
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1'); \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
hpc_interpol_to_epoch('${INPUT_FOLDER}'); \
exit;"

# To override thresholds, replace the last two lines above with:
# hpc_interpol_to_epoch('${INPUT_FOLDER}', ${VOLTAGE_DIFF}, ${VOLTAGE_ABS}); \
