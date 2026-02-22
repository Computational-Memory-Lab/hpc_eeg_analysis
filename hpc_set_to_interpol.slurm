#!/bin/bash
#SBATCH --account=def-jcaplan
#SBATCH --job-name=hpc_set_to_interpol
#SBATCH --output=/home/devon7y/scratch/devon7y/logs/%x_%j.out
#SBATCH --error=/home/devon7y/scratch/devon7y/logs/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G

# ============================================================
# INPUT: folder containing behavioral .set files to preprocess
# Output will be written to: INPUT_FOLDER/interpol/
# ============================================================
INPUT_FOLDER="/home/devon7y/scratch/devon7y/initial_set/behavioral_set"

module load matlab/2023b
module load StdEnv/2023

matlab -nodisplay -nosplash -r "\
restoredefaultpath; \
addpath('/home/devon7y/scratch/devon7y/eeglab2025.0.0'); \
eeglab nogui; \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
hpc_set_to_interpol('${INPUT_FOLDER}'); \
exit;"
