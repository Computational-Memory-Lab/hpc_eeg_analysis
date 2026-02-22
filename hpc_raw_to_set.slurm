#!/bin/bash
#SBATCH --account=def-jcaplan
#SBATCH --job-name=hpc_raw_to_set
#SBATCH --output=/home/devon7y/scratch/devon7y/logs/%x_%j.out
#SBATCH --error=/home/devon7y/scratch/devon7y/logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

# ============================================================
# INPUT: folder containing subject .raw + session .log + .eeglog files
# Output folders are created under the parent of INPUT_FOLDER:
#   <parent_of_INPUT_FOLDER>/initial_set/
#   <parent_of_INPUT_FOLDER>/initial_set/behavioral_set/
# ============================================================
INPUT_FOLDER="/home/devon7y/scratch/devon7y/mydata"

module load matlab/2023b
module load StdEnv/2023

matlab -nodisplay -nosplash -r "\
restoredefaultpath; \
addpath('/home/devon7y/scratch/devon7y/eeglab2025.0.0'); \
eeglab nogui; \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
hpc_raw_to_set('${INPUT_FOLDER}'); \
exit;"
