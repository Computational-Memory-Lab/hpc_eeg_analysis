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
#SBATCH --mem=256G

# ============================================================
# INPUT_FOLDER: a limo_second_level_<contrast_key>/ folder
# OUTPUT_DIR:   folder where .png files will be saved
# ============================================================
INPUT_FOLDER="${INPUT_FOLDER:-/home/devon7y/scratch/devon7y/initial_set/behavioral_set/interpol/epoch/limo_first_level/limo_second_level_test_hits_vs_test_misses}"
OUTPUT_DIR="${OUTPUT_DIR:-/home/devon7y/scratch/devon7y/plots}"

module load matlab/2023b

matlab -nodisplay -nosplash -r "\
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1'); \
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1/plugins/limo_tools-3.4/'); \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
hpc_limo_channel_time_plots('${INPUT_FOLDER}', '${OUTPUT_DIR}'); \
exit;"
