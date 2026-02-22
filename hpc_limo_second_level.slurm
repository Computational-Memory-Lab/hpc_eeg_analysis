#!/bin/bash
#SBATCH --account=def-jcaplan
#SBATCH --job-name=hpc_limo_second_level
#SBATCH --output=/home/devon7y/scratch/devon7y/logs/%x_%j.out
#SBATCH --error=/home/devon7y/scratch/devon7y/logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G

# ============================================================
# INPUT: the limo_first_level/ folder produced by hpc_limo_first_level
# P1, P2: condition indices to compare
#   1=Study_hits  2=Study_misses  3=Test_hits
#   4=Test_misses 5=Correct_rejections 6=False_alarms
# Output will be written to: INPUT_FOLDER/limo_second_level_<P1>_<P2>/
# ============================================================
INPUT_FOLDER="/home/devon7y/scratch/devon7y/mydata/interpol/epoch/limo_first_level"
P1=3
P2=4

module load matlab/2023b

matlab -nodisplay -nosplash -r "\
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1'); \
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis'); \
hpc_limo_second_level('${INPUT_FOLDER}', [${P1} ${P2}]); \
exit;"
