# Skill Workflow Notes

Use `eeg-hpc-pipeline` for EEG workflow requests (preprocessing, epoching, first-level LIMO, second-level LIMO, and channel-time plots):

`/Users/devon7y/.codex/skills/eeg-hpc-pipeline/SKILL.md`

Use `fir-hpc-workflow` for Alliance Fir execution mechanics (SSH/passcode handling, rsync, sbatch/squeue/sacct/scancel, GPU request syntax):

`/Users/devon7y/.codex/skills/fir-hpc-workflow/SKILL.md`

Role separation:

- `eeg-hpc-pipeline` decides what stages to run, from which start stage, with which parameters and validations.
- `fir-hpc-workflow` performs the cluster operations needed to execute that plan.
- For EEG processing on Fir, apply both skills together.
