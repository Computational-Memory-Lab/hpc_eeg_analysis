# Skill Workflow Notes

## HPC Skill

Use the Codex skill `fir-hpc-workflow` for tasks involving Alliance Canada Fir cluster operations.

Path:

`/Users/devon7y/.codex/skills/fir-hpc-workflow/SKILL.md`

Invoke this skill whenever work includes:

- SSH access via alias `fir`
- Passcode-based re-authentication handling
- rsync upload/download to `/home/devon7y/scratch/devon7y/`
- SLURM script creation, submission, monitoring, or cancellation
- GPU job requests (H100/MIG) on Fir

## EEG Pipeline Skill

Use the Codex skill `eeg-hpc-pipeline` whenever the user asks for EEG processing, epoching, first-level LIMO, second-level LIMO, or channel-time plot generation.

Path:

`/Users/devon7y/.codex/skills/eeg-hpc-pipeline/SKILL.md`

This includes requests such as:

- "Run LIMO on this interpol folder."
- "Do second-level Study_hits vs Study_misses and graph it."
- "I have raw files; run up to preprocessing."

## Role Separation

- `eeg-hpc-pipeline` handles EEG workflow logic: required inputs, start-stage detection, end-stage mapping, defaults/overrides, preflight checks, and partial pipeline orchestration.
- `fir-hpc-workflow` handles cluster mechanics: authentication, SSH/rsync, `sbatch` execution, queue/accounting monitoring, cancellation, and data transfer.
- For EEG jobs on Fir, use both skills together: EEG skill decides what to run, HPC skill executes it on cluster.
