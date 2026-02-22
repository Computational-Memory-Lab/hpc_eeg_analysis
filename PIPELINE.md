# EEG Analysis Pipeline

Scripts live in `/home/devon7y/scratch/devon7y/hpc_eeg_analysis/`.

This pipeline now supports general experiments by using:
- a precomputed `event_label` column in session logs
- configurable epoch trigger logic
- condition-name driven second-level contrasts
- subject-scoped execution for Stages 1-3 (for SLURM arrays)

---

## End-to-End Flow

```text
raw_input_folder (contains *.raw, session *.log, *.eeglog)
  -> hpc_raw_to_set(raw_input_folder)
      outputs in: <parent>/initial_set/
                  <parent>/initial_set/behavioral_set/

<parent>/initial_set/behavioral_set/
  -> hpc_set_to_interpol(behavioral_set_folder)
      output: behavioral_set_folder/interpol/

.../interpol/
  -> hpc_interpol_to_epoch(interpol_folder, ...)
      output: interpol_folder/epoch/

.../epoch/
  -> hpc_limo_first_level(epoch_folder, condition_order)
      output: epoch_folder/limo_first_level/

.../limo_first_level/
  -> hpc_limo_second_level(first_level_folder, contrast, output_tag)
      output: first_level_folder/limo_second_level_<output_tag>/

.../limo_second_level_<output_tag>/
  -> hpc_limo_channel_time_plots(second_level_folder, output_dir)
```

---

## Session Log Requirements

`hpc_raw_to_set.m` expects each session `.log` to have a tab-delimited header row.

Required columns:
- `timestamp` (aliases: `unixtime`, `unix_time`)
- `trigger_code` (aliases: `trigger`, `type`, `event_code`)
- `event_label` (aliases: `trial_type`, `condition`)

Optional columns:
- `block`
- `trial`
- `word1_id` (aliases: `w1id`, `word1`, `word1id`)
- `word2_id` (aliases: `w2id`, `word2`, `word2id`)
- `target`
- `response`
- `accuracy` (alias: `acc`)
- `rt` (alias: `rt_ms`)
- `include_in_analysis` (alias: `include`)

Notes:
- `event_label` must already be computed by the experiment code.
- Epoching can still be based on `trigger_code`, independent of `event_label`.
- Detailed migration instructions from your old session-log format:
  - `SESSION_LOG_MIGRATION.md`

---

## Script Interfaces

### 1) Raw -> Behavioral Set

```matlab
hpc_raw_to_set(input_folder)
hpc_raw_to_set(input_folder, subject_filter)
```

Input:
- folder with `.raw`, session `.log`, `.eeglog`
- optional `subject_filter` (numeric ID) to process one subject only
- `.log` and `.eeglog` filenames should include the subject ID as a standalone numeric token
  - examples that match subject `12`: `session_12.log`, `12.eeglog`, `sub-12_session.log`
  - avoids ambiguous matches like subject `12` accidentally matching subject `112`

Outputs (under parent of `input_folder`):
- `initial_set/<ID>.set`
- `initial_set/behavioral_set/processed_<ID>.set`
- `initial_set/behavioral_set/behavioral_<ID>.log`
- `initial_set/behavioral_set/EEGevents_<ID>.txt`
- `initial_set/behavioral_set/alignment_parameters_S<ID>.mat`

### 2) Behavioral Set -> Interpol

```matlab
hpc_set_to_interpol(input_folder)
hpc_set_to_interpol(input_folder, subject_filter)
```

Input:
- `initial_set/behavioral_set/`
- optional `subject_filter` (numeric ID) to process one subject only

Output:
- `initial_set/behavioral_set/interpol/`

### 3) Interpol -> Epoch

```matlab
hpc_interpol_to_epoch(input_folder)
hpc_interpol_to_epoch(input_folder, voltage_diff_threshold, voltage_abs_threshold)
hpc_interpol_to_epoch(input_folder, voltage_diff_threshold, voltage_abs_threshold, epoch_triggers, group_spec)
hpc_interpol_to_epoch(input_folder, voltage_diff_threshold, voltage_abs_threshold, epoch_triggers, group_spec, subject_filter)
```

Parameters:
- `epoch_triggers`: which trigger codes to epoch around (defines which trials are created), accepts cell/numeric/CSV
  - examples: `{'11','21','22'}` or `'11,21,22'`
- `group_spec`: how already-epoched trials are grouped for artifact rejection thresholds/counting
  - string form: `'SME:11;Test_Intact:21;Test_Recombined:22'`
  - cell form: `{'SME', {'11'}; 'Test', {'21','22'}}`
  - practical rule: group triggers should be a subset of `epoch_triggers`
- `subject_filter`: optional numeric ID to process one subject only

Behavior:
- epochs on `epoch_triggers` (`EEG.event.type`)
- assigns `trial_type` from imported `EEG.event.event_label`
- applies artifact rejection by groups defined in `group_spec`

Output:
- `interpol/epoch/<ID>_epoch.set`

### 4) LIMO First Level

```matlab
hpc_limo_first_level(input_folder)
hpc_limo_first_level(input_folder, condition_order)
```

`condition_order`:
- optional CSV/cell order for design mapping
- if omitted/empty, inferred from `trial_type`

Outputs:
- `epoch/limo_first_level/PairAsso_Epoched.study`
- `epoch/limo_first_level/derivatives/...`
- `epoch/limo_first_level/condition_order.mat`
- `epoch/limo_first_level/condition_order.txt`

### 5) LIMO Second Level

```matlab
hpc_limo_second_level(input_folder, [p1 p2])
hpc_limo_second_level(input_folder, {'ConditionA','ConditionB'})
hpc_limo_second_level(input_folder, contrast, output_tag)
```

Contrast modes:
- numeric indices `[p1 p2]` (backward compatible)
- condition labels `{'ConditionA','ConditionB'}` (preferred)

When labels are used, indices are resolved from `condition_order.mat/txt`.

Outputs:
- `limo_second_level_<output_tag>/`
- `contrast_metadata.txt`

### 6) Plot Generation

```matlab
hpc_limo_channel_time_plots(input_folder, output_dir)
```

Input:
- one second-level folder

Output:
- channel-time and topoplot `.png` files in `output_dir`

---

## SLURM Scripts

Each SLURM script accepts environment-variable overrides via `sbatch --export`.
By default, script headers write to `/home/devon7y/scratch/devon7y/logs/`.
For organized runs, override `--output` and `--error` at submission time.

- `hpc_raw_to_set.slurm`
  - `INPUT_FOLDER`
  - optional subject scope: `SUBJECT_ID`
  - array mode: `SUBJECTS_FILE` (one subject ID per line) + `--array`
- `hpc_set_to_interpol.slurm`
  - `INPUT_FOLDER`
  - optional subject scope: `SUBJECT_ID`
  - array mode: `SUBJECTS_FILE` + `--array`
- `hpc_interpol_to_epoch.slurm`
  - `INPUT_FOLDER`, `VOLTAGE_DIFF`, `VOLTAGE_ABS`, `EPOCH_TRIGGERS_CSV`, `EPOCH_GROUP_SPEC`
  - optional subject scope: `SUBJECT_ID`
  - array mode: `SUBJECTS_FILE` + `--array`
- `hpc_limo_first_level.slurm`
  - `INPUT_FOLDER`, `CONDITION_ORDER` (empty allowed)
- `hpc_limo_second_level.slurm`
  - label mode: `C1_LABEL`, `C2_LABEL`, `CONTRAST_KEY`
  - numeric mode: `P1`, `P2`, `CONTRAST_KEY`
- `hpc_limo_channel_time_plots.slurm`
  - `INPUT_FOLDER`, `OUTPUT_DIR`

---

## Full Pipeline Submission

Use:

```bash
bash /home/devon7y/scratch/devon7y/hpc_eeg_analysis/submit_pipeline.sh
```

`submit_pipeline.sh` now supports:
- configurable epoching (`EPOCH_TRIGGERS_CSV`, `EPOCH_GROUP_SPEC`)
- optional explicit `CONDITION_ORDER`
- label-based contrasts in `COMPARISONS` (`key|label1|label2|title`)
- subject-manifest generation from `RAW_INPUT/*.raw`
- Stage 1-3 submission as SLURM arrays (`--array`) with throttles:
  - `ARRAY_THROTTLE_STAGE1`
  - `ARRAY_THROTTLE_STAGE2`
  - `ARRAY_THROTTLE_STAGE3`
- stage-local log placement via `sbatch --output/--error`:
  - Stage 1 logs: `<initial_set>/logs/raw_to_set/`
  - Stage 2 logs: `<interpol>/logs/set_to_interpol/`
  - Stage 3 logs: `<epoch>/logs/interpol_to_epoch/`
  - Stage 4 logs: `<limo_first_level>/logs/limo_first_level/`
  - Stage 5 logs: `<limo_second_level_<key>>/logs/limo_second_level/`
  - Stage 6 logs: `<PLOTS>/logs/channel_time_plots/<key>/`

Monitor jobs:

```bash
squeue -u $USER
```
