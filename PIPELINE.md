# EEG Analysis Pipeline

Scripts live in `/home/devon7y/scratch/devon7y/hpc_eeg_analysis/`.

This pipeline now supports general experiments by using:
- a precomputed `event_label` column in session logs
- configurable epoch trigger logic
- an optional post-epoch ERP plotting branch
- condition-name driven second-level contrasts
- subject-scoped execution for Stages 1-3 (for SLURM arrays)

---

## End-to-End Flow

```text
raw_input_folder (contains *.raw, session *.log, *.eeglog)
  -> hpc_raw_to_set(raw_input_folder)
      outputs in: <parent>/initial_set/
                  <parent>/behavioral_set/

<parent>/behavioral_set/
  -> hpc_set_to_interpol(behavioral_set_folder)
      output: <parent>/interpol/

.../interpol/
  -> hpc_interpol_to_epoch(interpol_folder, ...)
      output: <parent>/epoch/

.../epoch/
  -> (Branch A) hpc_epoch_to_erp_plot(epoch_folder, trial_type_values, ...)
      output: <parent>/erp_plots/

.../epoch/
  -> (Branch B) hpc_limo_first_level(epoch_folder, condition_order)
      output: <parent>/limo_first_level/

.../limo_first_level/
  -> hpc_limo_second_level(first_level_folder, contrast, output_tag)
      output: <parent>/limo_second_level_<output_tag>/

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
- `behavioral_set/processed_<ID>.set`
- `behavioral_set/behavioral_<ID>.log`
- `behavioral_set/EEGevents_<ID>.txt`
- `behavioral_set/alignment_parameters_S<ID>.mat`

### 2) Behavioral Set -> Interpol

```matlab
hpc_set_to_interpol(input_folder)
hpc_set_to_interpol(input_folder, subject_filter)
```

Input:
- `behavioral_set/`
- optional `subject_filter` (numeric ID) to process one subject only

Output:
- `interpol/`

### 3) Interpol -> Epoch

```matlab
hpc_interpol_to_epoch(input_folder, epoch_window)
hpc_interpol_to_epoch(input_folder, epoch_window, voltage_diff_threshold, voltage_abs_threshold)
hpc_interpol_to_epoch(input_folder, epoch_window, voltage_diff_threshold, voltage_abs_threshold, epoch_triggers, group_spec)
hpc_interpol_to_epoch(input_folder, epoch_window, voltage_diff_threshold, voltage_abs_threshold, epoch_triggers, group_spec, subject_filter)
```

Parameters:
- `epoch_window`: REQUIRED `[start end]` window in seconds (no default)
  - examples: `[-0.1 1.5]`, `'-0.1,1.5'`, `'-0.1 1.5'`, `'-0.1;1.5'`
  - for `sbatch --export`, semicolon separator (`;`) is safer than comma
- `epoch_triggers`: which trigger codes to epoch around (defines which trials are created), accepts cell/numeric/CSV
  - examples: `{'11','21','22'}`, `'11,21,22'`, or `'11;21;22'`
  - for `sbatch --export`, semicolon separator (`;`) is safer than comma
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
- `epoch/<ID>_epoch.set`

### 4A) Epoch -> ERP Grand Average Plot (Branch)

```matlab
hpc_epoch_to_erp_plot(input_folder, trial_type_values, channels)
hpc_epoch_to_erp_plot(input_folder, trial_type_values, channels, output_dir)
hpc_epoch_to_erp_plot(input_folder, trial_type_values, channels, output_dir, figure_title)
hpc_epoch_to_erp_plot(input_folder, trial_type_values, channels, output_dir, figure_title, time_window_ms)
hpc_epoch_to_erp_plot(input_folder, trial_type_values, channels, output_dir, figure_title, time_window_ms, show_error_bars)
hpc_epoch_to_erp_plot(input_folder, trial_type_values, channels, output_dir, figure_title, time_window_ms, show_error_bars, plot_dimensions)
hpc_epoch_to_erp_plot(input_folder, trial_type_values, channels, output_dir, figure_title, time_window_ms, show_error_bars, plot_dimensions, x_axis_range_ms)
```

Inputs:
- `input_folder`: epoch folder containing `*_epoch.set`
- `trial_type_values`: labels to plot, accepts comma/semicolon list, string-array, or cell
  - examples: `{'Study_hits','Study_misses'}`, `'Study_hits,Study_misses'`, or `'Study_hits;Study_misses'`
- for `sbatch --export`, prefer semicolon separator (`;`) because commas also separate env vars
- required `channels`: channel indices to plot
  - accepts numeric vector or comma/semicolon/space-separated string (for `sbatch --export`, semicolon is safer)
  - if multiple channels are provided, the function produces one ERP figure per channel
- optional `output_dir`: default `<parent_of_epoch>/erp_plots`
- optional `figure_title`: full custom title
- optional `time_window_ms`: one or more `[start end]` windows
  - examples: `'300-500'`, `'300-500;600-800'`, or `[300 500; 600 800]`
  - when provided and exactly two conditions are requested:
    - the function computes subject-level mean voltage in each window for each condition
    - runs a paired t-test between the two condition window means (per window)
    - if only one subject has valid paired condition data, it falls back to a within-subject paired t-test across channel-time samples in the selected window
    - highlights each selected window on the ERP plot
    - annotates per-window significance on the plot
- optional `show_error_bars`: default `true`
  - accepts logical values (`true/false`, `1/0`, `yes/no`, `on/off`)
  - when `true`, plots basic within-subject SEM shading per condition
  - for single-subject runs, SEM shading falls back to trial-level SEM per condition/channel/timepoint
- optional `plot_dimensions`: figure position `[left bottom width height]`
  - default: `[100 100 1200 700]`
  - accepts numeric vectors or parseable strings
  - for `sbatch --export`, space-delimited values are safest (commas may need escaping)
- optional `x_axis_range_ms`: x-axis range in milliseconds
  - default: full timeline
  - accepts:
    - explicit range: `[start end]`, `'start-end'`, `'start,end'`, or `'start end'`
    - scalar duration in ms: for example `800` means from data start to data start + 800 ms
  - for `sbatch --export`, hyphen/space-delimited values are safest (commas may need escaping)
  - requested values are clamped to the available data range

Outputs:
- if one channel is requested:
  - `erp_plots/grand_average_erp_<trial_types>.png`
  - `erp_plots/grand_average_erp_<trial_types>.svg`
  - `erp_plots/grand_average_erp_<trial_types>.mat`
  - `erp_plots/grand_average_erp_<trial_types>_stats_<timestamp>.txt`
- if multiple channels are requested:
  - `erp_plots/grand_average_erp_<trial_types>_E<channel>.png`
  - `erp_plots/grand_average_erp_<trial_types>_E<channel>.svg`
  - `erp_plots/grand_average_erp_<trial_types>_E<channel>.mat`
  - `erp_plots/grand_average_erp_<trial_types>_E<channel>_stats_<timestamp>.txt`

Single-subject parameter-means variant (same plotting controls):

```matlab
hpc_epoch_to_erp_plot_single_subject(input_folder, trial_type_values, channels)
hpc_epoch_to_erp_plot_single_subject(input_folder, trial_type_values, channels, output_dir)
hpc_epoch_to_erp_plot_single_subject(input_folder, trial_type_values, channels, output_dir, figure_title)
hpc_epoch_to_erp_plot_single_subject(input_folder, trial_type_values, channels, output_dir, figure_title, time_window_ms)
hpc_epoch_to_erp_plot_single_subject(input_folder, trial_type_values, channels, output_dir, figure_title, time_window_ms, condition_order_source)
hpc_epoch_to_erp_plot_single_subject(input_folder, trial_type_values, channels, output_dir, figure_title, time_window_ms, condition_order_source, show_error_bars)
hpc_epoch_to_erp_plot_single_subject(input_folder, trial_type_values, channels, output_dir, figure_title, time_window_ms, condition_order_source, show_error_bars, plot_dimensions)
hpc_epoch_to_erp_plot_single_subject(input_folder, trial_type_values, channels, output_dir, figure_title, time_window_ms, condition_order_source, show_error_bars, plot_dimensions, x_axis_range_ms)
```

### 4B) LIMO First Level

```matlab
hpc_limo_first_level(input_folder)
hpc_limo_first_level(input_folder, condition_order)
```

`condition_order`:
- optional CSV/cell order for design mapping
- if omitted/empty, inferred from `trial_type`

Outputs:
- `limo_first_level/PairAsso_Epoched.study`
- `limo_first_level/derivatives/...`
- `limo_first_level/condition_order.mat`
- `limo_first_level/condition_order.txt`

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
hpc_limo_channel_time_plots(input_folder, output_dir, title_base)
hpc_limo_channel_time_plots(input_folder, output_dir, title_base, channel_time_title, topoplot_title)
hpc_limo_channel_time_plots(input_folder, output_dir, title_base, channel_time_title, topoplot_title, lr_title)
hpc_limo_channel_time_plots(input_folder, output_dir, title_base, channel_time_title, topoplot_title, lr_title, topoplot_layout_type)
```

Input:
- one second-level folder
- optional title overrides:
  - `title_base`: base contrast title used to build default titles
  - `channel_time_title`: full custom channel-time title
  - `topoplot_title`: full custom topoplot title
  - `lr_title`: full custom likelihood-ratio title
- optional topoplot layout override:
  - `topoplot_layout_type`: `'grid'` (default), `'zigzag'`, or `'line'`
  - `grid`: near-square arrangement
  - `zigzag`: staggered left-to-right zigzag line
  - `line`: single horizontal row

Output:
- channel-time and topoplot `.png` files in `output_dir`
- channel-time and topoplot `.svg` files in `output_dir`

---

## Manual Submission Preflight

Before launching Stage 5 or Stage 6 manually, validate the actual upstream folder:

- Stage 5 input must contain `*.study` and `derivatives/`.
- Stage 6 input must contain `LIMO.mat` and `paired_samples_ttest_parameter_*.mat`.
- Do not assume folder layout from memory; check with `ls`/`find` first.
- If Stage 3 is launched with `sbatch --export`, prefer semicolon separators for `EPOCH_WINDOW` and `EPOCH_TRIGGERS_CSV`.

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
  - `INPUT_FOLDER`, `EPOCH_WINDOW` (required), `VOLTAGE_DIFF`, `VOLTAGE_ABS`, `EPOCH_TRIGGERS_CSV`, `EPOCH_GROUP_SPEC`
  - `EPOCH_TRIGGERS_CSV` accepts comma or semicolon separators (semicolon preferred for `sbatch --export`)
  - escaped commas from `sbatch --export` are normalized by the wrapper
  - optional subject scope: `SUBJECT_ID`
  - array mode: `SUBJECTS_FILE` + `--array`
- `hpc_epoch_to_erp_plot.slurm`
  - `INPUT_FOLDER`, `TRIAL_TYPES_CSV`, `CHANNELS_CSV`
    - `TRIAL_TYPES_CSV` supports comma/semicolon separators (semicolon preferred for `sbatch --export`)
  - optional: `OUTPUT_DIR`, `FIGURE_TITLE`, `TIME_WINDOW_MS`, `SHOW_ERROR_BARS`, `PLOT_DIMENSIONS`, `X_AXIS_RANGE_MS`
- `hpc_limo_first_level.slurm`
  - `INPUT_FOLDER`, `CONDITION_ORDER` (empty allowed)
- `hpc_limo_second_level.slurm`
  - label mode: `C1_LABEL`, `C2_LABEL`, `CONTRAST_KEY`
  - numeric mode: `P1`, `P2`, `CONTRAST_KEY`
  - wrapper validates/resolves `INPUT_FOLDER` to a real first-level folder containing `*.study` + `derivatives/`
    - supports common layouts: `<run>/limo_first_level` and `<run>/epoch/limo_first_level`
- `hpc_limo_channel_time_plots.slurm`
  - `INPUT_FOLDER`, `OUTPUT_DIR`
  - optional title overrides:
    - `TITLE_BASE`
    - `CHANNEL_TIME_TITLE`
    - `TOPOPLOT_TITLE`
    - `LR_TITLE`
  - optional topoplot layout:
    - `TOPOPLOT_LAYOUT` = `grid` (default), `zigzag`, or `line`
  - wrapper validates/resolves `INPUT_FOLDER` to a real second-level folder containing `LIMO.mat` + `paired_samples_ttest_parameter_*.mat`
    - supports common layouts: `<run>/limo_second_level_<key>` and `<run>/epoch/limo_first_level/limo_second_level_<key>`

---

## Full Pipeline Submission

Use:

```bash
bash /home/devon7y/scratch/devon7y/hpc_eeg_analysis/submit_pipeline.sh
```

`submit_pipeline.sh` now supports:
- configurable epoching (`EPOCH_WINDOW` required, plus `EPOCH_TRIGGERS_CSV`, `EPOCH_GROUP_SPEC`)
  - `EPOCH_TRIGGERS_CSV` default is semicolon-separated (`11;21;22`) to avoid `sbatch --export` CSV splitting
  - commas in Stage 3 env values are escaped and normalized by wrappers
- optional explicit `CONDITION_ORDER`
- optional Stage 4A ERP branch:
  - `RUN_EPOCH_ERP_BRANCH` (`0` or `1`)
  - `TRIAL_TYPES_CSV` (comma/semicolon list; semicolon preferred)
  - `ERP_CHANNELS_CSV` (required when `RUN_EPOCH_ERP_BRANCH=1`)
  - `ERP_OUTPUT_DIR` (empty => `erp_plots`)
  - `ERP_FIGURE_TITLE`
  - `ERP_TIME_WINDOW_MS` (optional window(s) for stats/highlight, example: `"300-500;600-800"`)
  - `ERP_SHOW_ERROR_BARS` (optional; default `true`; accepts true/false, 1/0, yes/no, on/off)
  - `ERP_PLOT_DIMENSIONS` (optional figure position; default: `"100 100 1200 700"`)
  - `ERP_X_AXIS_RANGE_MS` (optional x-axis range/duration in ms; examples: `"-100-1000"` or `"1100"`)
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
  - Stage 4A logs (when enabled): `<erp_plots>/logs/epoch_to_erp_plot/`
  - Stage 4B logs: `<limo_first_level>/logs/limo_first_level/`
  - Stage 5 logs: `<limo_second_level_<key>>/logs/limo_second_level/`
  - Stage 6 logs: `<PLOTS>/logs/channel_time_plots/<key>/`

Note:
- `submit_pipeline.sh` always runs the LIMO branch after Stage 3.
- It can also run the ERP branch in parallel when `RUN_EPOCH_ERP_BRANCH=1`.
- You can still run only the ERP branch directly with `hpc_epoch_to_erp_plot.slurm`.

Monitor jobs:

```bash
squeue -u $USER
```
