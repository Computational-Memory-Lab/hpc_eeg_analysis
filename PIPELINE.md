# EEG Analysis Pipeline

Scripts are located in `/home/devon7y/scratch/devon7y/hpc_eeg_analysis/`.

Each script is a MATLAB **function** that takes a folder as input and writes output into a subfolder of that folder. The only thing that needs to change between runs is the arguments passed in the SLURM script.

---

## Pipeline Overview

```
Raw .set files  (input_folder)
        │
        ▼
hpc_raw_to_interpol(input_folder)
        │  output: input_folder/interpol/
        ▼
hpc_interpol_to_epoch(interpol_folder)
        │  output: interpol_folder/epoch/
        ▼
hpc_limo_first_level(epoch_folder)
        │  output: epoch_folder/limo_first_level/
        ▼
hpc_limo_second_level(limo_first_level_folder, [p1 p2])
        │  output: limo_first_level_folder/limo_second_level_<p1>_<p2>/
        ▼
hpc_limo_channel_time_plots(limo_second_level_folder, output_dir)
           output: output_dir/<test_name>_channel_time_plot.png
                   output_dir/<test_name>_channel_time_plot_LR.png
                   output_dir/<test_name>_topoplots.png
```

---

## Scripts

### 1. `hpc_raw_to_interpol.m`

**Purpose:** Filter and clean raw continuous EEG data.

```matlab
hpc_raw_to_interpol(input_folder)
```

| Argument | Description |
|---|---|
| `input_folder` | Folder containing raw `*.set` files |

**Output folder:** `input_folder/interpol/`

**Output files:**
- `preprocessed_full_<ID>.set` — cleaned EEG dataset per subject
- `preprocessing_summary_<ID>.txt` — text log of processing summary per subject

**Processing steps:**
1. Flatline removal (channels flat > 5 s)
2. High-pass filter (0.1 Hz)
3. Low-pass filter (50 Hz)
4. Line noise removal (60 Hz, 120 Hz — CleanLine)
5. Kurtosis-based channel rejection (Z = 2 SD)
6. Channel location assignment (256-ch HydroCel GSN)
7. RANSAC-based channel rejection (correlation threshold = 0.8)
8. ICA decomposition (Extended Infomax — runica)
9. IC classification (ICLabel)
10. IC rejection (Eye > 90%, Muscle > 90%, Heart > 90%)
11. ASR — Artifact Subspace Reconstruction (cutoff = 20 SD)
12. Bad window removal — **currently disabled**

---

### 2. `hpc_interpol_to_epoch.m`

**Purpose:** Epoch the continuous data and reject artifact trials.

```matlab
hpc_interpol_to_epoch(input_folder)
hpc_interpol_to_epoch(input_folder, voltage_diff_threshold, voltage_abs_threshold)
```

| Argument | Description | Default |
|---|---|---|
| `input_folder` | The `interpol/` folder from step 1 | — |
| `voltage_diff_threshold` | Max voltage step between samples (μV) | `20` |
| `voltage_abs_threshold` | Max absolute voltage (μV) | `1000` |

**Output folder:** `input_folder/epoch/`

**Output files:**
- `<ID>_epoch.set` — epoched, artifact-rejected EEG dataset per subject

**Processing steps:**
1. Epoch all trials together (triggers: `11`, `21`, `22`; window: −100 to 1500 ms)
2. Baseline correction (−100 to 0 ms)
3. Assign `trial_type` field to each epoch based on trigger + accuracy:

| trial_type | Trigger | Accuracy |
|---|---|---|
| `Study_hits` | 11 | 1 |
| `Study_misses` | 11 | 0 |
| `Test_hits` | 21 | 1 |
| `Test_misses` | 21 | 0 |
| `Correct_rejections` | 22 | 1 |
| `False_alarms` | 22 | 0 |

4. Group trials into artifact rejection groups: `SME` (trigger 11), `Test_Intact` (21), `Test_Recombined` (22)
5. Detect bad trials per group using voltage thresholds
6. Remove all bad trials in a single step (subjects excluded if > 56 trials rejected in any group)

---

### 3. `hpc_limo_first_level.m`

**Purpose:** Run LIMO EEG 1st level (per-subject) GLM analysis.

```matlab
hpc_limo_first_level(input_folder)
```

| Argument | Description |
|---|---|
| `input_folder` | The `epoch/` folder from step 2 |

**Output folder:** `input_folder/limo_first_level/`

**Output files:**
- `PairAsso_Epoched.study` — EEGLAB STUDY file
- `derivatives/sub-*/` — per-subject LIMO results including `Betas.mat`

**Analysis:**
- Method: Weighted Least Squares (WLS)
- Measure: ERP (daterp)
- 6 conditions modelled as categorical variable `trial_type`

**Expected filename format for epoch files:** `<ID>_epoch.set`

---

### 4. `hpc_limo_second_level.m`

**Purpose:** Run LIMO group-level paired t-test between two conditions.

```matlab
hpc_limo_second_level(input_folder, parameters)
```

| Argument | Description |
|---|---|
| `input_folder` | The `limo_first_level/` folder from step 3 |
| `parameters` | `[p1 p2]` — condition indices to compare (e.g. `[3 4]`) |

**Condition index reference:**

| Index | Condition |
|---|---|
| 1 | Study_hits |
| 2 | Study_misses |
| 3 | Test_hits |
| 4 | Test_misses |
| 5 | Correct_rejections |
| 6 | False_alarms |

**Output folder:** `input_folder/limo_second_level_<p1>_<p2>/`

**Output files:** LIMO paired t-test results including:
- `paired_samples_ttest_parameter_<p1p2>.mat`
- `LIMO.mat`
- `tfce/`, `H0/` subdirectories

**Analysis:**
- Paired t-test via `limo_random_select`
- Full scalp analysis
- TFCE cluster correction
- 1000 bootstrap iterations

**Example calls:**

```matlab
% Retrieval Success Effect (Test hits vs misses)
hpc_limo_second_level('/path/to/limo_first_level', [3 4])
% → output: limo_first_level/limo_second_level_3_4/

% Subsequent Memory Effect (Study hits vs misses)
hpc_limo_second_level('/path/to/limo_first_level', [1 2])
% → output: limo_first_level/limo_second_level_1_2/

% Old/New Effect (Test hits vs Correct rejections)
hpc_limo_second_level('/path/to/limo_first_level', [3 5])
% → output: limo_first_level/limo_second_level_3_5/
```

---

### 5. `hpc_limo_channel_time_plots.m`

**Purpose:** Generate channel-time heatmap, likelihood ratio plot, and topoplot grid from a second-level result folder.

```matlab
hpc_limo_channel_time_plots(input_folder, output_dir)
```

| Argument | Description |
|---|---|
| `input_folder` | A `limo_second_level_<p1>_<p2>/` folder from step 4 |
| `output_dir` | Folder where `.png` files are saved |

The `parameter_num` and `test_name` are derived automatically from the `input_folder` name — no manual configuration needed.

**Output files (saved to `output_dir`):**
- `<test_name>_channel_time_plot.png` — signed −log₁₀(p) channel × time heatmap with TFCE correction
- `<test_name>_channel_time_plot_LR.png` — log₁₀(Bayes Factor) channel × time plot *(only if likelihood data exists)*
- `<test_name>_topoplots.png` — topoplot grid at 200 ms intervals, masked by TFCE significance

**Example call:**

```matlab
hpc_limo_channel_time_plots( ...
    '/path/to/limo_first_level/limo_second_level_3_4', ...
    '/path/to/plots')
```

---

## Example: Full Run for One Dataset

```matlab
base = '/home/devon7y/scratch/devon7y/mydata';
plots = '/home/devon7y/scratch/devon7y/plots';

hpc_raw_to_interpol(base);
hpc_interpol_to_epoch(fullfile(base, 'interpol'));
hpc_limo_first_level(fullfile(base, 'interpol', 'epoch'));

first_level = fullfile(base, 'interpol', 'epoch', 'limo_first_level');

% Run multiple second-level comparisons
hpc_limo_second_level(first_level, [1 2]);  % SME
hpc_limo_second_level(first_level, [3 4]);  % RSE
hpc_limo_second_level(first_level, [3 5]);  % ONE

% Generate plots for each
hpc_limo_channel_time_plots(fullfile(first_level, 'limo_second_level_1_2'), plots);
hpc_limo_channel_time_plots(fullfile(first_level, 'limo_second_level_3_4'), plots);
hpc_limo_channel_time_plots(fullfile(first_level, 'limo_second_level_3_5'), plots);
```

---

## Running on the HPC (SLURM)

All SLURM scripts are in `hpc_eeg_analysis/`. Each script corresponds to one pipeline stage. They use chained job dependencies so each stage only starts after the previous one succeeds.

### SLURM scripts

| Script | Copied from | Time | CPUs | RAM |
|---|---|---|---|---|
| `hpc_raw_to_interpol.slurm` | `devon_preprocessing.slurm` | 12:00:00 | 64 | 128G |
| `hpc_interpol_to_epoch.slurm` | `Interpol_to_Epoch.slurm` | 00:05:00 | 64 | 128G |
| `hpc_limo_first_level.slurm` | `tamari_first_level.slurm` | 00:30:00 | 64 | 256G |
| `hpc_limo_second_level.slurm` | `combined_shits_smisses.slurm` | 01:00:00 | 64 | 128G |
| `hpc_limo_channel_time_plots.slurm` | `new_ct_plots.slurm` | 00:15:00 | 6 | 256G |

### Configuring a script

Each script has a short variables block at the top — that is the only section you need to edit:

**`hpc_raw_to_interpol.slurm`**
```bash
INPUT_FOLDER="/home/devon7y/scratch/devon7y/mydata"
```

**`hpc_interpol_to_epoch.slurm`**
```bash
INPUT_FOLDER="/home/devon7y/scratch/devon7y/mydata/interpol"
# Optional — uncomment to override defaults (20 uV diff, 1000 uV abs):
# VOLTAGE_DIFF=20
# VOLTAGE_ABS=1000
```

**`hpc_limo_first_level.slurm`**
```bash
INPUT_FOLDER="/home/devon7y/scratch/devon7y/mydata/interpol/epoch"
```

**`hpc_limo_second_level.slurm`**
```bash
INPUT_FOLDER="/home/devon7y/scratch/devon7y/mydata/interpol/epoch/limo_first_level"
P1=3   # condition index 1
P2=4   # condition index 2
```

**`hpc_limo_channel_time_plots.slurm`**
```bash
INPUT_FOLDER="/home/devon7y/scratch/devon7y/mydata/interpol/epoch/limo_first_level/limo_second_level_3_4"
OUTPUT_DIR="/home/devon7y/scratch/devon7y/plots"
```

### Submitting the full pipeline

`submit_pipeline.sh` submits all jobs at once with automatic chaining. Edit the paths and comparisons at the top, then run it once:

```bash
bash /home/devon7y/scratch/devon7y/hpc_eeg_analysis/submit_pipeline.sh
```

The script chains jobs using `--dependency=afterok`, and fans out stage 4 and 5 in parallel — one job per condition comparison:

```
[job1] hpc_raw_to_interpol
            │ afterok
[job2] hpc_interpol_to_epoch
            │ afterok
[job3] hpc_limo_first_level
            │ afterok (all three run in parallel)
     ┌──────┴──────┬──────────────┐
[job4a] 1_2    [job4b] 3_4    [job4c] 3_5
     │ afterok      │ afterok      │ afterok
[job5a] plots  [job5b] plots  [job5c] plots
```

To add or remove comparisons, edit the `COMPARISONS` array in `submit_pipeline.sh`:

```bash
declare -A COMPARISONS
COMPARISONS["1_2"]="Subsequent Memory Effect (Study Hits vs Misses)"
COMPARISONS["3_4"]="Retrieval Success Effect (Test Hits vs Misses)"
COMPARISONS["3_5"]="Old/New Effect (Test Hits vs Correct Rejections)"
```

Monitor submitted jobs with:
```bash
squeue -u $USER
```
