# EEG Analysis Pipeline

Scripts are located in `/home/devon7y/scratch/devon7y/hpc_eeg_analysis/`.

Each script is a MATLAB **function** that takes a folder as input and writes output into a deterministic subfolder path. Between runs, only SLURM arguments/paths need to change.

---

## Pipeline Overview

```
Raw folder (contains *.raw + *.log + *.eeglog)
        |
        v
hpc_raw_to_set(raw_input_folder)
        |  outputs (under parent of raw_input_folder):
        |    <parent>/initial_set/
        |    <parent>/initial_set/behavioral_set/
        v
hpc_set_to_interpol(behavioral_set_folder)
        |  output: behavioral_set_folder/interpol/
        v
hpc_interpol_to_epoch(interpol_folder)
        |  output: interpol_folder/epoch/
        v
hpc_limo_first_level(epoch_folder)
        |  output: epoch_folder/limo_first_level/
        v
hpc_limo_second_level(limo_first_level_folder, [p1 p2])
        |  output: limo_first_level_folder/limo_second_level_<p1>_<p2>/
        v
hpc_limo_channel_time_plots(limo_second_level_folder, output_dir)
           output: output_dir/<test_name>_channel_time_plot.png
                   output_dir/<test_name>_channel_time_plot_LR.png
                   output_dir/<test_name>_topoplots.png
```

---

## Scripts

### 1. `hpc_raw_to_set.m`

**Purpose:** Convert EGI `.raw` to initial `.set`, then import behavior/event timing into final processed `.set` files.

```matlab
hpc_raw_to_set(input_folder)
```

| Argument | Description |
|---|---|
| `input_folder` | Folder containing `.raw`, session `.log`, and `.eeglog` files |

**Output folders (created under the parent of `input_folder`):**
- `<parent_of_input_folder>/initial_set/`
- `<parent_of_input_folder>/initial_set/behavioral_set/`

**Output files:**
- `initial_set/<ID>.set` - raw EGI import output
- `initial_set/behavioral_set/processed_<ID>.set` - behavior-aligned dataset
- `initial_set/behavioral_set/behavioral_<ID>.log`
- `initial_set/behavioral_set/EEGevents_<ID>.txt`
- `initial_set/behavioral_set/alignment_parameters_S<ID>.mat`

**Processing steps:**
1. Parse subject ID from raw filename (e.g., `1007.raw`)
2. Convert raw to `.set` with `pop_readegi`
3. Parse behavioral events from session `.log`
4. Align behavior and EEG timing using `.eeglog` pulse stream
5. Import aligned events into EEGLAB dataset
6. Save processed set to `initial_set/behavioral_set/`

---

### 2. `hpc_set_to_interpol.m`

**Purpose:** Filter and clean behaviorally aligned continuous EEG data.

```matlab
hpc_set_to_interpol(input_folder)
```

| Argument | Description |
|---|---|
| `input_folder` | Folder containing processed behavior-aligned `*.set` files |

**Output folder:** `input_folder/interpol/`

**Output files:**
- `preprocessed_full_<ID>.set` - cleaned EEG dataset per subject
- `preprocessing_summary_<ID>.txt` - text summary log per subject

**Processing steps:**
1. Flatline removal (channels flat > 5 s)
2. High-pass filter (0.1 Hz)
3. Low-pass filter (50 Hz)
4. Line noise removal (60 Hz, 120 Hz - CleanLine)
5. Kurtosis-based channel rejection (Z = 2 SD)
6. Channel location assignment (256-ch HydroCel GSN)
7. RANSAC-based channel rejection (correlation threshold = 0.8)
8. ICA decomposition (Extended Infomax - runica)
9. IC classification (ICLabel)
10. IC rejection (Eye > 90%, Muscle > 90%, Heart > 90%)
11. ASR - Artifact Subspace Reconstruction (cutoff = 20 SD)
12. Bad window removal - **currently disabled**

---

### 3. `hpc_interpol_to_epoch.m`

**Purpose:** Epoch continuous data and reject artifact trials.

```matlab
hpc_interpol_to_epoch(input_folder)
hpc_interpol_to_epoch(input_folder, voltage_diff_threshold, voltage_abs_threshold)
```

| Argument | Description | Default |
|---|---|---|
| `input_folder` | The `interpol/` folder from step 2 | - |
| `voltage_diff_threshold` | Max voltage step between samples (uV) | `20` |
| `voltage_abs_threshold` | Max absolute voltage (uV) | `1000` |

**Output folder:** `input_folder/epoch/`

**Output files:**
- `<ID>_epoch.set` - epoched, artifact-rejected EEG dataset per subject

**Processing steps:**
1. Epoch all trials together (triggers: `11`, `21`, `22`; window: -100 to 1500 ms)
2. Baseline correction (-100 to 0 ms)
3. Assign `trial_type` based on trigger + accuracy
4. Group trials into artifact rejection groups: `SME`, `Test_Intact`, `Test_Recombined`
5. Detect bad trials per group using voltage thresholds
6. Remove all bad trials in a single step

---

### 4. `hpc_limo_first_level.m`

**Purpose:** Run LIMO EEG first-level (per-subject) GLM analysis.

```matlab
hpc_limo_first_level(input_folder)
```

| Argument | Description |
|---|---|
| `input_folder` | The `epoch/` folder from step 3 |

**Output folder:** `input_folder/limo_first_level/`

**Output files:**
- `PairAsso_Epoched.study` - EEGLAB STUDY file
- `derivatives/sub-*/` - per-subject LIMO outputs (including `Betas.mat`)

---

### 5. `hpc_limo_second_level.m`

**Purpose:** Run LIMO group-level paired t-test between two conditions.

```matlab
hpc_limo_second_level(input_folder, parameters)
```

| Argument | Description |
|---|---|
| `input_folder` | The `limo_first_level/` folder from step 4 |
| `parameters` | `[p1 p2]` condition indices to compare |

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

**Example calls:**

```matlab
hpc_limo_second_level('/path/to/limo_first_level', [1 2]);  % SME
hpc_limo_second_level('/path/to/limo_first_level', [3 4]);  % Retrieval success
hpc_limo_second_level('/path/to/limo_first_level', [3 5]);  % Old/New
```

---

### 6. `hpc_limo_channel_time_plots.m`

**Purpose:** Generate channel-time heatmaps and topoplot grids from second-level results.

```matlab
hpc_limo_channel_time_plots(input_folder, output_dir)
```

| Argument | Description |
|---|---|
| `input_folder` | A `limo_second_level_<p1>_<p2>/` folder from step 5 |
| `output_dir` | Folder where `.png` files are saved |

**Output files (saved to `output_dir`):**
- `<test_name>_channel_time_plot.png`
- `<test_name>_channel_time_plot_LR.png` (if likelihood data exists)
- `<test_name>_topoplots.png`

---

## Example: Full Run for One Dataset

```matlab
raw_folder = '/home/devon7y/scratch/devon7y/mydata';
plots = '/home/devon7y/scratch/devon7y/plots';

hpc_raw_to_set(raw_folder);

parent_folder = fileparts(raw_folder);
behavioral_set = fullfile(parent_folder, 'initial_set', 'behavioral_set');

hpc_set_to_interpol(behavioral_set);
hpc_interpol_to_epoch(fullfile(behavioral_set, 'interpol'));
hpc_limo_first_level(fullfile(behavioral_set, 'interpol', 'epoch'));

first_level = fullfile(behavioral_set, 'interpol', 'epoch', 'limo_first_level');

hpc_limo_second_level(first_level, [1 2]);
hpc_limo_second_level(first_level, [3 4]);
hpc_limo_second_level(first_level, [3 5]);

hpc_limo_channel_time_plots(fullfile(first_level, 'limo_second_level_1_2'), plots);
hpc_limo_channel_time_plots(fullfile(first_level, 'limo_second_level_3_4'), plots);
hpc_limo_channel_time_plots(fullfile(first_level, 'limo_second_level_3_5'), plots);
```

---

## Running on the HPC (SLURM)

All SLURM scripts are in `hpc_eeg_analysis/` and can be chained with dependencies.

### SLURM scripts

| Script | Time | CPUs | RAM | Purpose |
|---|---|---|---|---|
| `hpc_raw_to_set.slurm` | 04:00:00 | 16 | 64G | Raw conversion + behavioral import |
| `hpc_set_to_interpol.slurm` | 12:00:00 | 64 | 128G | Filtering/cleaning |
| `hpc_interpol_to_epoch.slurm` | 00:05:00 | 64 | 128G | Epoching + trial rejection |
| `hpc_limo_first_level.slurm` | 00:30:00 | 64 | 256G | First-level LIMO |
| `hpc_limo_second_level.slurm` | 01:00:00 | 64 | 128G | Second-level LIMO |
| `hpc_limo_channel_time_plots.slurm` | 00:15:00 | 6 | 256G | Plot generation |

### Configuring each script

**`hpc_raw_to_set.slurm`**
```bash
INPUT_FOLDER="/home/devon7y/scratch/devon7y/mydata"
```

**`hpc_set_to_interpol.slurm`**
```bash
INPUT_FOLDER="/home/devon7y/scratch/devon7y/initial_set/behavioral_set"
```

**`hpc_interpol_to_epoch.slurm`**
```bash
INPUT_FOLDER="/home/devon7y/scratch/devon7y/initial_set/behavioral_set/interpol"
# VOLTAGE_DIFF=20
# VOLTAGE_ABS=1000
```

**`hpc_limo_first_level.slurm`**
```bash
INPUT_FOLDER="/home/devon7y/scratch/devon7y/initial_set/behavioral_set/interpol/epoch"
```

**`hpc_limo_second_level.slurm`**
```bash
INPUT_FOLDER="/home/devon7y/scratch/devon7y/initial_set/behavioral_set/interpol/epoch/limo_first_level"
P1=3
P2=4
```

**`hpc_limo_channel_time_plots.slurm`**
```bash
INPUT_FOLDER="/home/devon7y/scratch/devon7y/initial_set/behavioral_set/interpol/epoch/limo_first_level/limo_second_level_3_4"
OUTPUT_DIR="/home/devon7y/scratch/devon7y/plots"
```

### Submitting the full chained pipeline

```bash
bash /home/devon7y/scratch/devon7y/hpc_eeg_analysis/submit_pipeline.sh
```

Dependency graph:

```
[job1] hpc_raw_to_set
            | afterok
[job2] hpc_set_to_interpol
            | afterok
[job3] hpc_interpol_to_epoch
            | afterok
[job4] hpc_limo_first_level
            | afterok (parallel)
     +------+------+--------------+
[job5a] 1_2    [job5b] 3_4    [job5c] 3_5
     | afterok      | afterok      | afterok
[job6a] plots  [job6b] plots  [job6c] plots
```

To change comparisons, edit `COMPARISONS` in `submit_pipeline.sh`:

```bash
declare -A COMPARISONS
COMPARISONS["1_2"]="Subsequent Memory Effect (Study Hits vs Misses)"
COMPARISONS["3_4"]="Retrieval Success Effect (Test Hits vs Misses)"
COMPARISONS["3_5"]="Old/New Effect (Test Hits vs Correct Rejections)"
```

Monitor jobs with:

```bash
squeue -u $USER
```
