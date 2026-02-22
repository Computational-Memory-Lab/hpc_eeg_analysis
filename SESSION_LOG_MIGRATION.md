# Session Log Migration Guide

This document explains exactly what to change in your session `.log` format so it works with `hpc_raw_to_set.m`.

It is based on your current format in:
- `/Users/devon7y/VS_Code/PyEPL3/data/12/session_12.log`

## Current Format (What You Have Now)

Your current logs are tab-delimited but mixed-layout (no header row, variable column counts):

- 4-column administrative rows:
  - `timestamp    0    LIST    P1`
  - `timestamp    0    B       Logging Begins`
  - `timestamp    0    E       Logging Ends`
- 9-column study rows (pair presentation):
  - `timestamp    0    list_id    trial    word1    word1_id    word2    word2_id    ipi`
- 11-column ITEM test rows:
  - `timestamp    0    list_id    trial    ITEM    word    word_id    target    response    accuracy    rt`
- 13-column ASSOC test rows:
  - `timestamp    0    list_id    trial    ASSOC    word1    word1_id    word2    word2_id    target    response    accuracy    rt`

Problem: the new pipeline no longer derives condition labels from this mixed structure. It expects precomputed labels in a normalized event table.

## Target Format (What You Need)

`hpc_raw_to_set.m` requires a tab-delimited header and event rows with at least:

- Required columns:
  - `timestamp`
  - `trigger_code`
  - `event_label`
- Recommended columns:
  - `block`
  - `trial`
  - `word1_id`
  - `word2_id`
  - `target`
  - `response`
  - `accuracy`
  - `rt`
  - `include_in_analysis` (1/0)

You may include extra columns; unknown columns are ignored.

## Required Changes

1. Add a header row.
- First analyzable header line must include the required columns above.
- Keep tab-delimited format.

2. Add `trigger_code` to analyzable event rows.
- This is the code transferred into `EEG.event.type`.
- `hpc_interpol_to_epoch.m` uses this for epoching.
- Choose trigger codes consistent with your epoch settings (`epoch_triggers`).

3. Add `event_label` to analyzable event rows.
- This must be the final analysis condition label.
- The pipeline now uses this label to create `trial_type` and LIMO conditions.
- Do not leave `event_label` blank for analyzable rows.

4. Keep a consistent schema for analyzable rows.
- No per-row column shifting for analyzable events.
- If a value is not applicable, use empty or `-1` but keep the column position.

5. Decide how to handle administrative/non-analyzable rows.
- Preferred: do not include them in the analyzable table.
- If you keep them, set `event_label` blank (or `include_in_analysis=0`) so they are ignored.

6. Ensure filename matching remains unambiguous.
- Keep subject ID in `.log` and `.eeglog` filenames as standalone numeric token.
- Good examples for subject `12`:
  - `session_12.log`
  - `12.eeglog`
  - `sub-12_session.log`

## Suggested Canonical Header

Use this as your standard schema:

```text
timestamp	block	list_id	trial	phase	word1	word2	word1_id	word2_id	target	response	accuracy	rt	trigger_code	event_label	include_in_analysis
```

Notes:
- `list_id`, `phase`, `word1`, `word2` are optional metadata (ignored by parser but useful for debugging).
- `block` can be numeric or strings like `P1`; parser can extract numeric token.

## Mapping From Your Current Rows

Use this mapping when writing the new log format.

1. Study rows (current 9-column rows):
- Current fields:
  - `timestamp, block_like_list, trial, word1, word1_id, word2, word2_id, ipi`
- New row:
  - Fill `trigger_code` for study events (for your current defaults usually `11`).
  - Fill `event_label` with your precomputed study condition (example: `Study_hits` / `Study_misses`).
  - Set non-applicable `target/response/accuracy` as `-1` or blank.

2. ASSOC rows (current 13-column rows):
- Current fields already contain `target`, `response`, `accuracy`, `rt`.
- Add:
  - `trigger_code` (for your current defaults usually `21` or `22`).
  - `event_label` (example labels: `Test_hits`, `Test_misses`, etc., based on your own design logic).

3. ITEM rows (current 11-column rows):
- Current fields contain one word ID plus `target`, `response`, `accuracy`, `rt`.
- Add:
  - `trigger_code` (example: `31`/`32` if used in your design).
  - `event_label` (example: `Correct_rejections`, `False_alarms`, etc., if applicable).
- For `word2_id`, set `-1` or blank.

4. LIST/B/E rows:
- Exclude from analyzable table, or keep with blank `event_label` / `include_in_analysis=0`.

## Trigger and Epoch Consistency

Important:
- `hpc_interpol_to_epoch.m` epochs by `trigger_code`, not by `event_label`.
- `group_spec` groups epoched trials by trigger code for artifact rejection.

So:
- if you want a row to be epoched, its `trigger_code` must be in `epoch_triggers`.
- `group_spec` trigger sets should be subsets of `epoch_triggers`.

## Example (Before -> After)

Before (current style, no header):

```text
1769382285047	0	P1	2	ASSOC	JAWBONE	1424	HATCHING	1265	1	1	1	1391
```

After (normalized row):

```text
timestamp	block	list_id	trial	phase	word1	word2	word1_id	word2_id	target	response	accuracy	rt	trigger_code	event_label	include_in_analysis
1769382285047	0	P1	2	ASSOC	JAWBONE	HATCHING	1424	1265	1	1	1	1391	21	Test_hits	1
```

## Validation Checklist Before Running Pipeline

1. Header exists and includes:
- `timestamp`, `trigger_code`, `event_label`

2. Every analyzable row has:
- numeric `timestamp`
- numeric `trigger_code`
- non-empty `event_label`

3. Labels are exactly what you want in LIMO:
- spelling/case are stable (`Study_hits` is different from `study_hits`)

4. Trigger configuration matches epoch step:
- `epoch_triggers` contains desired trigger codes
- `group_spec` is aligned with those triggers

5. Subject file naming is unambiguous:
- one session `.log` and one `.eeglog` per subject match

## Minimal Rule Set (If You Want the Short Version)

If you only do four things, do these:

1. Add a header row.
2. Add `trigger_code` column.
3. Add `event_label` column (precomputed final condition name).
4. Keep analyzable rows in one fixed tabular schema.

