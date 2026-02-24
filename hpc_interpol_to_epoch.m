function hpc_interpol_to_epoch(input_folder, voltage_diff_threshold, voltage_abs_threshold, epoch_triggers, group_spec, subject_filter)
% HPC_INTERPOL_TO_EPOCH - Epoch EEG data and apply artifact rejection
%
% Usage:
%   hpc_interpol_to_epoch(input_folder)
%   hpc_interpol_to_epoch(input_folder, voltage_diff_threshold, voltage_abs_threshold)
%   hpc_interpol_to_epoch(input_folder, voltage_diff_threshold, voltage_abs_threshold, epoch_triggers, group_spec)
%   hpc_interpol_to_epoch(input_folder, voltage_diff_threshold, voltage_abs_threshold, epoch_triggers, group_spec, subject_filter)
%
% Inputs:
%   input_folder           - Path to folder containing .set files
%   voltage_diff_threshold - Max voltage difference threshold in uV (default: 20)
%   voltage_abs_threshold  - Max absolute voltage threshold in uV (default: 1000)
%   epoch_triggers         - Trigger codes for epoching (default: {'11','21','22'})
%                            Accepts cell array, numeric array, or CSV string.
%   group_spec             - Artifact rejection groups by trigger code.
%                            Default: 'SME:11;Test_Intact:21;Test_Recombined:22'
%                            Accepts:
%                              1) Nx2 cell, e.g. {'SME', {'11'}; 'Test', {'21','22'}}
%                              2) spec string, e.g. 'SME:11;Test:21,22'
%   subject_filter         - Optional numeric subject ID. When provided,
%                            process only that subject.
%
% Notes:
%   - Epoching uses trigger codes (EEG.event.type), not event_label.
%   - trial_type is assigned from EEG.event.event_label after epoching.

if nargin < 2 || isempty(voltage_diff_threshold)
    voltage_diff_threshold = 20;
end
if nargin < 3 || isempty(voltage_abs_threshold)
    voltage_abs_threshold = 1000;
end
if nargin < 4 || isempty(epoch_triggers)
    epoch_triggers = {'11', '21', '22'};
end
if nargin < 5 || isempty(group_spec)
    group_spec = 'SME:11;Test_Intact:21;Test_Recombined:22';
end
if nargin < 6
    subject_filter = [];
else
    subject_filter = parse_optional_subject_filter(subject_filter);
end

epoch_triggers = normalize_trigger_list(epoch_triggers);
group_def = parse_group_spec(group_spec);
condition_names = {'Study_hits', 'Study_misses', 'Test_hits', 'Test_misses', 'Correct_rejections', 'False_alarms'};

% Create output folder and summary file path.
epoch_folder = fullfile(input_folder, 'epoch');
if ~exist(epoch_folder, 'dir')
    mkdir(epoch_folder);
end
summary_txt_file = fullfile(epoch_folder, ...
    sprintf('Processing_Summary_%s.txt', datestr(now, 'yyyy-mm-dd_HH-MM-SS')));

fprintf('\n==============================================\n');
fprintf('  HPC INTERPOL -> EPOCH\n');
fprintf('==============================================\n');
fprintf('Input folder:                %s\n', input_folder);
fprintf('Voltage diff threshold (uV): %.1f\n', voltage_diff_threshold);
fprintf('Voltage abs threshold  (uV): %.1f\n', voltage_abs_threshold);
fprintf('Epoch triggers:              %s\n', strjoin(epoch_triggers, ','));
fprintf('Artifact groups:\n');
for g = 1:size(group_def, 1)
    fprintf('  - %s: %s\n', group_def{g, 1}, strjoin(group_def{g, 2}, ','));
end
fprintf('==============================================\n\n');

% Initialize EEGLAB
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1');
eeglab nogui;

% Find input files
set_files = dir(fullfile(input_folder, '*.set'));
if isempty(set_files)
    error('No .set files found in: %s', input_folder);
end

filtered_files = struct('name', {}, 'folder', {}, 'date', {}, 'bytes', {}, 'isdir', {}, 'datenum', {});
for i = 1:length(set_files)
    fname = set_files(i).name;
    if ~startsWith(fname, 'before_step9') && ~contains(fname, '_epoch')
        filtered_files(end+1) = set_files(i); %#ok<AGROW>
    end
end
if isempty(filtered_files)
    error('No valid input .set files found in: %s', input_folder);
end

file_info = struct('filename', {}, 'filepath', {}, 'subject_id', {});
for i = 1:length(filtered_files)
    fname = filtered_files(i).name;
    [~, basename, ~] = fileparts(fname);
    tokens = regexp(basename, '(\d+)', 'tokens');
    if ~isempty(tokens)
        sid = str2double(tokens{1}{1});
        if ~isempty(subject_filter) && sid ~= subject_filter
            continue;
        end
        file_info(end+1).filename = fname; %#ok<AGROW>
        file_info(end).filepath = filtered_files(i).folder;
        file_info(end).subject_id = sid;
    else
        fprintf('WARNING: Could not extract subject ID from "%s", skipping\n', fname);
    end
end
if isempty(file_info)
    if ~isempty(subject_filter)
        error('Subject %d not found in interpol folder: %s', subject_filter, input_folder);
    else
        error('No .set files with parseable subject IDs found in: %s', input_folder);
    end
end

fprintf('Files to process: %d\n', length(file_info));
if ~isempty(subject_filter)
    fprintf('Subject filter: %d\n', subject_filter);
end

% Processing settings
max_rejected_trials_per_group = 56;

% Legacy-compatible summary containers (replaces old Processing_Summary_*.mat content).
stats = initialize_epoch_stats(length(file_info), voltage_diff_threshold, ...
    voltage_abs_threshold, max_rejected_trials_per_group, condition_names);
results = initialize_epoch_results(length(file_info), condition_names);

for file_idx = 1:length(file_info)
    current_file = file_info(file_idx);
    s = current_file.subject_id;

    output_filename = sprintf('%d_epoch.set', s);
    output_filepath = fullfile(epoch_folder, output_filename);
    results.subject_id{file_idx} = sprintf('%d', s);

    fprintf('\n----------------------------------------------\n');
    fprintf('Subject %d (%d/%d)\n', s, file_idx, length(file_info));

    if exist(output_filepath, 'file')
        fprintf('Output already exists, skipping: %s\n', output_filename);
        results.skipped(file_idx) = true;
        results.processed(file_idx) = true;
        continue;
    end

    EEG = pop_loadset('filename', current_file.filename, 'filepath', current_file.filepath);

    if ~isfield(EEG.event, 'event_label')
        error(['Subject %d is missing EEG.event.event_label. ' ...
               'Ensure hpc_raw_to_set imported event_label from the session log.'], s);
    end

    if isfield(EEG.event, 'trial_type')
        EEG.event = rmfield(EEG.event, 'trial_type');
    end

    % STEP 1: epoch by trigger codes
    EEG = pop_epoch(EEG, epoch_triggers, [-0.1 1.5], ...
        'newname', sprintf('Subject %d - all epoched', s), 'epochinfo', 'yes');
    EEG = pop_rmbase(EEG, [-100 0]);
    results.original_trials(file_idx) = EEG.trials;
    fprintf('Epoched %d trials\n', EEG.trials);

    % STEP 2: set trial_type directly from event_label
    [EEG, cond_counts] = assign_trial_type_from_event_label(EEG, epoch_triggers);
    subject_condition_original = condition_struct_from_map(cond_counts, condition_names);
    fprintf('Condition counts (pre-rejection):\n');
    for c = 1:numel(condition_names)
        cond_name = condition_names{c};
        fprintf('  %-30s %d\n', cond_name, subject_condition_original.(cond_name));
    end

    % STEP 3: assign epochs to artifact groups
    trial_indices_by_group = assign_epochs_to_groups(EEG, epoch_triggers, group_def);

    % STEP 4: detect bad trials per group
    all_bad_trials_absolute = [];
    exclude_subject = false;
    exclusion_reason = '';
    all_channels = 1:EEG.nbchan;
    subject_stats = struct('condition1_total', 0, 'condition2_total', 0, ...
        'SME_original', 0, 'SME_removed', 0, 'SME_final', 0, ...
        'Test_Intact_original', 0, 'Test_Intact_removed', 0, 'Test_Intact_final', 0, ...
        'Test_Recombined_original', 0, 'Test_Recombined_removed', 0, 'Test_Recombined_final', 0);

    for g = 1:size(group_def, 1)
        group_name = group_def{g, 1};
        group_trial_indices = trial_indices_by_group{g};

        fprintf('Group %s: %d trials\n', group_name, numel(group_trial_indices));
        if isempty(group_trial_indices)
            continue;
        end

        group_original_trials = numel(group_trial_indices);
        if strcmp(group_name, 'SME')
            subject_stats.SME_original = group_original_trials;
        elseif strcmp(group_name, 'Test_Intact')
            subject_stats.Test_Intact_original = group_original_trials;
        elseif strcmp(group_name, 'Test_Recombined')
            subject_stats.Test_Recombined_original = group_original_trials;
        end

        EEG_group = pop_select(EEG, 'trial', group_trial_indices);
        [~, bad_trials_relative, cond1_bad, cond2_bad] = get_bad_trials_v2(EEG_group, ...
            voltage_diff_threshold, voltage_abs_threshold, all_channels);
        fprintf('  Bad trials in %s: %d (diff=%d, abs=%d)\n', ...
            group_name, numel(bad_trials_relative), numel(cond1_bad), numel(cond2_bad));

        subject_stats.condition1_total = subject_stats.condition1_total + numel(cond1_bad);
        subject_stats.condition2_total = subject_stats.condition2_total + numel(cond2_bad);

        if numel(bad_trials_relative) > max_rejected_trials_per_group
            exclude_subject = true;
            exclusion_reason = sprintf('Too many trials rejected in %s group (%d > %d)', ...
                group_name, numel(bad_trials_relative), max_rejected_trials_per_group);
        end

        if strcmp(group_name, 'SME')
            subject_stats.SME_removed = numel(bad_trials_relative);
        elseif strcmp(group_name, 'Test_Intact')
            subject_stats.Test_Intact_removed = numel(bad_trials_relative);
        elseif strcmp(group_name, 'Test_Recombined')
            subject_stats.Test_Recombined_removed = numel(bad_trials_relative);
        end

        bad_trials_absolute = group_trial_indices(bad_trials_relative);
        all_bad_trials_absolute = [all_bad_trials_absolute(:); bad_trials_absolute(:)]; %#ok<AGROW>
    end

    for c = 1:numel(condition_names)
        cond_name = condition_names{c};
        results.(['cond_orig_' cond_name])(file_idx) = subject_condition_original.(cond_name);
    end

    if exclude_subject
        results.excluded(file_idx) = true;
        results.exclusion_reason{file_idx} = sprintf('Subject %d: %s', s, exclusion_reason);
        fprintf('EXCLUDED subject %d: %s\n', s, exclusion_reason);
        for c = 1:numel(condition_names)
            cond_name = condition_names{c};
            results.(['cond_final_' cond_name])(file_idx) = 0;
            results.(['cond_removed_' cond_name])(file_idx) = subject_condition_original.(cond_name);
        end
        continue;
    end

    % STEP 5: reject bad trials and save
    final_bad_trials = unique(all_bad_trials_absolute);
    fprintf('Total unique bad trials removed: %d\n', numel(final_bad_trials));

    if ~isempty(final_bad_trials)
        EEG = pop_select(EEG, 'notrial', final_bad_trials);
    end

    [EEG, cond_counts_final] = assign_trial_type_from_event_label(EEG, epoch_triggers);
    subject_condition_final = condition_struct_from_map(cond_counts_final, condition_names);
    for c = 1:numel(condition_names)
        cond_name = condition_names{c};
        results.(['cond_final_' cond_name])(file_idx) = subject_condition_final.(cond_name);
        results.(['cond_removed_' cond_name])(file_idx) = ...
            subject_condition_original.(cond_name) - subject_condition_final.(cond_name);
    end

    subject_stats.SME_final = max(subject_stats.SME_original - subject_stats.SME_removed, 0);
    subject_stats.Test_Intact_final = max(subject_stats.Test_Intact_original - subject_stats.Test_Intact_removed, 0);
    subject_stats.Test_Recombined_final = max(subject_stats.Test_Recombined_original - subject_stats.Test_Recombined_removed, 0);

    results.final_trials(file_idx) = EEG.trials;
    results.removed_trials(file_idx) = numel(final_bad_trials);
    results.condition1_removals(file_idx) = subject_stats.condition1_total;
    results.condition2_removals(file_idx) = subject_stats.condition2_total;
    results.SME_original(file_idx) = subject_stats.SME_original;
    results.SME_removals(file_idx) = subject_stats.SME_removed;
    results.SME_final(file_idx) = subject_stats.SME_final;
    results.Test_Intact_original(file_idx) = subject_stats.Test_Intact_original;
    results.Test_Intact_removals(file_idx) = subject_stats.Test_Intact_removed;
    results.Test_Intact_final(file_idx) = subject_stats.Test_Intact_final;
    results.Test_Recombined_original(file_idx) = subject_stats.Test_Recombined_original;
    results.Test_Recombined_removals(file_idx) = subject_stats.Test_Recombined_removed;
    results.Test_Recombined_final(file_idx) = subject_stats.Test_Recombined_final;

    EEG = eeg_checkset(EEG, 'eventconsistency');
    EEG.setname = sprintf('Subject %d - All Conditions Clean', s);
    EEG.filename = output_filename;
    pop_saveset(EEG, 'filename', output_filename, 'filepath', epoch_folder, 'savemode', 'twofiles');

    results.processed(file_idx) = true;
    fprintf('Saved: %s\n', output_filepath);
end

stats = aggregate_epoch_stats(stats, results, condition_names);
processing_summary = build_processing_summary(stats, condition_names);

fprintf('\n==============================================\n');
fprintf('  EPOCHING COMPLETE\n');
fprintf('==============================================\n');
fprintf('Processed subjects: %d\n', stats.processed_subjects);
fprintf('Skipped subjects:   %d\n', sum(results.skipped));
fprintf('Excluded subjects:  %d\n', stats.excluded_subjects);
if stats.excluded_subjects > 0
    fprintf('Exclusions:\n');
    for i = 1:numel(stats.exclusion_reasons)
        fprintf('  - %s\n', stats.exclusion_reasons{i});
    end
end

try
    write_summary_text_file(summary_txt_file, processing_summary, stats);
    fprintf('Processing summary: %s\n', summary_txt_file);
catch ME
    fprintf('WARNING: Could not write summary text file (%s)\n', ME.message);
end

fprintf('Output folder: %s\n', epoch_folder);
fprintf('==============================================\n\n');

end

function stats = initialize_epoch_stats(total_subjects, voltage_diff_threshold, voltage_abs_threshold, max_rejected_trials_per_group, condition_names)
stats = struct();
stats.total_subjects = total_subjects;
stats.processed_subjects = 0;
stats.excluded_subjects = 0;
stats.excluded_subject_ids = {};
stats.voltage_diff_threshold = voltage_diff_threshold;
stats.voltage_abs_threshold = voltage_abs_threshold;
stats.max_rejected_trials_per_group = max_rejected_trials_per_group;
stats.total_original_trials = 0;
stats.total_final_trials = 0;
stats.total_removed_trials = 0;
stats.SME_removals = [];
stats.SME_original_trials = [];
stats.SME_final_trials = [];
stats.Test_Intact_removals = [];
stats.Test_Intact_original_trials = [];
stats.Test_Intact_final_trials = [];
stats.Test_Recombined_removals = [];
stats.Test_Recombined_original_trials = [];
stats.Test_Recombined_final_trials = [];
stats.condition1_removals = [];
stats.condition2_removals = [];
stats.total_removals_per_subject = [];
stats.subject_ids = [];
stats.subject_original_trials = [];
stats.subject_final_trials = [];
stats.subject_removed_trials = [];
stats.condition_original_trials = struct();
stats.condition_final_trials = struct();
stats.condition_removed_trials = struct();
for c = 1:numel(condition_names)
    cond_name = condition_names{c};
    stats.condition_original_trials.(cond_name) = [];
    stats.condition_final_trials.(cond_name) = [];
    stats.condition_removed_trials.(cond_name) = [];
end
stats.exclusion_reasons = {};
end

function results = initialize_epoch_results(total_subjects, condition_names)
results = struct();
results.processed = false(total_subjects, 1);
results.skipped = false(total_subjects, 1);
results.excluded = false(total_subjects, 1);
results.exclusion_reason = cell(total_subjects, 1);
results.subject_id = cell(total_subjects, 1);
results.original_trials = zeros(total_subjects, 1);
results.final_trials = zeros(total_subjects, 1);
results.removed_trials = zeros(total_subjects, 1);
results.condition1_removals = zeros(total_subjects, 1);
results.condition2_removals = zeros(total_subjects, 1);
results.SME_removals = zeros(total_subjects, 1);
results.SME_original = zeros(total_subjects, 1);
results.SME_final = zeros(total_subjects, 1);
results.Test_Intact_removals = zeros(total_subjects, 1);
results.Test_Intact_original = zeros(total_subjects, 1);
results.Test_Intact_final = zeros(total_subjects, 1);
results.Test_Recombined_removals = zeros(total_subjects, 1);
results.Test_Recombined_original = zeros(total_subjects, 1);
results.Test_Recombined_final = zeros(total_subjects, 1);
for c = 1:numel(condition_names)
    cond_name = condition_names{c};
    results.(['cond_orig_' cond_name]) = zeros(total_subjects, 1);
    results.(['cond_final_' cond_name]) = zeros(total_subjects, 1);
    results.(['cond_removed_' cond_name]) = zeros(total_subjects, 1);
end
end

function out = condition_struct_from_map(cond_counts, condition_names)
out = struct();
for c = 1:numel(condition_names)
    out.(condition_names{c}) = 0;
end
if isempty(cond_counts)
    return;
end
map_keys = cond_counts.keys;
for i = 1:numel(map_keys)
    key = map_keys{i};
    if isfield(out, key)
        out.(key) = cond_counts(key);
    end
end
end

function stats = aggregate_epoch_stats(stats, results, condition_names)
stats.processed_subjects = sum(results.processed);
stats.excluded_subjects = sum(results.excluded);

excluded_indices = find(results.excluded);
stats.excluded_subject_ids = results.subject_id(excluded_indices);
stats.exclusion_reasons = results.exclusion_reason(excluded_indices);

processed_not_excluded = results.processed & ~results.excluded & ~results.skipped;

stats.total_original_trials = sum(results.original_trials(processed_not_excluded));
stats.total_final_trials = sum(results.final_trials(processed_not_excluded));
stats.total_removed_trials = sum(results.removed_trials(processed_not_excluded));

stats.condition1_removals = results.condition1_removals(processed_not_excluded);
stats.condition2_removals = results.condition2_removals(processed_not_excluded);
stats.total_removals_per_subject = results.removed_trials(processed_not_excluded);

stats.SME_removals = results.SME_removals(processed_not_excluded);
stats.SME_original_trials = results.SME_original(processed_not_excluded);
stats.SME_final_trials = results.SME_final(processed_not_excluded);
stats.Test_Intact_removals = results.Test_Intact_removals(processed_not_excluded);
stats.Test_Intact_original_trials = results.Test_Intact_original(processed_not_excluded);
stats.Test_Intact_final_trials = results.Test_Intact_final(processed_not_excluded);
stats.Test_Recombined_removals = results.Test_Recombined_removals(processed_not_excluded);
stats.Test_Recombined_original_trials = results.Test_Recombined_original(processed_not_excluded);
stats.Test_Recombined_final_trials = results.Test_Recombined_final(processed_not_excluded);

for c = 1:numel(condition_names)
    cond_name = condition_names{c};
    stats.condition_original_trials.(cond_name) = results.(['cond_orig_' cond_name])(processed_not_excluded);
    stats.condition_final_trials.(cond_name) = results.(['cond_final_' cond_name])(processed_not_excluded);
    stats.condition_removed_trials.(cond_name) = results.(['cond_removed_' cond_name])(processed_not_excluded);
end

valid_subject_ids = results.subject_id(processed_not_excluded);
if isempty(valid_subject_ids)
    stats.subject_ids = [];
else
    stats.subject_ids = cellfun(@str2double, valid_subject_ids);
end
stats.subject_original_trials = results.original_trials(processed_not_excluded);
stats.subject_final_trials = results.final_trials(processed_not_excluded);
stats.subject_removed_trials = results.removed_trials(processed_not_excluded);
end

function processing_summary = build_processing_summary(stats, condition_names)
processing_summary = struct();

processing_summary.metadata = struct();
processing_summary.metadata.processing_date = datestr(now);
processing_summary.metadata.total_subjects_selected = stats.total_subjects;
processing_summary.metadata.successfully_processed = stats.processed_subjects;
processing_summary.metadata.excluded_subjects = stats.excluded_subjects;
processing_summary.metadata.excluded_subject_ids = stats.excluded_subject_ids;
processing_summary.metadata.exclusion_reasons = stats.exclusion_reasons;
processing_summary.metadata.processing_success_rate = safe_percent(stats.processed_subjects, stats.total_subjects);

processing_summary.parameters = struct();
processing_summary.parameters.voltage_diff_threshold = stats.voltage_diff_threshold;
processing_summary.parameters.voltage_abs_threshold = stats.voltage_abs_threshold;
processing_summary.parameters.max_rejected_trials_per_group = stats.max_rejected_trials_per_group;
processing_summary.parameters.epoch_window = [-0.1, 1.5];
processing_summary.parameters.baseline_window = [-100, 0];

processing_summary.overall = struct();
processing_summary.overall.total_original_trials = stats.total_original_trials;
processing_summary.overall.total_final_trials = stats.total_final_trials;
processing_summary.overall.total_removed_trials = stats.total_removed_trials;
processing_summary.overall.overall_removal_rate = safe_percent(stats.total_removed_trials, stats.total_original_trials);
processing_summary.overall.data_retention_rate = safe_percent(stats.total_final_trials, stats.total_original_trials);

processing_summary.condition_based = struct();
processing_summary.condition_based.condition1_removals = stats.condition1_removals;
processing_summary.condition_based.condition2_removals = stats.condition2_removals;
processing_summary.condition_based.condition1_mean = safe_mean(stats.condition1_removals);
processing_summary.condition_based.condition1_std = safe_std(stats.condition1_removals);
processing_summary.condition_based.condition2_mean = safe_mean(stats.condition2_removals);
processing_summary.condition_based.condition2_std = safe_std(stats.condition2_removals);
processing_summary.condition_based.overlap_trials = ...
    sum(stats.condition1_removals) + sum(stats.condition2_removals) - stats.total_removed_trials;

processing_summary.group_based = struct();
processing_summary.group_based.SME = build_group_summary(stats.SME_original_trials, stats.SME_removals, stats.SME_final_trials);
processing_summary.group_based.Test_Intact = build_group_summary(stats.Test_Intact_original_trials, stats.Test_Intact_removals, stats.Test_Intact_final_trials);
processing_summary.group_based.Test_Recombined = build_group_summary(stats.Test_Recombined_original_trials, stats.Test_Recombined_removals, stats.Test_Recombined_final_trials);

processing_summary.experimental_conditions = struct();
for c = 1:numel(condition_names)
    cond_name = condition_names{c};
    original_trials = stats.condition_original_trials.(cond_name);
    removed_trials = stats.condition_removed_trials.(cond_name);
    final_trials = stats.condition_final_trials.(cond_name);

    processing_summary.experimental_conditions.(cond_name) = struct();
    processing_summary.experimental_conditions.(cond_name).original_trials = original_trials;
    processing_summary.experimental_conditions.(cond_name).removed_trials = removed_trials;
    processing_summary.experimental_conditions.(cond_name).final_trials = final_trials;
    processing_summary.experimental_conditions.(cond_name).mean_original = safe_mean(original_trials);
    processing_summary.experimental_conditions.(cond_name).mean_removed = safe_mean(removed_trials);
    processing_summary.experimental_conditions.(cond_name).mean_final = safe_mean(final_trials);
    processing_summary.experimental_conditions.(cond_name).std_original = safe_std(original_trials);
    processing_summary.experimental_conditions.(cond_name).std_removed = safe_std(removed_trials);
    processing_summary.experimental_conditions.(cond_name).std_final = safe_std(final_trials);
    processing_summary.experimental_conditions.(cond_name).removal_rate = safe_percent(sum(removed_trials), sum(original_trials));
end

processing_summary.subject_based = struct();
processing_summary.subject_based.subject_ids = stats.subject_ids;
processing_summary.subject_based.original_trials = stats.subject_original_trials;
processing_summary.subject_based.removed_trials = stats.subject_removed_trials;
processing_summary.subject_based.final_trials = stats.subject_final_trials;
processing_summary.subject_based.mean_original = safe_mean(stats.subject_original_trials);
processing_summary.subject_based.mean_removed = safe_mean(stats.subject_removed_trials);
processing_summary.subject_based.mean_final = safe_mean(stats.subject_final_trials);
processing_summary.subject_based.std_original = safe_std(stats.subject_original_trials);
processing_summary.subject_based.std_removed = safe_std(stats.subject_removed_trials);
processing_summary.subject_based.std_final = safe_std(stats.subject_final_trials);
processing_summary.subject_based.removal_rate = safe_percent(sum(stats.subject_removed_trials), sum(stats.subject_original_trials));
processing_summary.subject_based.range_original = [safe_min(stats.subject_original_trials), safe_max(stats.subject_original_trials)];
processing_summary.subject_based.range_removed = [safe_min(stats.subject_removed_trials), safe_max(stats.subject_removed_trials)];
processing_summary.subject_based.range_final = [safe_min(stats.subject_final_trials), safe_max(stats.subject_final_trials)];

processing_summary.quality_assessment = struct();
processing_summary.quality_assessment.retention_rate = safe_percent(stats.total_final_trials, stats.total_original_trials);
if processing_summary.quality_assessment.retention_rate >= 85
    processing_summary.quality_assessment.quality_rating = 'EXCELLENT';
elseif processing_summary.quality_assessment.retention_rate >= 75
    processing_summary.quality_assessment.quality_rating = 'GOOD';
elseif processing_summary.quality_assessment.retention_rate >= 65
    processing_summary.quality_assessment.quality_rating = 'ACCEPTABLE';
else
    processing_summary.quality_assessment.quality_rating = 'CONCERNING';
end
processing_summary.quality_assessment.subjects_approaching_threshold = sum(stats.total_removals_per_subject > 40);
end

function group_summary = build_group_summary(original_trials, removed_trials, final_trials)
group_summary = struct();
group_summary.original_trials = original_trials;
group_summary.removed_trials = removed_trials;
group_summary.final_trials = final_trials;
group_summary.mean_original = safe_mean(original_trials);
group_summary.mean_removed = safe_mean(removed_trials);
group_summary.mean_final = safe_mean(final_trials);
group_summary.std_original = safe_std(original_trials);
group_summary.std_removed = safe_std(removed_trials);
group_summary.std_final = safe_std(final_trials);
group_summary.removal_rate = safe_percent(sum(removed_trials), sum(original_trials));
end

function write_summary_text_file(summary_txt_file, processing_summary, stats)
[fid, msg] = fopen(summary_txt_file, 'w');
if fid == -1
    error('Unable to open summary file: %s', msg);
end
cleanup_fid = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'Processing Summary Export\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));
fprintf(fid, '[processing_summary]\n');
write_struct_lines(fid, 'processing_summary', processing_summary);
fprintf(fid, '\n[stats]\n');
write_struct_lines(fid, 'stats', stats);
end

function write_struct_lines(fid, prefix, value)
if isstruct(value)
    fields = fieldnames(value);
    for i = 1:numel(fields)
        field_name = fields{i};
        write_struct_lines(fid, [prefix '.' field_name], value.(field_name));
    end
    return;
end
fprintf(fid, '%s = %s\n', prefix, format_value_for_text(value));
end

function out = format_value_for_text(value)
if isnumeric(value) || islogical(value)
    out = mat2str(value);
elseif ischar(value)
    out = value;
elseif isstring(value)
    if isscalar(value)
        out = char(value);
    else
        out = strjoin(cellstr(value), ', ');
    end
elseif iscell(value)
    if isempty(value)
        out = '{}';
    else
        parts = cell(1, numel(value));
        for i = 1:numel(value)
            parts{i} = format_value_for_text(value{i});
        end
        out = ['{' strjoin(parts, ', ') '}'];
    end
else
    out = '<unsupported>';
end
end

function value = safe_mean(x)
if isempty(x)
    value = NaN;
else
    value = mean(x);
end
end

function value = safe_std(x)
if isempty(x)
    value = NaN;
else
    value = std(x);
end
end

function value = safe_min(x)
if isempty(x)
    value = NaN;
else
    value = min(x);
end
end

function value = safe_max(x)
if isempty(x)
    value = NaN;
else
    value = max(x);
end
end

function value = safe_percent(numerator, denominator)
if denominator <= 0
    value = 0;
else
    value = (numerator / denominator) * 100;
end
end

function out = parse_optional_subject_filter(value)
if isnumeric(value) && isscalar(value) && isfinite(value) && mod(value, 1) == 0
    out = double(value);
    return;
end

if (ischar(value) && ~isempty(strtrim(value))) || (isstring(value) && isscalar(value))
    n = str2double(strtrim(char(value)));
    if isfinite(n) && mod(n, 1) == 0
        out = double(n);
        return;
    end
end

error('subject_filter must be a finite integer subject ID.');
end

function triggers = normalize_trigger_list(value)
if ischar(value) || (isstring(value) && isscalar(value))
    parts = strsplit(char(value), ',');
    triggers = cellfun(@strtrim, parts, 'UniformOutput', false);
elseif isnumeric(value)
    triggers = arrayfun(@(x) num2str(x), value(:)', 'UniformOutput', false);
elseif iscell(value)
    triggers = cell(1, numel(value));
    for i = 1:numel(value)
        if isnumeric(value{i})
            triggers{i} = num2str(value{i});
        else
            triggers{i} = strtrim(char(value{i}));
        end
    end
else
    error('Unsupported trigger list input type.');
end

triggers = triggers(~cellfun(@isempty, triggers));
if isempty(triggers)
    error('No valid trigger codes provided.');
end

triggers = unique(triggers, 'stable');
end

function group_def = parse_group_spec(group_spec)
if ischar(group_spec) || (isstring(group_spec) && isscalar(group_spec))
    spec = char(group_spec);
    entries = strsplit(spec, ';');
    group_def = cell(0, 2);

    for i = 1:numel(entries)
        entry = strtrim(entries{i});
        if isempty(entry)
            continue;
        end
        parts = strsplit(entry, ':');
        if numel(parts) ~= 2
            error('Invalid group spec entry: %s. Expected format "Group:11,21"', entry);
        end

        group_name = strtrim(parts{1});
        group_triggers = normalize_trigger_list(parts{2});
        group_def(end+1, :) = {group_name, group_triggers}; %#ok<AGROW>
    end

elseif iscell(group_spec)
    if size(group_spec, 2) ~= 2
        error('Cell group_spec must be Nx2: {group_name, trigger_list}.');
    end
    group_def = cell(size(group_spec, 1), 2);
    for i = 1:size(group_spec, 1)
        group_def{i, 1} = char(group_spec{i, 1});
        group_def{i, 2} = normalize_trigger_list(group_spec{i, 2});
    end
else
    error('Unsupported group_spec input type.');
end

if isempty(group_def)
    error('No artifact groups were defined.');
end
end

function [EEG, condition_counts] = assign_trial_type_from_event_label(EEG, epoch_triggers)
[EEG.event.trial_type] = deal('');
condition_counts = containers.Map('KeyType', 'char', 'ValueType', 'double');

num_events = length(EEG.event);
event_epochs = zeros(1, num_events);
event_latencies = zeros(1, num_events);
event_types = cell(1, num_events);
event_labels = cell(1, num_events);

for ev_idx = 1:num_events
    event_epochs(ev_idx) = EEG.event(ev_idx).epoch;
    event_latencies(ev_idx) = EEG.event(ev_idx).latency;
    event_types{ev_idx} = event_type_to_string(EEG.event(ev_idx).type);
    event_labels{ev_idx} = normalize_event_label_value(EEG.event(ev_idx).event_label);
end

for epoch_num = 1:EEG.trials
    epoch_event_indices = find(event_epochs == epoch_num);
    if isempty(epoch_event_indices)
        continue;
    end

    selected_idx = [];
    earliest_latency = inf;
    for idx = epoch_event_indices
        if ismember(event_types{idx}, epoch_triggers)
            if event_latencies(idx) < earliest_latency
                earliest_latency = event_latencies(idx);
                selected_idx = idx;
            end
        end
    end

    if isempty(selected_idx)
        continue;
    end

    trial_label = event_labels{selected_idx};
    if isempty(trial_label)
        for idx = epoch_event_indices
            if ~isempty(event_labels{idx})
                trial_label = event_labels{idx};
                break;
            end
        end
    end

    if isempty(trial_label)
        continue;
    end

    for idx = epoch_event_indices
        EEG.event(idx).trial_type = trial_label;
    end

    if isKey(condition_counts, trial_label)
        condition_counts(trial_label) = condition_counts(trial_label) + 1;
    else
        condition_counts(trial_label) = 1;
    end
end
end

function trial_indices_by_group = assign_epochs_to_groups(EEG, epoch_triggers, group_def)
num_groups = size(group_def, 1);
trial_indices_by_group = cell(num_groups, 1);
for g = 1:num_groups
    trial_indices_by_group{g} = [];
end

num_events = length(EEG.event);
event_epochs = zeros(1, num_events);
event_latencies = zeros(1, num_events);
event_types = cell(1, num_events);

for ev_idx = 1:num_events
    event_epochs(ev_idx) = EEG.event(ev_idx).epoch;
    event_latencies(ev_idx) = EEG.event(ev_idx).latency;
    event_types{ev_idx} = event_type_to_string(EEG.event(ev_idx).type);
end

epoch_trigger = cell(1, EEG.trials);
for epoch_num = 1:EEG.trials
    epoch_event_indices = find(event_epochs == epoch_num);
    if isempty(epoch_event_indices)
        continue;
    end

    selected_trigger = '';
    earliest_latency = inf;
    for idx = epoch_event_indices
        t = event_types{idx};
        if ismember(t, epoch_triggers) && event_latencies(idx) < earliest_latency
            earliest_latency = event_latencies(idx);
            selected_trigger = t;
        end
    end
    epoch_trigger{epoch_num} = selected_trigger;
end

for epoch_num = 1:EEG.trials
    t = epoch_trigger{epoch_num};
    if isempty(t)
        continue;
    end

    for g = 1:num_groups
        g_triggers = group_def{g, 2};
        if ismember(t, g_triggers)
            trial_indices_by_group{g}(end+1) = epoch_num; %#ok<AGROW>
            break;
        end
    end
end
end

function out = event_type_to_string(value)
if isnumeric(value)
    out = num2str(value);
elseif isstring(value)
    out = strtrim(char(value));
else
    out = strtrim(char(value));
end
end

function out = normalize_event_label_value(value)
if iscell(value)
    if isempty(value)
        out = '';
    else
        out = strtrim(char(value{1}));
    end
elseif isstring(value)
    out = strtrim(char(value));
elseif ischar(value)
    out = strtrim(value);
else
    out = '';
end
end

function [clean_EEG, bad_trials, condition1_bad, condition2_bad] = get_bad_trials_v2(EEG, voltage_diff_threshold, voltage_abs_threshold, all_channels)
bad_trials_condition1 = [];
bad_trials_condition2 = [];
for epoch_idx = 1:size(EEG.data, 3)
    for ch = all_channels
        voltage_data = squeeze(EEG.data(ch, :, epoch_idx));
        if max(abs(diff(voltage_data))) > voltage_diff_threshold
            bad_trials_condition1 = [bad_trials_condition1, epoch_idx]; %#ok<AGROW>
        end
        if max(abs(voltage_data)) > voltage_abs_threshold
            bad_trials_condition2 = [bad_trials_condition2, epoch_idx]; %#ok<AGROW>
        end
    end
end
condition1_bad = unique(bad_trials_condition1);
condition2_bad = unique(bad_trials_condition2);
bad_trials = unique([condition1_bad, condition2_bad]);
clean_EEG = EEG;
end
