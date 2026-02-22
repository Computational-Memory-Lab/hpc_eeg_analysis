function hpc_interpol_to_epoch(input_folder, voltage_diff_threshold, voltage_abs_threshold, epoch_triggers, group_spec)
% HPC_INTERPOL_TO_EPOCH - Epoch EEG data and apply artifact rejection
%
% Usage:
%   hpc_interpol_to_epoch(input_folder)
%   hpc_interpol_to_epoch(input_folder, voltage_diff_threshold, voltage_abs_threshold)
%   hpc_interpol_to_epoch(input_folder, voltage_diff_threshold, voltage_abs_threshold, epoch_triggers, group_spec)
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

epoch_triggers = normalize_trigger_list(epoch_triggers);
group_def = parse_group_spec(group_spec);

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
        file_info(end+1).filename = fname; %#ok<AGROW>
        file_info(end).filepath = filtered_files(i).folder;
        file_info(end).subject_id = str2double(tokens{1}{1});
    else
        fprintf('WARNING: Could not extract subject ID from "%s", skipping\n', fname);
    end
end
if isempty(file_info)
    error('No .set files with parseable subject IDs found in: %s', input_folder);
end

fprintf('Files to process: %d\n', length(file_info));

% Create output directory
epoch_folder = fullfile(input_folder, 'epoch');
if ~exist(epoch_folder, 'dir')
    mkdir(epoch_folder);
end

% Processing settings
max_rejected_trials_per_group = 56;
processed_count = 0;
skipped_count = 0;
excluded_subjects = {};

for file_idx = 1:length(file_info)
    current_file = file_info(file_idx);
    s = current_file.subject_id;

    output_filename = sprintf('%d_epoch.set', s);
    output_filepath = fullfile(epoch_folder, output_filename);

    fprintf('\n----------------------------------------------\n');
    fprintf('Subject %d (%d/%d)\n', s, file_idx, length(file_info));

    if exist(output_filepath, 'file')
        fprintf('Output already exists, skipping: %s\n', output_filename);
        skipped_count = skipped_count + 1;
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
    fprintf('Epoched %d trials\n', EEG.trials);

    % STEP 2: set trial_type directly from event_label
    [EEG, cond_counts] = assign_trial_type_from_event_label(EEG, epoch_triggers);
    condition_names = sort(cond_counts.keys);
    fprintf('Condition counts (pre-rejection):\n');
    for c = 1:numel(condition_names)
        fprintf('  %-30s %d\n', condition_names{c}, cond_counts(condition_names{c}));
    end

    % STEP 3: assign epochs to artifact groups
    trial_indices_by_group = assign_epochs_to_groups(EEG, epoch_triggers, group_def);

    % STEP 4: detect bad trials per group
    all_bad_trials_absolute = [];
    exclude_subject = false;
    exclusion_reason = '';
    all_channels = 1:EEG.nbchan;

    for g = 1:size(group_def, 1)
        group_name = group_def{g, 1};
        group_trial_indices = trial_indices_by_group{g};

        fprintf('Group %s: %d trials\n', group_name, numel(group_trial_indices));
        if isempty(group_trial_indices)
            continue;
        end

        EEG_group = pop_select(EEG, 'trial', group_trial_indices);
        [~, bad_trials_relative, cond1_bad, cond2_bad] = get_bad_trials_v2(EEG_group, ...
            voltage_diff_threshold, voltage_abs_threshold, all_channels);
        fprintf('  Bad trials in %s: %d (diff=%d, abs=%d)\n', ...
            group_name, numel(bad_trials_relative), numel(cond1_bad), numel(cond2_bad));

        if numel(bad_trials_relative) > max_rejected_trials_per_group
            exclude_subject = true;
            exclusion_reason = sprintf('Too many trials rejected in %s (%d > %d)', ...
                group_name, numel(bad_trials_relative), max_rejected_trials_per_group);
        end

        bad_trials_absolute = group_trial_indices(bad_trials_relative);
        all_bad_trials_absolute = [all_bad_trials_absolute(:); bad_trials_absolute(:)]; %#ok<AGROW>
    end

    if exclude_subject
        excluded_subjects{end+1} = sprintf('%d (%s)', s, exclusion_reason); %#ok<AGROW>
        fprintf('EXCLUDED subject %d: %s\n', s, exclusion_reason);
        continue;
    end

    % STEP 5: reject bad trials and save
    final_bad_trials = unique(all_bad_trials_absolute);
    fprintf('Total unique bad trials removed: %d\n', numel(final_bad_trials));

    if ~isempty(final_bad_trials)
        EEG = pop_select(EEG, 'notrial', final_bad_trials);
    end

    EEG = eeg_checkset(EEG, 'eventconsistency');
    EEG.setname = sprintf('Subject %d - All Conditions Clean', s);
    EEG.filename = output_filename;
    pop_saveset(EEG, 'filename', output_filename, 'filepath', epoch_folder, 'savemode', 'twofiles');

    processed_count = processed_count + 1;
    fprintf('Saved: %s\n', output_filepath);
end

fprintf('\n==============================================\n');
fprintf('  EPOCHING COMPLETE\n');
fprintf('==============================================\n');
fprintf('Processed subjects: %d\n', processed_count);
fprintf('Skipped subjects:   %d\n', skipped_count);
fprintf('Excluded subjects:  %d\n', numel(excluded_subjects));
if ~isempty(excluded_subjects)
    fprintf('Exclusions:\n');
    for i = 1:numel(excluded_subjects)
        fprintf('  - %s\n', excluded_subjects{i});
    end
end
fprintf('Output folder: %s\n', epoch_folder);
fprintf('==============================================\n\n');

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
