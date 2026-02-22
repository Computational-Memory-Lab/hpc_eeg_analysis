function hpc_interpol_to_epoch(input_folder, voltage_diff_threshold, voltage_abs_threshold)
% HPC_INTERPOL_TO_EPOCH - Epoch EEG data and apply artifact rejection
%
% Usage:
%   hpc_interpol_to_epoch(input_folder)
%   hpc_interpol_to_epoch(input_folder, voltage_diff_threshold, voltage_abs_threshold)
%
% Inputs:
%   input_folder           - Path to the interpol folder containing .set files
%   voltage_diff_threshold - (Optional) Max voltage difference threshold in uV (default: 20)
%   voltage_abs_threshold  - (Optional) Max absolute voltage threshold in uV (default: 1000)
%
% Outputs:
%   - <input_folder>/epoch/<ID>_epoch.set  Epoched and artifact-rejected EEG datasets
%
% Workflow:
%   1. Epoch all trials together (triggers: 11, 21, 22), window [-100, 1500] ms
%   2. Add custom 'trial_type' field based on trigger and accuracy
%   3. Identify trial indices for each artifact rejection group (SME, Test_Intact, Test_Recombined)
%   4. Detect bad trials per group using voltage thresholds
%   5. Remove all bad trials in a single step
%   6. Save the clean dataset

% Set defaults for optional arguments
if nargin < 2 || isempty(voltage_diff_threshold)
    voltage_diff_threshold = 20;
end
if nargin < 3 || isempty(voltage_abs_threshold)
    voltage_abs_threshold = 1000;
end

% ==================== INITIALIZE EEGLAB (HEADLESS MODE) ====================
eeglab_path = '/home/devon7y/scratch/devon7y/eeglab2022.1';
addpath(eeglab_path);
eeglab nogui;
fprintf('EEGLAB initialized in headless mode\n');
fprintf('Voltage difference threshold: %.1f uV\n', voltage_diff_threshold);
fprintf('Absolute voltage threshold:   %.1f uV\n', voltage_abs_threshold);

% ==================== FIND INPUT FILES ====================
set_files = dir(fullfile(input_folder, '*.set'));
if isempty(set_files)
    error('No .set files found in: %s', input_folder);
end

% Filter out checkpoint files and already-epoched files
filtered_files = struct('name', {}, 'folder', {}, 'date', {}, 'bytes', {}, 'isdir', {}, 'datenum', {});
for i = 1:length(set_files)
    fname = set_files(i).name;
    if ~startsWith(fname, 'before_step9') && ~contains(fname, '_epoch')
        filtered_files(end+1) = set_files(i);
    end
end

if isempty(filtered_files)
    error('No valid input .set files found in: %s (after filtering checkpoints/epoch files)', input_folder);
end

fprintf('Found %d .set files to process\n', length(filtered_files));

% Build file_info struct
file_info = struct('filename', {}, 'filepath', {}, 'subject_id', {});
for i = 1:length(filtered_files)
    fname = filtered_files(i).name;
    [~, basename, ~] = fileparts(fname);
    tokens = regexp(basename, '(\d+)', 'tokens');
    if ~isempty(tokens)
        subject_id = str2double(tokens{1}{1});
        file_info(end+1).filename = fname;
        file_info(end).filepath = filtered_files(i).folder;
        file_info(end).subject_id = subject_id;
    else
        fprintf('WARNING: Could not extract subject ID from "%s", skipping\n', fname);
    end
end

if isempty(file_info)
    error('No .set files with parseable subject IDs found in: %s', input_folder);
end

fprintf('Total files to process: %d\n', length(file_info));

% ==================== CREATE OUTPUT DIRECTORY ====================
epoch_folder = fullfile(input_folder, 'epoch');
if ~exist(epoch_folder, 'dir'), mkdir(epoch_folder); end
fprintf('Output will be saved to: %s\n', epoch_folder);

% ==================== INITIALIZE STATISTICS TRACKING ====================
condition_names = {'Study_hits', 'Study_misses', 'Test_hits', 'Test_misses', 'Correct_rejections', 'False_alarms'};

stats = struct();
stats.total_subjects = length(file_info);
stats.processed_subjects = 0;
stats.excluded_subjects = 0;
stats.excluded_subject_ids = {};
stats.voltage_diff_threshold = voltage_diff_threshold;
stats.voltage_abs_threshold = voltage_abs_threshold;
stats.max_rejected_trials_per_group = 56;
stats.total_original_trials = 0;
stats.total_final_trials = 0;
stats.total_removed_trials = 0;
stats.SME_removals = [];
stats.Test_Intact_removals = [];
stats.Test_Recombined_removals = [];
stats.condition1_removals = [];
stats.condition2_removals = [];
stats.total_removals_per_subject = [];
stats.SME_original_trials = [];
stats.SME_final_trials = [];
stats.Test_Intact_original_trials = [];
stats.Test_Intact_final_trials = [];
stats.Test_Recombined_original_trials = [];
stats.Test_Recombined_final_trials = [];
stats.subject_ids = [];
stats.subject_original_trials = [];
stats.subject_final_trials = [];
stats.subject_removed_trials = [];
stats.condition_original_trials = struct();
stats.condition_final_trials = struct();
stats.condition_removed_trials = struct();
for cond = condition_names
    stats.condition_original_trials.(cond{1}) = [];
    stats.condition_final_trials.(cond{1}) = [];
    stats.condition_removed_trials.(cond{1}) = [];
end
stats.exclusion_reasons = {};

% ==================== PRE-ALLOCATE RESULTS STORAGE ====================
num_files = length(file_info);
parallel_results = struct();
parallel_results.processed = false(num_files, 1);
parallel_results.skipped = false(num_files, 1);
parallel_results.excluded = false(num_files, 1);
parallel_results.exclusion_reason = cell(num_files, 1);
parallel_results.subject_id = cell(num_files, 1);
parallel_results.original_trials = zeros(num_files, 1);
parallel_results.final_trials = zeros(num_files, 1);
parallel_results.removed_trials = zeros(num_files, 1);
parallel_results.condition1_removals = zeros(num_files, 1);
parallel_results.condition2_removals = zeros(num_files, 1);
parallel_results.SME_removals = zeros(num_files, 1);
parallel_results.SME_original = zeros(num_files, 1);
parallel_results.SME_final = zeros(num_files, 1);
parallel_results.Test_Intact_removals = zeros(num_files, 1);
parallel_results.Test_Intact_original = zeros(num_files, 1);
parallel_results.Test_Intact_final = zeros(num_files, 1);
parallel_results.Test_Recombined_removals = zeros(num_files, 1);
parallel_results.Test_Recombined_original = zeros(num_files, 1);
parallel_results.Test_Recombined_final = zeros(num_files, 1);
for cond = condition_names
    parallel_results.(['cond_orig_' cond{1}]) = zeros(num_files, 1);
    parallel_results.(['cond_final_' cond{1}]) = zeros(num_files, 1);
    parallel_results.(['cond_removed_' cond{1}]) = zeros(num_files, 1);
end

% ==================== PROCESS EACH SUBJECT ====================
for file_idx = 1:length(file_info)
    current_file = file_info(file_idx);
    s = current_file.subject_id;
    input_filename = current_file.filename;
    input_filepath = current_file.filepath;

    output_filename = sprintf('%d_epoch.set', s);
    output_filepath = fullfile(epoch_folder, output_filename);

    parallel_results.subject_id{file_idx} = sprintf('%d', s);

    if exist(output_filepath, 'file')
        fprintf('\nSubject %d already processed (file exists: %s). Skipping...\n', s, output_filename);
        parallel_results.skipped(file_idx) = true;
        parallel_results.processed(file_idx) = true;
        continue;
    end

    EEG = pop_loadset('filename', input_filename, 'filepath', input_filepath);
    fprintf('\nProcessing subject %d...\n', s);

    if isfield(EEG.event, 'trial_type')
        EEG.event = rmfield(EEG.event, 'trial_type');
    end

    % --- STEP 1: Epoch all trials together ---
    all_triggers = {'11', '21', '22'};
    EEG = pop_epoch(EEG, all_triggers, [-0.1 1.5], 'newname', sprintf('Subject %d - all epoched', s), 'epochinfo', 'yes');
    EEG = pop_rmbase(EEG, [-100 0]);
    fprintf('  Epoched all %d trials together.\n', EEG.trials);
    parallel_results.original_trials(file_idx) = EEG.trials;

    % --- STEP 2: Add the custom 'trial_type' field ---
    fprintf('  Adding custom ''trial_type'' field to events...\n');

    condition_def = { ...
        'Study_hits',         '11', 1; ...
        'Study_misses',       '11', 0;
        'Test_hits',          '21', 1;
        'Test_misses',        '21', 0;
        'Correct_rejections', '22', 1;
        'False_alarms',       '22', 0};

    [EEG.event.trial_type] = deal('');

    subject_condition_counts = struct();
    subject_condition_counts.Study_hits = 0;
    subject_condition_counts.Study_misses = 0;
    subject_condition_counts.Test_hits = 0;
    subject_condition_counts.Test_misses = 0;
    subject_condition_counts.Correct_rejections = 0;
    subject_condition_counts.False_alarms = 0;

    epoch_condition_map = cell(1, EEG.trials);

    for i = 1:length(EEG.event)
        event_type_str = num2str(EEG.event(i).type);
        for c = 1:size(condition_def, 1)
            cond_name = condition_def{c, 1};
            cond_trigger = condition_def{c, 2};
            cond_acc = condition_def{c, 3};
            if isfield(EEG.event(i), 'Acc') && strcmp(event_type_str, cond_trigger) && EEG.event(i).Acc == cond_acc
                EEG.event(i).trial_type = cond_name;
                epoch_num = EEG.event(i).epoch;
                if isempty(epoch_condition_map{epoch_num})
                    epoch_condition_map{epoch_num} = cond_name;
                    subject_condition_counts.(cond_name) = subject_condition_counts.(cond_name) + 1;
                end
                break;
            end
        end
    end

    parallel_results.cond_orig_Study_hits(file_idx) = subject_condition_counts.Study_hits;
    parallel_results.cond_orig_Study_misses(file_idx) = subject_condition_counts.Study_misses;
    parallel_results.cond_orig_Test_hits(file_idx) = subject_condition_counts.Test_hits;
    parallel_results.cond_orig_Test_misses(file_idx) = subject_condition_counts.Test_misses;
    parallel_results.cond_orig_Correct_rejections(file_idx) = subject_condition_counts.Correct_rejections;
    parallel_results.cond_orig_False_alarms(file_idx) = subject_condition_counts.False_alarms;

    % --- STEP 3: Identify trial indices per artifact rejection group ---
    fprintf('  Identifying trials for artifact rejection groups...\n');
    trial_indices_by_group = assign_epochs_to_groups_fixed(EEG);
    fprintf('  Identified trials per group: SME(%d), Test_Intact(%d), Test_Recombined(%d)\n', ...
        length(trial_indices_by_group.SME), length(trial_indices_by_group.Test_Intact), length(trial_indices_by_group.Test_Recombined));

    total_assigned = length(trial_indices_by_group.SME) + length(trial_indices_by_group.Test_Intact) + length(trial_indices_by_group.Test_Recombined);
    if total_assigned ~= EEG.trials
        fprintf('  WARNING: Group assignment mismatch! Total epochs: %d, Assigned: %d\n', EEG.trials, total_assigned);
    else
        fprintf('  VALIDATION PASSED: Each epoch assigned to exactly one group (Total: %d)\n', total_assigned);
    end

    % --- STEP 4: Find bad trials per group ---
    all_bad_trials_absolute = [];
    subject_should_be_skipped = false;
    exclusion_reason = '';
    all_channels = 1:EEG.nbchan;

    subject_stats = struct();
    subject_stats.condition1_total = 0;
    subject_stats.condition2_total = 0;
    subject_stats.SME_original = 0;
    subject_stats.SME_removed = 0;
    subject_stats.SME_final = 0;
    subject_stats.Test_Intact_original = 0;
    subject_stats.Test_Intact_removed = 0;
    subject_stats.Test_Intact_final = 0;
    subject_stats.Test_Recombined_original = 0;
    subject_stats.Test_Recombined_removed = 0;
    subject_stats.Test_Recombined_final = 0;

    group_def = {'SME',             {'11'};
                 'Test_Intact',     {'21'};
                 'Test_Recombined', {'22'}};

    for g = 1:size(group_def, 1)
        group_name = group_def{g, 1};
        group_trial_indices = trial_indices_by_group.(group_name);

        if isempty(group_trial_indices)
            fprintf('    No trials found for group %s, skipping artifact detection.\n', group_name);
            continue;
        end

        EEG_group = pop_select(EEG, 'trial', group_trial_indices);
        group_original_trials = length(group_trial_indices);

        if strcmp(group_name, 'SME')
            subject_stats.SME_original = group_original_trials;
        elseif strcmp(group_name, 'Test_Intact')
            subject_stats.Test_Intact_original = group_original_trials;
        elseif strcmp(group_name, 'Test_Recombined')
            subject_stats.Test_Recombined_original = group_original_trials;
        end

        [~, bad_trials_relative, cond1_bad, cond2_bad] = get_bad_trials_v2(EEG_group, stats.voltage_diff_threshold, stats.voltage_abs_threshold, all_channels);
        fprintf('    Found %d bad trials in %s group.\n', length(bad_trials_relative), group_name);

        subject_stats.condition1_total = subject_stats.condition1_total + length(cond1_bad);
        subject_stats.condition2_total = subject_stats.condition2_total + length(cond2_bad);

        if length(bad_trials_relative) > stats.max_rejected_trials_per_group
            fprintf('    WARNING: Subject %d excluded. Too many trials rejected in %s group (%d > %d).\n', ...
                s, group_name, length(bad_trials_relative), stats.max_rejected_trials_per_group);
            subject_should_be_skipped = true;
            exclusion_reason = sprintf('Too many trials rejected in %s group (%d > %d)', ...
                group_name, length(bad_trials_relative), stats.max_rejected_trials_per_group);
        end

        bad_trials_absolute = group_trial_indices(bad_trials_relative);
        all_bad_trials_absolute = [all_bad_trials_absolute(:); bad_trials_absolute(:)];

        if strcmp(group_name, 'SME')
            subject_stats.SME_removed = length(bad_trials_relative);
        elseif strcmp(group_name, 'Test_Intact')
            subject_stats.Test_Intact_removed = length(bad_trials_relative);
        elseif strcmp(group_name, 'Test_Recombined')
            subject_stats.Test_Recombined_removed = length(bad_trials_relative);
        end

        parallel_results.([group_name '_removals'])(file_idx) = length(bad_trials_relative);
    end

    if subject_should_be_skipped
        parallel_results.excluded(file_idx) = true;
        parallel_results.exclusion_reason{file_idx} = sprintf('Subject %d: %s', s, exclusion_reason);
        parallel_results.cond_final_Study_hits(file_idx) = 0;
        parallel_results.cond_final_Study_misses(file_idx) = 0;
        parallel_results.cond_final_Test_hits(file_idx) = 0;
        parallel_results.cond_final_Test_misses(file_idx) = 0;
        parallel_results.cond_final_Correct_rejections(file_idx) = 0;
        parallel_results.cond_final_False_alarms(file_idx) = 0;
        parallel_results.cond_removed_Study_hits(file_idx) = subject_condition_counts.Study_hits;
        parallel_results.cond_removed_Study_misses(file_idx) = subject_condition_counts.Study_misses;
        parallel_results.cond_removed_Test_hits(file_idx) = subject_condition_counts.Test_hits;
        parallel_results.cond_removed_Test_misses(file_idx) = subject_condition_counts.Test_misses;
        parallel_results.cond_removed_Correct_rejections(file_idx) = subject_condition_counts.Correct_rejections;
        parallel_results.cond_removed_False_alarms(file_idx) = subject_condition_counts.False_alarms;
        continue;
    end

    % --- STEP 5: Reject all bad trials at once ---
    final_bad_trials = unique(all_bad_trials_absolute);
    fprintf('  Total unique bad trials to remove: %d\n', length(final_bad_trials));

    if ~isempty(final_bad_trials)
        EEG = pop_select(EEG, 'notrial', final_bad_trials);
    end

    fprintf('  After rejection, %d trials remain.\n', EEG.trials);
    parallel_results.removed_trials(file_idx) = length(final_bad_trials);
    parallel_results.final_trials(file_idx) = EEG.trials;

    subject_stats.SME_final = subject_stats.SME_original - subject_stats.SME_removed;
    subject_stats.Test_Intact_final = subject_stats.Test_Intact_original - subject_stats.Test_Intact_removed;
    subject_stats.Test_Recombined_final = subject_stats.Test_Recombined_original - subject_stats.Test_Recombined_removed;

    % Count final trials per condition
    subject_final_condition_counts = struct();
    subject_final_condition_counts.Study_hits = 0;
    subject_final_condition_counts.Study_misses = 0;
    subject_final_condition_counts.Test_hits = 0;
    subject_final_condition_counts.Test_misses = 0;
    subject_final_condition_counts.Correct_rejections = 0;
    subject_final_condition_counts.False_alarms = 0;

    final_epoch_condition_map = cell(1, EEG.trials);
    for i = 1:length(EEG.event)
        if ~isempty(EEG.event(i).trial_type)
            epoch_num = EEG.event(i).epoch;
            trial_type = EEG.event(i).trial_type;
            if isempty(final_epoch_condition_map{epoch_num})
                final_epoch_condition_map{epoch_num} = trial_type;
                subject_final_condition_counts.(trial_type) = subject_final_condition_counts.(trial_type) + 1;
            end
        end
    end

    parallel_results.cond_final_Study_hits(file_idx) = subject_final_condition_counts.Study_hits;
    parallel_results.cond_final_Study_misses(file_idx) = subject_final_condition_counts.Study_misses;
    parallel_results.cond_final_Test_hits(file_idx) = subject_final_condition_counts.Test_hits;
    parallel_results.cond_final_Test_misses(file_idx) = subject_final_condition_counts.Test_misses;
    parallel_results.cond_final_Correct_rejections(file_idx) = subject_final_condition_counts.Correct_rejections;
    parallel_results.cond_final_False_alarms(file_idx) = subject_final_condition_counts.False_alarms;

    parallel_results.cond_removed_Study_hits(file_idx) = subject_condition_counts.Study_hits - subject_final_condition_counts.Study_hits;
    parallel_results.cond_removed_Study_misses(file_idx) = subject_condition_counts.Study_misses - subject_final_condition_counts.Study_misses;
    parallel_results.cond_removed_Test_hits(file_idx) = subject_condition_counts.Test_hits - subject_final_condition_counts.Test_hits;
    parallel_results.cond_removed_Test_misses(file_idx) = subject_condition_counts.Test_misses - subject_final_condition_counts.Test_misses;
    parallel_results.cond_removed_Correct_rejections(file_idx) = subject_condition_counts.Correct_rejections - subject_final_condition_counts.Correct_rejections;
    parallel_results.cond_removed_False_alarms(file_idx) = subject_condition_counts.False_alarms - subject_final_condition_counts.False_alarms;

    % --- STEP 6: Final validation and save ---
    EEG = eeg_checkset(EEG, 'eventconsistency');
    EEG.setname = sprintf('Subject %d - All Conditions Clean', s);
    EEG.filename = sprintf('%d_epoch.set', s);

    fprintf('  Saving final clean dataset for subject %d...\n', s);
    pop_saveset(EEG, 'filename', EEG.filename, 'filepath', epoch_folder, 'savemode', 'twofiles');
    parallel_results.processed(file_idx) = true;

    parallel_results.condition1_removals(file_idx) = subject_stats.condition1_total;
    parallel_results.condition2_removals(file_idx) = subject_stats.condition2_total;
    parallel_results.SME_original(file_idx) = subject_stats.SME_original;
    parallel_results.SME_final(file_idx) = subject_stats.SME_final;
    parallel_results.Test_Intact_original(file_idx) = subject_stats.Test_Intact_original;
    parallel_results.Test_Intact_final(file_idx) = subject_stats.Test_Intact_final;
    parallel_results.Test_Recombined_original(file_idx) = subject_stats.Test_Recombined_original;
    parallel_results.Test_Recombined_final(file_idx) = subject_stats.Test_Recombined_final;
end

% ==================== AGGREGATE RESULTS ====================
fprintf('\nAggregating results...\n');

stats.processed_subjects = sum(parallel_results.processed);
stats.excluded_subjects = sum(parallel_results.excluded);

excluded_indices = find(parallel_results.excluded);
stats.excluded_subject_ids = parallel_results.subject_id(excluded_indices);
stats.exclusion_reasons = parallel_results.exclusion_reason(excluded_indices);

processed_not_excluded = parallel_results.processed & ~parallel_results.excluded & ~parallel_results.skipped;
stats.total_original_trials = sum(parallel_results.original_trials(processed_not_excluded));
stats.total_final_trials = sum(parallel_results.final_trials(processed_not_excluded));
stats.total_removed_trials = sum(parallel_results.removed_trials(processed_not_excluded));
stats.condition1_removals = parallel_results.condition1_removals(processed_not_excluded);
stats.condition2_removals = parallel_results.condition2_removals(processed_not_excluded);
stats.total_removals_per_subject = parallel_results.removed_trials(processed_not_excluded);
stats.SME_removals = parallel_results.SME_removals(processed_not_excluded);
stats.SME_original_trials = parallel_results.SME_original(processed_not_excluded);
stats.SME_final_trials = parallel_results.SME_final(processed_not_excluded);
stats.Test_Intact_removals = parallel_results.Test_Intact_removals(processed_not_excluded);
stats.Test_Intact_original_trials = parallel_results.Test_Intact_original(processed_not_excluded);
stats.Test_Intact_final_trials = parallel_results.Test_Intact_final(processed_not_excluded);
stats.Test_Recombined_removals = parallel_results.Test_Recombined_removals(processed_not_excluded);
stats.Test_Recombined_original_trials = parallel_results.Test_Recombined_original(processed_not_excluded);
stats.Test_Recombined_final_trials = parallel_results.Test_Recombined_final(processed_not_excluded);
for cond = condition_names
    stats.condition_original_trials.(cond{1}) = parallel_results.(['cond_orig_' cond{1}])(processed_not_excluded);
    stats.condition_final_trials.(cond{1}) = parallel_results.(['cond_final_' cond{1}])(processed_not_excluded);
    stats.condition_removed_trials.(cond{1}) = parallel_results.(['cond_removed_' cond{1}])(processed_not_excluded);
end

valid_subject_ids = parallel_results.subject_id(processed_not_excluded);
stats.subject_ids = cellfun(@(x) str2double(x), valid_subject_ids);
stats.subject_original_trials = parallel_results.original_trials(processed_not_excluded);
stats.subject_final_trials = parallel_results.final_trials(processed_not_excluded);
stats.subject_removed_trials = parallel_results.removed_trials(processed_not_excluded);

% ==================== DETAILED STATISTICS SUMMARY ====================
fprintf('\n==================== DETAILED STATISTICS SUMMARY ====================\n');
fprintf('Subject Processing Summary:\n');
fprintf('  Total subjects selected: %d\n', stats.total_subjects);
fprintf('  Successfully processed: %d\n', stats.processed_subjects);
fprintf('  Excluded subjects: %d\n', stats.excluded_subjects);
if stats.excluded_subjects > 0
    fprintf('  Excluded subject IDs: %s\n', strjoin(stats.excluded_subject_ids, ', '));
    fprintf('  Exclusion reasons:\n');
    for i = 1:length(stats.exclusion_reasons)
        fprintf('    %s\n', stats.exclusion_reasons{i});
    end
end
fprintf('  Processing success rate: %.1f%%\n', (stats.processed_subjects/stats.total_subjects)*100);

fprintf('\nOverall Trial Statistics:\n');
fprintf('  Total original trials: %d\n', stats.total_original_trials);
fprintf('  Total final trials: %d\n', stats.total_final_trials);
fprintf('  Total removed trials: %d\n', stats.total_removed_trials);
if stats.total_original_trials > 0
    fprintf('  Overall removal rate: %.1f%%\n', (stats.total_removed_trials/stats.total_original_trials)*100);
end

if stats.processed_subjects > 0
    fprintf('\nThreshold Parameters Used:\n');
    fprintf('  Voltage difference threshold: %.1f uV\n', stats.voltage_diff_threshold);
    fprintf('  Absolute voltage threshold:   %.1f uV\n', stats.voltage_abs_threshold);

    fprintf('\nData Quality Assessment:\n');
    avg_retention_rate = (stats.total_final_trials / stats.total_original_trials) * 100;
    fprintf('  Overall data retention rate: %.1f%%\n', avg_retention_rate);
    if avg_retention_rate >= 85
        fprintf('  Data quality: EXCELLENT (>85%% retention)\n');
    elseif avg_retention_rate >= 75
        fprintf('  Data quality: GOOD (75-85%% retention)\n');
    elseif avg_retention_rate >= 65
        fprintf('  Data quality: ACCEPTABLE (65-75%% retention)\n');
    else
        fprintf('  Data quality: CONCERNING (<65%% retention)\n');
    end
end

fprintf('\n==================== PROCESSING COMPLETE ====================\n');

end


%% ================== HELPER FUNCTIONS ==================

function [clean_EEG, bad_trials, condition1_bad, condition2_bad] = get_bad_trials_v2(EEG, voltage_diff_threshold, voltage_abs_threshold, all_channels)
    bad_trials_condition1 = [];
    bad_trials_condition2 = [];
    for epoch_idx = 1:size(EEG.data, 3)
        for ch = all_channels
            voltage_data = squeeze(EEG.data(ch, :, epoch_idx));
            if max(abs(diff(voltage_data))) > voltage_diff_threshold
                bad_trials_condition1 = [bad_trials_condition1, epoch_idx];
            end
            if max(abs(voltage_data)) > voltage_abs_threshold
                bad_trials_condition2 = [bad_trials_condition2, epoch_idx];
            end
        end
    end
    condition1_bad = unique(bad_trials_condition1);
    condition2_bad = unique(bad_trials_condition2);
    bad_trials = unique([condition1_bad, condition2_bad]);
    clean_EEG = EEG;
end

function trial_indices_by_group = assign_epochs_to_groups_fixed(EEG)
    trial_indices_by_group = struct();
    trial_indices_by_group.SME = [];
    trial_indices_by_group.Test_Intact = [];
    trial_indices_by_group.Test_Recombined = [];

    num_events = length(EEG.event);
    event_epochs = zeros(1, num_events);
    event_latencies = zeros(1, num_events);
    event_types_cell = cell(1, num_events);

    for ev_idx = 1:num_events
        event_epochs(ev_idx) = EEG.event(ev_idx).epoch;
        event_latencies(ev_idx) = EEG.event(ev_idx).latency;
        if isnumeric(EEG.event(ev_idx).type)
            event_types_cell{ev_idx} = num2str(EEG.event(ev_idx).type);
        else
            event_types_cell{ev_idx} = EEG.event(ev_idx).type;
        end
    end

    epoch_to_trigger = zeros(1, EEG.trials);
    for epoch_num = 1:EEG.trials
        epoch_event_indices = find(event_epochs == epoch_num);
        if ~isempty(epoch_event_indices)
            triggering_event_idx = [];
            earliest_latency = inf;
            for idx = epoch_event_indices
                event_type_str = event_types_cell{idx};
                if ismember(event_type_str, {'11', '21', '22'})
                    if event_latencies(idx) < earliest_latency
                        earliest_latency = event_latencies(idx);
                        triggering_event_idx = idx;
                    end
                end
            end
            if ~isempty(triggering_event_idx)
                epoch_to_trigger(epoch_num) = str2double(event_types_cell{triggering_event_idx});
            end
        end
    end

    trial_indices_by_group.SME = find(epoch_to_trigger == 11);
    trial_indices_by_group.Test_Intact = find(epoch_to_trigger == 21);
    trial_indices_by_group.Test_Recombined = find(epoch_to_trigger == 22);
end
