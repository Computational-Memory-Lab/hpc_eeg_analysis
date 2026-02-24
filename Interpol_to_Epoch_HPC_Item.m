%% EEG DATA EPOCHING AND PROCESSING SCRIPT - VERSION 3 ITEM (Robust Workflow)
% This script processes EEG data files with condition-specific trial removal
% for the legacy OldNew item-recognition structure (type/testacc/testcresp).
% It uses a robust, standard EEGLAB workflow that avoids data corruption.
%
% WORKFLOW:
% 1. Epoch all trials for all conditions of interest together into one dataset.
% 2. Add a custom 'trial_type' field to each event based on its properties.
% 3. Identify the trial indices corresponding to each experimental group for artifact rejection.
% 4. For each group, detect artifacts and collect the indices of "bad" trials.
% 5. After all groups are checked, remove all bad trials from the main dataset in a single step.
% 6. Save the resulting clean, consistent, and valid dataset.

% ==================== INITIALIZE EEGLAB (HEADLESS MODE) ====================
% Add EEGLAB to path and initialize without GUI
eeglab_path = '/home/devon7y/scratch/devon7y/eeglab2022.1';
addpath(eeglab_path);
eeglab nogui;
fprintf('EEGLAB initialized in headless mode\n');

% ==================== FILE SELECTION (HEADLESS MODE - MERGED INPUT) ====================
% Optional environment override:
%   INPUT_FOLDER=/absolute/path/to/folder
input_folder_override = strtrim(getenv('INPUT_FOLDER'));
if ~isempty(input_folder_override)
    input_paths = {input_folder_override};
else
    input_paths = {'/home/devon7y/scratch/devon7y/combined/interpol_merged_renum2_20260221_231258'};
end
input_pattern = '*.set'; % Match all .set files

% Since we're processing unclustered data, set clustering type
clustering_type = 'unclustered';
fprintf('Processing unclustered data (no clustering applied)\n');
fprintf('Input root(s): %s\n', strjoin(input_paths, ', '));

% Collect all .set files from both directories
input_files = struct('name', {}, 'folder', {}, 'date', {}, 'bytes', {}, 'isdir', {}, 'datenum', {});
all_files_count = 0;

for path_idx = 1:length(input_paths)
    current_path = input_paths{path_idx};
    fprintf('Scanning directory: %s\n', current_path);

    % Find all .set files in the current directory
    files_in_dir = dir(fullfile(current_path, input_pattern));
    all_files_count = all_files_count + length(files_in_dir);

    % Filter out LIMO-generated files and already-epoched files
    for i = 1:length(files_in_dir)
        filename = files_in_dir(i).name;
        % Exclude files that start with "sub-" or contain "_ses-" or "_design" or "_epoch"
        if ~startsWith(filename, 'sub-') && ~contains(filename, '_ses-') && ...
           ~contains(filename, '_design') && ~contains(filename, '_epoch')
            input_files(end+1) = files_in_dir(i);
        end
    end
end

fprintf('Found %d .set files total (%d after excluding LIMO-generated and epoch files)\n', all_files_count, length(input_files));

% Create a file list with study tracking
file_info = struct('filename', {}, 'filepath', {}, 'subject_id', {}, 'study_id', {}, 'clustering_type', {});

% Process all files and extract subject ID and study name
for i = 1:length(input_files)
    filename = input_files(i).name;
    filepath = input_files(i).folder; % Use the actual folder from the file struct

    % Extract first numeric token from filename as subject ID.
    id_tokens = regexp(filename, '(\d+)', 'tokens', 'once');
    if isempty(id_tokens)
        fprintf('Skipping file without numeric subject token: %s\n', filename);
        continue;
    end
    subject_id = str2double(id_tokens{1});
    if isnan(subject_id)
        fprintf('Skipping file with invalid subject token: %s\n', filename);
        continue;
    end

    % Infer source study from filename/path conventions used in this merged set.
    lower_filename = lower(filename);
    lower_filepath = lower(filepath);
    if contains(lower_filename, 'interpol_resampled') || subject_id >= 2000 || contains(lower_filepath, 'tamari')
        study_name = 'tamari';
    elseif contains(lower_filename, 'no_reference') || contains(lower_filepath, 'yvonne')
        study_name = 'yvonne';
    else
        study_name = 'unknown';
    end

    file_info(end+1).filename = filename;
    file_info(end).filepath = filepath;
    file_info(end).subject_id = subject_id;
    file_info(end).study_id = study_name;
    file_info(end).clustering_type = clustering_type;
end

if isempty(file_info)
    error('No .set files found matching expected patterns in input directories');
end

subject_ids = [file_info.subject_id];
[unique_ids, ~, ic] = unique(subject_ids);
id_counts = accumarray(ic, 1);
dup_ids = unique_ids(id_counts > 1);
if ~isempty(dup_ids)
    dup_text = strjoin(arrayfun(@(x) num2str(x), dup_ids, 'UniformOutput', false), ',');
    error('Duplicate subject IDs detected after parsing: %s', dup_text);
end

fprintf('Total files to process: %d\n', length(file_info));

% ==================== CREATE OUTPUT DIRECTORY ====================
% Optional environment override:
%   OUTPUT_FOLDER=/absolute/path/to/output/folder
output_folder_override = strtrim(getenv('OUTPUT_FOLDER'));
if ~isempty(output_folder_override)
    epoched_dir = output_folder_override;
else
    epoched_dir = fullfile(input_paths{1}, 'epoch');
end
if ~exist(epoched_dir, 'dir'), mkdir(epoched_dir); end
fprintf('Output will be saved to: %s\n', epoched_dir);

% ==================== INITIALIZE STATISTICS TRACKING ====================
stats = struct();
stats.total_subjects = length(file_info);
stats.processed_subjects = 0;
stats.excluded_subjects = 0;
stats.excluded_subject_ids = {}; % cell array to store 'subjectID_study' format
stats.voltage_diff_threshold = 20;
stats.voltage_abs_threshold = 1000;
stats.max_rejected_trials_per_group = 56;
stats.total_original_trials = 0;
stats.total_final_trials = 0;
stats.total_removed_trials = 0;
stats.SME_removals = [];
stats.Test_Intact_removals = [];
stats.Test_Recombined_removals = [];

% Condition-based removal statistics
stats.condition1_removals = []; % voltage difference removals per subject
stats.condition2_removals = []; % absolute voltage removals per subject
stats.total_removals_per_subject = []; % total removals per subject

% Group trial count statistics
stats.SME_original_trials = []; % original trials per subject
stats.SME_final_trials = []; % final trials per subject
stats.Test_Intact_original_trials = [];
stats.Test_Intact_final_trials = [];
stats.Test_Recombined_original_trials = [];
stats.Test_Recombined_final_trials = [];

% NEW: Subject-based statistics (per processed subject)
stats.subject_ids = []; % track which subjects were processed
stats.subject_original_trials = []; % total original trials per subject
stats.subject_final_trials = []; % total final trials per subject
stats.subject_removed_trials = []; % total removed trials per subject

% NEW: Trial type statistics (6 conditions)
stats.condition_original_trials = struct(); % trials per condition across all subjects
stats.condition_final_trials = struct();
stats.condition_removed_trials = struct();
condition_names = {'Study_hits', 'Study_misses', 'Test_hits', 'Test_misses', 'Correct_rejections', 'False_alarms'};
for cond = condition_names
    stats.condition_original_trials.(cond{1}) = [];
    stats.condition_final_trials.(cond{1}) = [];
    stats.condition_removed_trials.(cond{1}) = [];
end

% Exclusion tracking with reasons
stats.exclusion_reasons = {}; % cell array storing why each subject was excluded

% ==================== PRE-ALLOCATE PARALLEL RESULTS STORAGE ====================
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

% Condition counts (6 conditions)
for cond = condition_names
    parallel_results.(['cond_orig_' cond{1}]) = zeros(num_files, 1);
    parallel_results.(['cond_final_' cond{1}]) = zeros(num_files, 1);
    parallel_results.(['cond_removed_' cond{1}]) = zeros(num_files, 1);
end

% ==================== PROCESS EACH SUBJECT (SERIAL) ====================
for file_idx = 1:length(file_info)
    % Get current file information
    current_file = file_info(file_idx);
    s = current_file.subject_id;
    study_id = current_file.study_id;
    cluster_type = current_file.clustering_type;
    input_filename = current_file.filename;
    input_filepath = current_file.filepath;

    % Use canonical filename expected by hpc_limo_first_level.m
    output_filename = sprintf('%d_epoch.set', s);
    output_filepath = fullfile(epoched_dir, output_filename);

    % Store subject identifier
    parallel_results.subject_id{file_idx} = sprintf('%d_%s', s, study_id);

    if exist(output_filepath, 'file')
        fprintf('\nSubject %d (%s study) already processed (file exists: %s). Skipping...\n', s, study_id, output_filename);
        parallel_results.skipped(file_idx) = true;
        parallel_results.processed(file_idx) = true;
        continue;
    end

    EEG = pop_loadset('filename', input_filename, 'filepath', input_filepath);
    fprintf('\nProcessing subject %d (%s study)...\n', s, study_id);
    
    % Clear any existing trial_type field before processing
    if isfield(EEG.event, 'trial_type')
        EEG.event = rmfield(EEG.event, 'trial_type');
    end

    % =========================================================================
    % ===== NEW ROBUST WORKFLOW STARTS HERE ===================================
    % =========================================================================

    % --- STEP 1: Epoch all trials for all conditions of interest together ---
    % Legacy item-recognition structure:
    %   study onset = type 2
    %   test onset  = type 7
    all_triggers = {'2', '7'};
    EEG = pop_epoch(EEG, all_triggers, [-0.1 1.5], 'newname', sprintf('Subject %d (%s) - all epoched', s, study_id), 'epochinfo', 'yes');
    EEG = pop_rmbase(EEG, [-100 0]);
    fprintf('  Epoched all %d trials together.\n', EEG.trials);
    parallel_results.original_trials(file_idx) = EEG.trials;

    % --- STEP 2: Add the custom 'trial_type' field based on 6 conditions ---
    fprintf('  Adding custom ''trial_type'' field to events...\n');

    % Initialize the new field for all events
    [EEG.event.trial_type] = deal('');

    % Count original trials per condition for this subject
    subject_condition_counts = struct();
    subject_condition_counts.Study_hits = 0;
    subject_condition_counts.Study_misses = 0;
    subject_condition_counts.Test_hits = 0;
    subject_condition_counts.Test_misses = 0;
    subject_condition_counts.Correct_rejections = 0;
    subject_condition_counts.False_alarms = 0;

    % Per-epoch trigger summary for classification.
    [epoch_trigger_type, epoch_trigger_latency, epoch_testacc, epoch_testcresp] = summarize_epoch_trigger_info_item(EEG);
    epoch_trial_type = repmat({''}, 1, EEG.trials);

    % First pass: classify test epochs directly from testacc/testcresp.
    for epoch_num = 1:EEG.trials
        if strcmp(epoch_trigger_type{epoch_num}, '7')
            cond_name = classify_test_condition(epoch_testacc(epoch_num), epoch_testcresp(epoch_num));
            if ~isempty(cond_name)
                epoch_trial_type{epoch_num} = cond_name;
            end
        end
    end

    % Second pass: define Study_hits/Study_misses.
    % We map study onsets (type 2) to target-old test outcomes
    % (type 7 with testcresp==1) by chronological order.
    study_epochs = find(strcmp(epoch_trigger_type, '2'));
    target_test_epochs = find(strcmp(epoch_trigger_type, '7') & epoch_testcresp == 1 & ...
                              ismember(epoch_trial_type, {'Test_hits', 'Test_misses'}));

    [~, study_sort_idx] = sort(epoch_trigger_latency(study_epochs), 'ascend');
    [~, test_sort_idx] = sort(epoch_trigger_latency(target_test_epochs), 'ascend');
    study_epochs = study_epochs(study_sort_idx);
    target_test_epochs = target_test_epochs(test_sort_idx);

    if numel(study_epochs) ~= numel(target_test_epochs)
        fprintf('  WARNING: Study/Test target count mismatch (study=%d, target test=%d). Pairing by min count.\n', ...
                numel(study_epochs), numel(target_test_epochs));
    end

    n_pairs = min(numel(study_epochs), numel(target_test_epochs));
    for k = 1:n_pairs
        paired_test_type = epoch_trial_type{target_test_epochs(k)};
        if strcmp(paired_test_type, 'Test_hits')
            epoch_trial_type{study_epochs(k)} = 'Study_hits';
        elseif strcmp(paired_test_type, 'Test_misses')
            epoch_trial_type{study_epochs(k)} = 'Study_misses';
        end
    end

    % Propagate epoch-level trial_type onto all events in that epoch.
    for i = 1:length(EEG.event)
        if isfield(EEG.event(i), 'epoch')
            epoch_num = EEG.event(i).epoch;
            if epoch_num >= 1 && epoch_num <= EEG.trials
                EEG.event(i).trial_type = epoch_trial_type{epoch_num};
            end
        end
    end

    % Count original trials per condition (one count per epoch).
    for epoch_num = 1:EEG.trials
        cond_name = epoch_trial_type{epoch_num};
        if ~isempty(cond_name)
            subject_condition_counts.(cond_name) = subject_condition_counts.(cond_name) + 1;
        end
    end

    % Store original condition counts for this subject
    parallel_results.cond_orig_Study_hits(file_idx) = subject_condition_counts.Study_hits;
    parallel_results.cond_orig_Study_misses(file_idx) = subject_condition_counts.Study_misses;
    parallel_results.cond_orig_Test_hits(file_idx) = subject_condition_counts.Test_hits;
    parallel_results.cond_orig_Test_misses(file_idx) = subject_condition_counts.Test_misses;
    parallel_results.cond_orig_Correct_rejections(file_idx) = subject_condition_counts.Correct_rejections;
    parallel_results.cond_orig_False_alarms(file_idx) = subject_condition_counts.False_alarms;

    % --- STEP 3: Identify trial indices for each artifact rejection group ---
    fprintf('  Identifying trials for artifact rejection groups...\n');

    % Call helper function (defined at end of file) - simpler for parfor
    trial_indices_by_group = assign_epochs_to_groups_fixed(EEG);

    fprintf('  Identified trials per group: SME(%d), Test_Intact(%d), Test_Recombined(%d)\n', ...
        length(trial_indices_by_group.SME), length(trial_indices_by_group.Test_Intact), length(trial_indices_by_group.Test_Recombined));

    % VALIDATION: Check that each epoch belongs to exactly one group
    total_assigned = length(trial_indices_by_group.SME) + length(trial_indices_by_group.Test_Intact) + length(trial_indices_by_group.Test_Recombined);
    if total_assigned ~= EEG.trials
        fprintf('  WARNING: Group assignment mismatch! Total epochs: %d, Assigned: %d\n', EEG.trials, total_assigned);
    else
        fprintf('  VALIDATION PASSED: Each epoch assigned to exactly one group (Total: %d)\n', total_assigned);
    end

    % --- STEP 4: Find bad trials per group without deleting them yet ---
    all_bad_trials_absolute = [];
    subject_should_be_skipped = false;
    exclusion_reason = '';
    all_channels = 1:EEG.nbchan;

    % Initialize per-subject tracking
    subject_stats = struct();
    subject_stats.condition1_total = 0; % total trials removed by voltage difference
    subject_stats.condition2_total = 0; % total trials removed by absolute voltage
    subject_stats.SME_original = 0;
    subject_stats.SME_removed = 0;
    subject_stats.SME_final = 0;
    subject_stats.Test_Intact_original = 0;
    subject_stats.Test_Intact_removed = 0;
    subject_stats.Test_Intact_final = 0;
    subject_stats.Test_Recombined_original = 0;
    subject_stats.Test_Recombined_removed = 0;
    subject_stats.Test_Recombined_final = 0;
    total_removed_trials = 0;

    % Define the 3 artifact rejection groups (by derived trial_type)
    group_def = {'SME';
                 'Test_Intact';
                 'Test_Recombined'};

    for g = 1:size(group_def, 1)
        group_name = group_def{g, 1};
        group_trial_indices = trial_indices_by_group.(group_name);

        if isempty(group_trial_indices)
            fprintf('    No trials found for group %s, skipping artifact detection.\n', group_name);
            continue;
        end

        % Create a temporary dataset for this group only
        EEG_group = pop_select(EEG, 'trial', group_trial_indices);

        % Store original trial count for this group
        group_original_trials = length(group_trial_indices);
        if strcmp(group_name, 'SME')
            subject_stats.SME_original = group_original_trials;
        elseif strcmp(group_name, 'Test_Intact')
            subject_stats.Test_Intact_original = group_original_trials;
        elseif strcmp(group_name, 'Test_Recombined')
            subject_stats.Test_Recombined_original = group_original_trials;
        end

        % Get bad trials (indices are relative to the temporary EEG_group)
        [~, bad_trials_relative, cond1_bad, cond2_bad] = get_bad_trials_v2(EEG_group, stats.voltage_diff_threshold, stats.voltage_abs_threshold, all_channels);
        
        fprintf('    Found %d bad trials in %s group.\n', length(bad_trials_relative), group_name);

        % Track condition-specific removals for this subject
        subject_stats.condition1_total = subject_stats.condition1_total + length(cond1_bad);
        subject_stats.condition2_total = subject_stats.condition2_total + length(cond2_bad);

        % Check for exclusion
        if length(bad_trials_relative) > stats.max_rejected_trials_per_group
            fprintf('    WARNING: Subject %d (%s study) excluded. Too many trials rejected in %s group (%d > %d).\n', s, study_id, group_name, length(bad_trials_relative), stats.max_rejected_trials_per_group);
            subject_should_be_skipped = true;
            exclusion_reason = sprintf('Too many trials rejected in %s group (%d > %d)', group_name, length(bad_trials_relative), stats.max_rejected_trials_per_group);
        end

        % IMPORTANT: Map relative indices back to absolute indices in the main EEG set
        bad_trials_absolute = group_trial_indices(bad_trials_relative);
        all_bad_trials_absolute = [all_bad_trials_absolute(:); bad_trials_absolute(:)];
        
        % Store group-specific stats
        if strcmp(group_name, 'SME')
            subject_stats.SME_removed = length(bad_trials_relative);
        elseif strcmp(group_name, 'Test_Intact')
            subject_stats.Test_Intact_removed = length(bad_trials_relative);
        elseif strcmp(group_name, 'Test_Recombined')
            subject_stats.Test_Recombined_removed = length(bad_trials_relative);
        end
        
        % Store stats in parallel_results
        parallel_results.([group_name '_removals'])(file_idx) = length(bad_trials_relative);
    end

    if subject_should_be_skipped
        parallel_results.excluded(file_idx) = true;
        parallel_results.exclusion_reason{file_idx} = sprintf('Subject %d (%s): %s', s, study_id, exclusion_reason);

        % Store zero final trials and all original as removed for excluded subjects
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
        continue; % Skip to the next subject
    end

    % --- STEP 5: Reject all bad trials from all groups at once ---
    final_bad_trials = unique(all_bad_trials_absolute);
    fprintf('  Total unique bad trials to remove: %d\n', length(final_bad_trials));
    
    if ~isempty(final_bad_trials)
        EEG = pop_select(EEG, 'notrial', final_bad_trials);
    end
    
    fprintf('  After rejection, %d trials remain.\n', EEG.trials);
    total_removed_trials = length(final_bad_trials);
    parallel_results.removed_trials(file_idx) = total_removed_trials;
    parallel_results.final_trials(file_idx) = EEG.trials;

    % Update final trial counts for each group
    subject_stats.SME_final = subject_stats.SME_original - subject_stats.SME_removed;
    subject_stats.Test_Intact_final = subject_stats.Test_Intact_original - subject_stats.Test_Intact_removed;
    subject_stats.Test_Recombined_final = subject_stats.Test_Recombined_original - subject_stats.Test_Recombined_removed;

    % NEW: Count final trials per condition for this subject
    subject_final_condition_counts = struct();
    subject_final_condition_counts.Study_hits = 0;
    subject_final_condition_counts.Study_misses = 0;
    subject_final_condition_counts.Test_hits = 0;
    subject_final_condition_counts.Test_misses = 0;
    subject_final_condition_counts.Correct_rejections = 0;
    subject_final_condition_counts.False_alarms = 0;

    % Track which epochs have been counted per condition (to avoid double-counting)
    final_epoch_condition_map = cell(1, EEG.trials);

    for i = 1:length(EEG.event)
        if ~isempty(EEG.event(i).trial_type)
            epoch_num = EEG.event(i).epoch;
            trial_type = EEG.event(i).trial_type;

            % Only count each epoch once per condition
            if isempty(final_epoch_condition_map{epoch_num})
                final_epoch_condition_map{epoch_num} = trial_type;
                subject_final_condition_counts.(trial_type) = subject_final_condition_counts.(trial_type) + 1;
            end
        end
    end

    % Store final condition counts and calculate removed counts for this subject
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
    EEG.setname = sprintf('Subject %d (%s) - All Conditions Clean', s, study_id);
    EEG.filename = sprintf('%d_epoch.set', s);

    fprintf('  Saving final clean dataset for subject %d (%s study)...\n', s, study_id);
    pop_saveset(EEG, 'filename', EEG.filename, 'filepath', epoched_dir, 'savemode', 'twofiles');
    parallel_results.processed(file_idx) = true;

    % Store per-subject statistics in parallel_results
    parallel_results.condition1_removals(file_idx) = subject_stats.condition1_total;
    parallel_results.condition2_removals(file_idx) = subject_stats.condition2_total;

    % Store group-based statistics in parallel_results
    parallel_results.SME_original(file_idx) = subject_stats.SME_original;
    parallel_results.SME_final(file_idx) = subject_stats.SME_final;
    parallel_results.Test_Intact_original(file_idx) = subject_stats.Test_Intact_original;
    parallel_results.Test_Intact_final(file_idx) = subject_stats.Test_Intact_final;
    parallel_results.Test_Recombined_original(file_idx) = subject_stats.Test_Recombined_original;
    parallel_results.Test_Recombined_final(file_idx) = subject_stats.Test_Recombined_final;

    % =========================================================================
    % ===== END OF NEW WORKFLOW ===============================================
    % =========================================================================
end

% ==================== AGGREGATE PARALLEL RESULTS ====================
fprintf('\nAggregating results from parallel processing...\n');

% Count processed, excluded, and skipped subjects
stats.processed_subjects = sum(parallel_results.processed);
stats.excluded_subjects = sum(parallel_results.excluded);

% Collect exclusion info
excluded_indices = find(parallel_results.excluded);
stats.excluded_subject_ids = parallel_results.subject_id(excluded_indices);
stats.exclusion_reasons = parallel_results.exclusion_reason(excluded_indices);

% Aggregate trial counts (only from processed, non-excluded subjects)
processed_not_excluded = parallel_results.processed & ~parallel_results.excluded & ~parallel_results.skipped;
stats.total_original_trials = sum(parallel_results.original_trials(processed_not_excluded));
stats.total_final_trials = sum(parallel_results.final_trials(processed_not_excluded));
stats.total_removed_trials = sum(parallel_results.removed_trials(processed_not_excluded));

% Aggregate condition-based removals
stats.condition1_removals = parallel_results.condition1_removals(processed_not_excluded);
stats.condition2_removals = parallel_results.condition2_removals(processed_not_excluded);
stats.total_removals_per_subject = parallel_results.removed_trials(processed_not_excluded);

% Aggregate group-based statistics
stats.SME_removals = parallel_results.SME_removals(processed_not_excluded);
stats.SME_original_trials = parallel_results.SME_original(processed_not_excluded);
stats.SME_final_trials = parallel_results.SME_final(processed_not_excluded);
stats.Test_Intact_removals = parallel_results.Test_Intact_removals(processed_not_excluded);
stats.Test_Intact_original_trials = parallel_results.Test_Intact_original(processed_not_excluded);
stats.Test_Intact_final_trials = parallel_results.Test_Intact_final(processed_not_excluded);
stats.Test_Recombined_removals = parallel_results.Test_Recombined_removals(processed_not_excluded);
stats.Test_Recombined_original_trials = parallel_results.Test_Recombined_original(processed_not_excluded);
stats.Test_Recombined_final_trials = parallel_results.Test_Recombined_final(processed_not_excluded);

% Aggregate condition statistics
for cond = condition_names
    stats.condition_original_trials.(cond{1}) = parallel_results.(['cond_orig_' cond{1}])(processed_not_excluded);
    stats.condition_final_trials.(cond{1}) = parallel_results.(['cond_final_' cond{1}])(processed_not_excluded);
    stats.condition_removed_trials.(cond{1}) = parallel_results.(['cond_removed_' cond{1}])(processed_not_excluded);
end

% Aggregate subject-based statistics
valid_subject_ids = parallel_results.subject_id(processed_not_excluded);
stats.subject_ids = cellfun(@(x) str2double(regexp(x, '^\d+', 'match', 'once')), valid_subject_ids);
stats.subject_original_trials = parallel_results.original_trials(processed_not_excluded);
stats.subject_final_trials = parallel_results.final_trials(processed_not_excluded);
stats.subject_removed_trials = parallel_results.removed_trials(processed_not_excluded);

fprintf('Aggregation complete.\n');

% ==================== DETAILED STATISTICS SUMMARY ====================
fprintf('\n==================== DETAILED STATISTICS SUMMARY ====================\n');

% Subject processing summary
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

% Overall trial statistics
fprintf('\nOverall Trial Statistics:\n');
fprintf('  Total original trials: %d\n', stats.total_original_trials);
fprintf('  Total final trials: %d\n', stats.total_final_trials);
fprintf('  Total removed trials: %d\n', stats.total_removed_trials);
fprintf('  Overall removal rate: %.1f%%\n', (stats.total_removed_trials/stats.total_original_trials)*100);

% Condition-based removal statistics (only for processed subjects)
if stats.processed_subjects > 0
    fprintf('\nCondition-Based Removal Statistics (averaged across %d processed subjects):\n', stats.processed_subjects);
    fprintf('  Condition 1 (Voltage Difference > %.1f μV):\n', stats.voltage_diff_threshold);
    fprintf('    Average trials removed per subject: %.1f (±%.1f)\n', mean(stats.condition1_removals), std(stats.condition1_removals));
    fprintf('    Range: %d - %d trials\n', min(stats.condition1_removals), max(stats.condition1_removals));
    fprintf('    Total trials removed by condition 1: %d\n', sum(stats.condition1_removals));
    
    fprintf('  Condition 2 (Absolute Voltage > %.1f μV):\n', stats.voltage_abs_threshold);
    fprintf('    Average trials removed per subject: %.1f (±%.1f)\n', mean(stats.condition2_removals), std(stats.condition2_removals));
    fprintf('    Range: %d - %d trials\n', min(stats.condition2_removals), max(stats.condition2_removals));
    fprintf('    Total trials removed by condition 2: %d\n', sum(stats.condition2_removals));
    
    % Note about overlap
    total_cond1_cond2 = sum(stats.condition1_removals) + sum(stats.condition2_removals);
    overlap_trials = total_cond1_cond2 - stats.total_removed_trials;
    fprintf('  Overlap (trials meeting both conditions): %d\n', overlap_trials);
    
    % Group-based statistics
    fprintf('\nGroup-Based Trial Statistics (averaged across %d processed subjects):\n', stats.processed_subjects);
    
    fprintf('  SME Group (Study Memory Effects):\n');
    fprintf('    Average original trials per subject: %.1f (±%.1f)\n', mean(stats.SME_original_trials), std(stats.SME_original_trials));
    fprintf('    Average trials removed per subject: %.1f (±%.1f)\n', mean(stats.SME_removals), std(stats.SME_removals));
    fprintf('    Average final trials per subject: %.1f (±%.1f)\n', mean(stats.SME_final_trials), std(stats.SME_final_trials));
    fprintf('    Group removal rate: %.1f%%\n', (sum(stats.SME_removals)/sum(stats.SME_original_trials))*100);
    
    fprintf('  Test_Intact Group:\n');
    fprintf('    Average original trials per subject: %.1f (±%.1f)\n', mean(stats.Test_Intact_original_trials), std(stats.Test_Intact_original_trials));
    fprintf('    Average trials removed per subject: %.1f (±%.1f)\n', mean(stats.Test_Intact_removals), std(stats.Test_Intact_removals));
    fprintf('    Average final trials per subject: %.1f (±%.1f)\n', mean(stats.Test_Intact_final_trials), std(stats.Test_Intact_final_trials));
    fprintf('    Group removal rate: %.1f%%\n', (sum(stats.Test_Intact_removals)/sum(stats.Test_Intact_original_trials))*100);
    
    fprintf('  Test_Recombined Group:\n');
    fprintf('    Average original trials per subject: %.1f (±%.1f)\n', mean(stats.Test_Recombined_original_trials), std(stats.Test_Recombined_original_trials));
    fprintf('    Average trials removed per subject: %.1f (±%.1f)\n', mean(stats.Test_Recombined_removals), std(stats.Test_Recombined_removals));
    fprintf('    Average final trials per subject: %.1f (±%.1f)\n', mean(stats.Test_Recombined_final_trials), std(stats.Test_Recombined_final_trials));
    fprintf('    Group removal rate: %.1f%%\n', (sum(stats.Test_Recombined_removals)/sum(stats.Test_Recombined_original_trials))*100);
    
    % NEW: Experimental Condition-Based Trial Statistics
    fprintf('\nExperimental Condition-Based Trial Statistics (averaged across %d processed subjects):\n', stats.processed_subjects);
    
    fprintf('  Study_hits:\n');
    fprintf('    Average original trials per subject: %.1f (±%.1f)\n', mean(stats.condition_original_trials.Study_hits), std(stats.condition_original_trials.Study_hits));
    fprintf('    Average trials removed per subject: %.1f (±%.1f)\n', mean(stats.condition_removed_trials.Study_hits), std(stats.condition_removed_trials.Study_hits));
    fprintf('    Average final trials per subject: %.1f (±%.1f)\n', mean(stats.condition_final_trials.Study_hits), std(stats.condition_final_trials.Study_hits));
    fprintf('    Condition removal rate: %.1f%%\n', (sum(stats.condition_removed_trials.Study_hits)/sum(stats.condition_original_trials.Study_hits))*100);
    
    fprintf('  Study_misses:\n');
    fprintf('    Average original trials per subject: %.1f (±%.1f)\n', mean(stats.condition_original_trials.Study_misses), std(stats.condition_original_trials.Study_misses));
    fprintf('    Average trials removed per subject: %.1f (±%.1f)\n', mean(stats.condition_removed_trials.Study_misses), std(stats.condition_removed_trials.Study_misses));
    fprintf('    Average final trials per subject: %.1f (±%.1f)\n', mean(stats.condition_final_trials.Study_misses), std(stats.condition_final_trials.Study_misses));
    fprintf('    Condition removal rate: %.1f%%\n', (sum(stats.condition_removed_trials.Study_misses)/sum(stats.condition_original_trials.Study_misses))*100);
    
    fprintf('  Test_hits:\n');
    fprintf('    Average original trials per subject: %.1f (±%.1f)\n', mean(stats.condition_original_trials.Test_hits), std(stats.condition_original_trials.Test_hits));
    fprintf('    Average trials removed per subject: %.1f (±%.1f)\n', mean(stats.condition_removed_trials.Test_hits), std(stats.condition_removed_trials.Test_hits));
    fprintf('    Average final trials per subject: %.1f (±%.1f)\n', mean(stats.condition_final_trials.Test_hits), std(stats.condition_final_trials.Test_hits));
    fprintf('    Condition removal rate: %.1f%%\n', (sum(stats.condition_removed_trials.Test_hits)/sum(stats.condition_original_trials.Test_hits))*100);
    
    fprintf('  Test_misses:\n');
    fprintf('    Average original trials per subject: %.1f (±%.1f)\n', mean(stats.condition_original_trials.Test_misses), std(stats.condition_original_trials.Test_misses));
    fprintf('    Average trials removed per subject: %.1f (±%.1f)\n', mean(stats.condition_removed_trials.Test_misses), std(stats.condition_removed_trials.Test_misses));
    fprintf('    Average final trials per subject: %.1f (±%.1f)\n', mean(stats.condition_final_trials.Test_misses), std(stats.condition_final_trials.Test_misses));
    fprintf('    Condition removal rate: %.1f%%\n', (sum(stats.condition_removed_trials.Test_misses)/sum(stats.condition_original_trials.Test_misses))*100);
    
    fprintf('  Correct_rejections:\n');
    fprintf('    Average original trials per subject: %.1f (±%.1f)\n', mean(stats.condition_original_trials.Correct_rejections), std(stats.condition_original_trials.Correct_rejections));
    fprintf('    Average trials removed per subject: %.1f (±%.1f)\n', mean(stats.condition_removed_trials.Correct_rejections), std(stats.condition_removed_trials.Correct_rejections));
    fprintf('    Average final trials per subject: %.1f (±%.1f)\n', mean(stats.condition_final_trials.Correct_rejections), std(stats.condition_final_trials.Correct_rejections));
    fprintf('    Condition removal rate: %.1f%%\n', (sum(stats.condition_removed_trials.Correct_rejections)/sum(stats.condition_original_trials.Correct_rejections))*100);
    
    fprintf('  False_alarms:\n');
    fprintf('    Average original trials per subject: %.1f (±%.1f)\n', mean(stats.condition_original_trials.False_alarms), std(stats.condition_original_trials.False_alarms));
    fprintf('    Average trials removed per subject: %.1f (±%.1f)\n', mean(stats.condition_removed_trials.False_alarms), std(stats.condition_removed_trials.False_alarms));
    fprintf('    Average final trials per subject: %.1f (±%.1f)\n', mean(stats.condition_final_trials.False_alarms), std(stats.condition_final_trials.False_alarms));
    fprintf('    Condition removal rate: %.1f%%\n', (sum(stats.condition_removed_trials.False_alarms)/sum(stats.condition_original_trials.False_alarms))*100);
    
    % NEW: Subject-Based Trial Statistics
    fprintf('\nSubject-Based Trial Statistics (individual subject performance):\n');
    fprintf('    Average original trials per subject: %.1f (±%.1f)\n', mean(stats.subject_original_trials), std(stats.subject_original_trials));
    fprintf('    Average trials removed per subject: %.1f (±%.1f)\n', mean(stats.subject_removed_trials), std(stats.subject_removed_trials));
    fprintf('    Average final trials per subject: %.1f (±%.1f)\n', mean(stats.subject_final_trials), std(stats.subject_final_trials));
    fprintf('    Subject-level removal rate: %.1f%%\n', (sum(stats.subject_removed_trials)/sum(stats.subject_original_trials))*100);
    fprintf('    Range of original trials: %d - %d trials\n', min(stats.subject_original_trials), max(stats.subject_original_trials));
    fprintf('    Range of removed trials: %d - %d trials\n', min(stats.subject_removed_trials), max(stats.subject_removed_trials));
    fprintf('    Range of final trials: %d - %d trials\n', min(stats.subject_final_trials), max(stats.subject_final_trials));
    
    % Per-subject removal distribution
    fprintf('\nPer-Subject Total Removal Distribution:\n');
    fprintf('  Average total trials removed per subject: %.1f (±%.1f)\n', mean(stats.total_removals_per_subject), std(stats.total_removals_per_subject));
    fprintf('  Range: %d - %d trials\n', min(stats.total_removals_per_subject), max(stats.total_removals_per_subject));
    
    % Subjects approaching exclusion threshold
    approaching_threshold = sum(stats.total_removals_per_subject > 40); % 40 as a warning threshold
    if approaching_threshold > 0
        fprintf('  WARNING: %d subjects had >40 total trials removed (approaching exclusion risk)\n', approaching_threshold);
    end
    
    % Quality assessment
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

% ==================== EXPORT STATISTICS TO TEXT FILE ====================
if stats.processed_subjects > 0
    fprintf('\n==================== EXPORTING STATISTICS ====================\n');
    
    % Create a comprehensive summary structure for export
    processing_summary = struct();
    
    % Basic processing info
    processing_summary.metadata = struct();
    processing_summary.metadata.processing_date = datestr(now);
    processing_summary.metadata.total_subjects_selected = stats.total_subjects;
    processing_summary.metadata.successfully_processed = stats.processed_subjects;
    processing_summary.metadata.excluded_subjects = stats.excluded_subjects;
    processing_summary.metadata.excluded_subject_ids = stats.excluded_subject_ids;
    processing_summary.metadata.exclusion_reasons = stats.exclusion_reasons;
    processing_summary.metadata.processing_success_rate = (stats.processed_subjects/stats.total_subjects)*100;
    
    % Processing parameters
    processing_summary.parameters = struct();
    processing_summary.parameters.voltage_diff_threshold = stats.voltage_diff_threshold;
    processing_summary.parameters.voltage_abs_threshold = stats.voltage_abs_threshold;
    processing_summary.parameters.max_rejected_trials_per_group = stats.max_rejected_trials_per_group;
    processing_summary.parameters.epoch_window = [-0.1, 1.5];
    processing_summary.parameters.baseline_window = [-100, 0];
    
    % Overall trial statistics
    processing_summary.overall = struct();
    processing_summary.overall.total_original_trials = stats.total_original_trials;
    processing_summary.overall.total_final_trials = stats.total_final_trials;
    processing_summary.overall.total_removed_trials = stats.total_removed_trials;
    processing_summary.overall.overall_removal_rate = (stats.total_removed_trials/stats.total_original_trials)*100;
    processing_summary.overall.data_retention_rate = (stats.total_final_trials / stats.total_original_trials) * 100;
    
    % Condition-based removal statistics
    processing_summary.condition_based = struct();
    processing_summary.condition_based.condition1_removals = stats.condition1_removals;
    processing_summary.condition_based.condition2_removals = stats.condition2_removals;
    processing_summary.condition_based.condition1_mean = mean(stats.condition1_removals);
    processing_summary.condition_based.condition1_std = std(stats.condition1_removals);
    processing_summary.condition_based.condition2_mean = mean(stats.condition2_removals);
    processing_summary.condition_based.condition2_std = std(stats.condition2_removals);
    processing_summary.condition_based.overlap_trials = sum(stats.condition1_removals) + sum(stats.condition2_removals) - stats.total_removed_trials;
    
    % Group-based statistics
    processing_summary.group_based = struct();
    processing_summary.group_based.SME = struct();
    processing_summary.group_based.SME.original_trials = stats.SME_original_trials;
    processing_summary.group_based.SME.removed_trials = stats.SME_removals;
    processing_summary.group_based.SME.final_trials = stats.SME_final_trials;
    processing_summary.group_based.SME.mean_original = mean(stats.SME_original_trials);
    processing_summary.group_based.SME.mean_removed = mean(stats.SME_removals);
    processing_summary.group_based.SME.mean_final = mean(stats.SME_final_trials);
    processing_summary.group_based.SME.std_original = std(stats.SME_original_trials);
    processing_summary.group_based.SME.std_removed = std(stats.SME_removals);
    processing_summary.group_based.SME.std_final = std(stats.SME_final_trials);
    processing_summary.group_based.SME.removal_rate = (sum(stats.SME_removals)/sum(stats.SME_original_trials))*100;
    
    processing_summary.group_based.Test_Intact = struct();
    processing_summary.group_based.Test_Intact.original_trials = stats.Test_Intact_original_trials;
    processing_summary.group_based.Test_Intact.removed_trials = stats.Test_Intact_removals;
    processing_summary.group_based.Test_Intact.final_trials = stats.Test_Intact_final_trials;
    processing_summary.group_based.Test_Intact.mean_original = mean(stats.Test_Intact_original_trials);
    processing_summary.group_based.Test_Intact.mean_removed = mean(stats.Test_Intact_removals);
    processing_summary.group_based.Test_Intact.mean_final = mean(stats.Test_Intact_final_trials);
    processing_summary.group_based.Test_Intact.std_original = std(stats.Test_Intact_original_trials);
    processing_summary.group_based.Test_Intact.std_removed = std(stats.Test_Intact_removals);
    processing_summary.group_based.Test_Intact.std_final = std(stats.Test_Intact_final_trials);
    processing_summary.group_based.Test_Intact.removal_rate = (sum(stats.Test_Intact_removals)/sum(stats.Test_Intact_original_trials))*100;
    
    processing_summary.group_based.Test_Recombined = struct();
    processing_summary.group_based.Test_Recombined.original_trials = stats.Test_Recombined_original_trials;
    processing_summary.group_based.Test_Recombined.removed_trials = stats.Test_Recombined_removals;
    processing_summary.group_based.Test_Recombined.final_trials = stats.Test_Recombined_final_trials;
    processing_summary.group_based.Test_Recombined.mean_original = mean(stats.Test_Recombined_original_trials);
    processing_summary.group_based.Test_Recombined.mean_removed = mean(stats.Test_Recombined_removals);
    processing_summary.group_based.Test_Recombined.mean_final = mean(stats.Test_Recombined_final_trials);
    processing_summary.group_based.Test_Recombined.std_original = std(stats.Test_Recombined_original_trials);
    processing_summary.group_based.Test_Recombined.std_removed = std(stats.Test_Recombined_removals);
    processing_summary.group_based.Test_Recombined.std_final = std(stats.Test_Recombined_final_trials);
    processing_summary.group_based.Test_Recombined.removal_rate = (sum(stats.Test_Recombined_removals)/sum(stats.Test_Recombined_original_trials))*100;
    
    % Experimental condition-based statistics
    processing_summary.experimental_conditions = struct();
    for cond = condition_names
        cond_name = cond{1};
        processing_summary.experimental_conditions.(cond_name) = struct();
        processing_summary.experimental_conditions.(cond_name).original_trials = stats.condition_original_trials.(cond_name);
        processing_summary.experimental_conditions.(cond_name).removed_trials = stats.condition_removed_trials.(cond_name);
        processing_summary.experimental_conditions.(cond_name).final_trials = stats.condition_final_trials.(cond_name);
        processing_summary.experimental_conditions.(cond_name).mean_original = mean(stats.condition_original_trials.(cond_name));
        processing_summary.experimental_conditions.(cond_name).mean_removed = mean(stats.condition_removed_trials.(cond_name));
        processing_summary.experimental_conditions.(cond_name).mean_final = mean(stats.condition_final_trials.(cond_name));
        processing_summary.experimental_conditions.(cond_name).std_original = std(stats.condition_original_trials.(cond_name));
        processing_summary.experimental_conditions.(cond_name).std_removed = std(stats.condition_removed_trials.(cond_name));
        processing_summary.experimental_conditions.(cond_name).std_final = std(stats.condition_final_trials.(cond_name));
        if sum(stats.condition_original_trials.(cond_name)) > 0
            processing_summary.experimental_conditions.(cond_name).removal_rate = (sum(stats.condition_removed_trials.(cond_name))/sum(stats.condition_original_trials.(cond_name)))*100;
        else
            processing_summary.experimental_conditions.(cond_name).removal_rate = 0;
        end
    end
    
    % Subject-based statistics
    processing_summary.subject_based = struct();
    processing_summary.subject_based.subject_ids = stats.subject_ids;
    processing_summary.subject_based.original_trials = stats.subject_original_trials;
    processing_summary.subject_based.removed_trials = stats.subject_removed_trials;
    processing_summary.subject_based.final_trials = stats.subject_final_trials;
    processing_summary.subject_based.mean_original = mean(stats.subject_original_trials);
    processing_summary.subject_based.mean_removed = mean(stats.subject_removed_trials);
    processing_summary.subject_based.mean_final = mean(stats.subject_final_trials);
    processing_summary.subject_based.std_original = std(stats.subject_original_trials);
    processing_summary.subject_based.std_removed = std(stats.subject_removed_trials);
    processing_summary.subject_based.std_final = std(stats.subject_final_trials);
    processing_summary.subject_based.removal_rate = (sum(stats.subject_removed_trials)/sum(stats.subject_original_trials))*100;
    processing_summary.subject_based.range_original = [min(stats.subject_original_trials), max(stats.subject_original_trials)];
    processing_summary.subject_based.range_removed = [min(stats.subject_removed_trials), max(stats.subject_removed_trials)];
    processing_summary.subject_based.range_final = [min(stats.subject_final_trials), max(stats.subject_final_trials)];
    
    % Quality assessment
    processing_summary.quality_assessment = struct();
    processing_summary.quality_assessment.retention_rate = (stats.total_final_trials / stats.total_original_trials) * 100;
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
    
    % Save the summary to a text file
    summary_filename = sprintf('Processing_Summary_%s.txt', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
    summary_filepath = fullfile(epoched_dir, summary_filename);
    
    try
        write_summary_text_file(summary_filepath, processing_summary, stats);
        fprintf('Statistics summary exported to: %s\n', summary_filepath);
        fprintf('Contains: processing_summary and stats as flattened key/value text.\n');
    catch ME
        fprintf('Error saving statistics file: %s\n', ME.message);
        fprintf('Statistics could not be exported, but processing completed successfully.\n');
    end
else
    fprintf('\nNo statistics to export - no subjects were successfully processed.\n');
end

fprintf('\n==================== PROCESSING COMPLETE ====================\n');


%% ================== HELPER FUNCTIONS (copied from original script) ==================
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
    clean_EEG = EEG; % Return original EEG unchanged
end

function trial_indices_by_group = assign_epochs_to_groups_fixed(EEG)
    % Helper function to assign epochs to artifact-rejection groups.
    % Primary path: use derived trial_type labels.
    % Fallback path: infer from trigger type + testcresp.

    trial_indices_by_group = struct();
    trial_indices_by_group.SME = [];
    trial_indices_by_group.Test_Intact = [];
    trial_indices_by_group.Test_Recombined = [];

    epoch_trial_type = repmat({''}, 1, EEG.trials);
    has_trial_type = false;

    if isfield(EEG.event, 'trial_type')
        has_trial_type = true;
        for ev_idx = 1:length(EEG.event)
            if ~isfield(EEG.event(ev_idx), 'epoch')
                continue;
            end
            ep = EEG.event(ev_idx).epoch;
            if ep < 1 || ep > EEG.trials
                continue;
            end
            tt = EEG.event(ev_idx).trial_type;
            if isstring(tt)
                tt = char(tt);
            end
            if isempty(tt)
                continue;
            end
            if isempty(epoch_trial_type{ep})
                epoch_trial_type{ep} = tt;
            end
        end
    end

    if has_trial_type && any(~cellfun(@isempty, epoch_trial_type))
        trial_indices_by_group.SME = find(ismember(epoch_trial_type, {'Study_hits', 'Study_misses'}));
        trial_indices_by_group.Test_Intact = find(ismember(epoch_trial_type, {'Test_hits', 'Test_misses'}));
        trial_indices_by_group.Test_Recombined = find(ismember(epoch_trial_type, {'Correct_rejections', 'False_alarms'}));
        return;
    end

    [epoch_trigger_type, ~, ~, epoch_testcresp] = summarize_epoch_trigger_info_item(EEG);
    trial_indices_by_group.SME = find(strcmp(epoch_trigger_type, '2'));
    trial_indices_by_group.Test_Intact = find(strcmp(epoch_trigger_type, '7') & epoch_testcresp == 1);
    trial_indices_by_group.Test_Recombined = find(strcmp(epoch_trigger_type, '7') & epoch_testcresp == 2);
end

function [epoch_trigger_type, epoch_trigger_latency, epoch_testacc, epoch_testcresp] = summarize_epoch_trigger_info_item(EEG)
    % Summarize, for each epoch, the earliest trigger event and associated
    % test fields used for old/new condition classification.
    epoch_trigger_type = repmat({''}, 1, EEG.trials);
    epoch_trigger_latency = inf(1, EEG.trials);
    epoch_testacc = nan(1, EEG.trials);
    epoch_testcresp = nan(1, EEG.trials);

    for ev_idx = 1:length(EEG.event)
        if ~isfield(EEG.event(ev_idx), 'epoch')
            continue;
        end
        ep = EEG.event(ev_idx).epoch;
        if ep < 1 || ep > EEG.trials
            continue;
        end

        if ~isfield(EEG.event(ev_idx), 'type')
            continue;
        end
        event_type_str = event_type_as_char(EEG.event(ev_idx).type);

        if ~ismember(event_type_str, {'2', '7'})
            continue;
        end

        if isfield(EEG.event(ev_idx), 'latency')
            this_latency = double(EEG.event(ev_idx).latency);
        else
            this_latency = inf;
        end

        if this_latency < epoch_trigger_latency(ep)
            epoch_trigger_latency(ep) = this_latency;
            epoch_trigger_type{ep} = event_type_str;

            if isfield(EEG.event(ev_idx), 'testacc') && isnumeric(EEG.event(ev_idx).testacc) && isscalar(EEG.event(ev_idx).testacc)
                epoch_testacc(ep) = double(EEG.event(ev_idx).testacc);
            else
                epoch_testacc(ep) = nan;
            end

            if isfield(EEG.event(ev_idx), 'testcresp') && isnumeric(EEG.event(ev_idx).testcresp) && isscalar(EEG.event(ev_idx).testcresp)
                epoch_testcresp(ep) = double(EEG.event(ev_idx).testcresp);
            else
                epoch_testcresp(ep) = nan;
            end
        end
    end
end

function cond_name = classify_test_condition(testacc, testcresp)
    % Test-phase condition mapping:
    %   Hits:               testacc==1 && testcresp==1
    %   Misses:             testacc==0 && testcresp==1
    %   Correct rejections: testacc==1 && testcresp==2
    %   False alarms:       testacc==0 && testcresp==2
    cond_name = '';

    if isnan(testacc) || isnan(testcresp)
        return;
    end

    if testcresp == 1 && testacc == 1
        cond_name = 'Test_hits';
    elseif testcresp == 1 && testacc == 0
        cond_name = 'Test_misses';
    elseif testcresp == 2 && testacc == 1
        cond_name = 'Correct_rejections';
    elseif testcresp == 2 && testacc == 0
        cond_name = 'False_alarms';
    end
end

function out = event_type_as_char(event_type)
    if isnumeric(event_type)
        out = num2str(event_type);
    elseif isstring(event_type)
        out = char(event_type);
    else
        out = char(event_type);
    end
    out = strtrim(out);
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
