% Calculate and plot mean voltage across all participants in an EEGLAB STUDY
% with bad trial detection and removal
% Configuration
hardcoded_event_field = 'trial_type';
hardcoded_event_values = {'Shits', 'Smiss'};
hardcoded_channels = [21];

% Bad trial detection parameters (same as original code)
voltage_diff_threshold = 25;    % μV - maximum allowed difference between consecutive samples
voltage_abs_threshold = 1000;   % μV - maximum allowed absolute voltage
max_rejected_trials = 56;      % maximum trials that can be rejected per subject

% Check if STUDY and ALLEEG exist
if ~exist('STUDY', 'var') || ~exist('ALLEEG', 'var')
    error('STUDY and ALLEEG must be loaded in workspace');
end

% Get number of datasets
n_datasets = length(STUDY.datasetinfo);

% Load first dataset to get structure info
dataset_idx = STUDY.datasetinfo(1).index;
filepath = ALLEEG(dataset_idx).filepath;
filename = ALLEEG(dataset_idx).filename;
EEG_first = pop_loadset('filename', filename, 'filepath', filepath);

% Initialize storage
n_selected_values = length(hardcoded_event_values);
all_data_by_value = cell(n_selected_values, 1);

% Initialize bad trial tracking
bad_trial_summary = struct();
subjects_with_excessive_rejection = [];
total_trials_rejected_condition1 = 0;
total_trials_rejected_condition2 = 0;
total_trials_rejected_both = 0;

fprintf('=== BAD TRIAL DETECTION SUMMARY ===\n');
fprintf('Condition 1: Voltage difference > %.1f μV\n', voltage_diff_threshold);
fprintf('Condition 2: Absolute voltage > %.1f μV\n', voltage_abs_threshold);
fprintf('=========================================\n\n');

% Process each dataset
for i = 1:n_datasets
    dataset_idx = STUDY.datasetinfo(i).index;
    filepath = ALLEEG(dataset_idx).filepath;
    filename = ALLEEG(dataset_idx).filename;
    EEG = pop_loadset('filename', filename, 'filepath', filepath);
    
    % Get subject ID from filename or dataset info
    if isfield(STUDY.datasetinfo(i), 'subject')
        subject_id = STUDY.datasetinfo(i).subject;
    else
        [~, subject_id, ~] = fileparts(filename);
    end
    
    % Initialize storage on first iteration
    if i == 1
        n_timepoints = size(EEG.data, 2);
        n_channels = length(hardcoded_channels);
        for v = 1:n_selected_values
            all_data_by_value{v} = zeros(n_channels, n_timepoints, n_datasets);
        end
    end
    
    % Initialize bad trial tracking for this subject
    bad_trial_summary.(sprintf('subject_%d', i)) = struct();
    bad_trial_summary.(sprintf('subject_%d', i)).subject_id = subject_id;
    bad_trial_summary.(sprintf('subject_%d', i)).condition1_trials = [];
    bad_trial_summary.(sprintf('subject_%d', i)).condition2_trials = [];
    bad_trial_summary.(sprintf('subject_%d', i)).total_bad_trials = [];
    
    % Process each event value
    for v = 1:n_selected_values
        current_value = hardcoded_event_values{v};
        
        % Find matching events
        if isfield(EEG, 'event') && ~isempty(EEG.event)
            event_field_values = {EEG.event.(hardcoded_event_field)};
            matching_events = find(strcmp(event_field_values, current_value));
            
            if ~isempty(matching_events) && ndims(EEG.data) == 3
                % Get epochs for matching events
                matching_epochs = [EEG.event(matching_events).epoch];
                matching_epochs = unique(matching_epochs);
                
                if ~isempty(matching_epochs)
                    epoch_data = EEG.data(:, :, matching_epochs);
                    
                    % BAD TRIAL DETECTION
                    bad_trials_condition1 = [];
                    bad_trials_condition2 = [];
                    
                    % Check each epoch for bad trials
                    for epoch_idx = 1:length(matching_epochs)
                        epoch_num = matching_epochs(epoch_idx);
                        
                        % Check specified channels for artifacts
                        for ch = hardcoded_channels
                            % Get voltage data for this channel and epoch
                            voltage_data = squeeze(EEG.data(ch, :, epoch_num));
                            
                            % CONDITION 1: Check voltage differences between consecutive samples
                            voltage_diffs = abs(diff(voltage_data));
                            max_voltage_diff = max(voltage_diffs);
                            if max_voltage_diff > voltage_diff_threshold
                                bad_trials_condition1 = [bad_trials_condition1, epoch_num];
                            end
                            
                            % CONDITION 2: Check absolute voltage threshold
                            max_abs_voltage = max(abs(voltage_data));
                            if max_abs_voltage > voltage_abs_threshold
                                bad_trials_condition2 = [bad_trials_condition2, epoch_num];
                            end
                        end
                    end
                    
                    % Get unique bad trials
                    bad_trials_condition1 = unique(bad_trials_condition1);
                    bad_trials_condition2 = unique(bad_trials_condition2);
                    all_bad_trials = unique([bad_trials_condition1, bad_trials_condition2]);
                    
                    % Store bad trial info
                    bad_trial_summary.(sprintf('subject_%d', i)).condition1_trials = ...
                        [bad_trial_summary.(sprintf('subject_%d', i)).condition1_trials, bad_trials_condition1];
                    bad_trial_summary.(sprintf('subject_%d', i)).condition2_trials = ...
                        [bad_trial_summary.(sprintf('subject_%d', i)).condition2_trials, bad_trials_condition2];
                    bad_trial_summary.(sprintf('subject_%d', i)).total_bad_trials = ...
                        [bad_trial_summary.(sprintf('subject_%d', i)).total_bad_trials, all_bad_trials];
                    
                    % Remove bad trials from matching_epochs
                    clean_epochs = setdiff(matching_epochs, all_bad_trials);
                    
                    % Calculate participant mean from clean epochs
                    if ~isempty(clean_epochs)
                        clean_epoch_data = EEG.data(:, :, clean_epochs);
                        participant_mean = mean(clean_epoch_data, 3);
                    else
                        participant_mean = zeros(size(EEG.data, 1), size(EEG.data, 2));
                    end
                    
                    % Print subject-specific results
                    if ~isempty(all_bad_trials)
                        fprintf('Subject %s: %d bad trials removed (Condition 1: %d, Condition 2: %d)\n', ...
                            subject_id, length(all_bad_trials), length(bad_trials_condition1), length(bad_trials_condition2));
                        
                        % Update totals
                        total_trials_rejected_condition1 = total_trials_rejected_condition1 + length(bad_trials_condition1);
                        total_trials_rejected_condition2 = total_trials_rejected_condition2 + length(bad_trials_condition2);
                        total_trials_rejected_both = total_trials_rejected_both + length(all_bad_trials);
                        
                        % Check if too many trials rejected
                        if length(all_bad_trials) > max_rejected_trials
                            subjects_with_excessive_rejection = [subjects_with_excessive_rejection, i];
                            fprintf('  WARNING: Subject %s has too many rejected trials (%d > %d)\n', ...
                                subject_id, length(all_bad_trials), max_rejected_trials);
                        end
                    end
                    
                else
                    participant_mean = zeros(size(EEG.data, 1), size(EEG.data, 2));
                end
            else
                participant_mean = zeros(size(EEG.data, 1), size(EEG.data, 2));
            end
        else
            participant_mean = zeros(size(EEG.data, 1), size(EEG.data, 2));
        end
        
        % Select specified channels and store
        all_data_by_value{v}(:, :, i) = participant_mean(hardcoded_channels, :);
    end
end

% Print summary statistics
fprintf('\n=== SUMMARY STATISTICS ===\n');
fprintf('Total trials rejected by Condition 1 (voltage diff > %.1f μV): %d\n', voltage_diff_threshold, total_trials_rejected_condition1);
fprintf('Total trials rejected by Condition 2 (abs voltage > %.1f μV): %d\n', voltage_abs_threshold, total_trials_rejected_condition2);
fprintf('Total unique trials rejected (both conditions): %d\n', total_trials_rejected_both);
fprintf('Number of subjects with excessive trial rejection (>%d): %d\n', max_rejected_trials, length(subjects_with_excessive_rejection));

if ~isempty(subjects_with_excessive_rejection)
    fprintf('Subject IDs with excessive rejection: ');
    for idx = subjects_with_excessive_rejection
        subject_id = bad_trial_summary.(sprintf('subject_%d', idx)).subject_id;
        fprintf('%s ', subject_id);
    end
    fprintf('\n');
end

% Calculate grand average
mean_voltage_by_value = cell(n_selected_values, 1);
for v = 1:n_selected_values
    mean_voltage_by_value{v} = mean(all_data_by_value{v}, 3);
end

% Plot results
figure;
colors = lines(n_selected_values);

% Get time vector
if isfield(EEG, 'times') && length(EEG.times) == size(mean_voltage_by_value{1}, 2)
    time_vector = EEG.times;
    xlabel_text = 'Time (ms)';
else
    time_vector = 1:size(mean_voltage_by_value{1}, 2);
    xlabel_text = 'Time Points';
end

% Plot each condition
hold on;
legend_entries = {};
for v = 1:n_selected_values
    current_data = mean_voltage_by_value{v};
    % Average across channels if multiple
    if size(current_data, 1) > 1
        plot_data = mean(current_data, 1);
    else
        plot_data = current_data;
    end
    plot(time_vector, plot_data, 'Color', colors(v, :), 'LineWidth', 2);
    legend_entries{end+1} = hardcoded_event_values{v};
end

% Format plot
xlabel(xlabel_text);
ylabel('Voltage (µV)');
title('ERP Plot (After Bad Trial Removal)');
legend(legend_entries, 'Location', 'best');
grid on;
hold off;

fprintf('\n=== ANALYSIS COMPLETE ===\n');
fprintf('Clean data plotted with bad trials removed.\n');