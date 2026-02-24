function results = hpc_epoch_to_erp_plot(input_folder, trial_type_values, channels, output_dir, figure_title)
% HPC_EPOCH_TO_ERP_PLOT - Plot grand-average ERP curves by trial_type.
%
% Usage:
%   hpc_epoch_to_erp_plot(input_folder, trial_type_values)
%   hpc_epoch_to_erp_plot(input_folder, trial_type_values, channels)
%   hpc_epoch_to_erp_plot(input_folder, trial_type_values, channels, output_dir)
%   hpc_epoch_to_erp_plot(input_folder, trial_type_values, channels, output_dir, figure_title)
%
% Inputs:
%   input_folder      - Path to epoch folder containing *_epoch.set files.
%   trial_type_values - Trial type labels to plot. Accepts CSV/semicolon list, string array,
%                       or cell array (e.g., {'Study_hits','Study_misses'}).
%   channels          - Optional channel indices (default: 21).
%                       Accepts numeric vector or CSV string (e.g., '21,22').
%   output_dir        - Optional output folder (default: <input_folder>/erp_plots).
%   figure_title      - Optional custom figure title.
%
% Outputs:
%   - <output_dir>/grand_average_erp_<trial_types>.png
%   - <output_dir>/grand_average_erp_<trial_types>.svg (best effort on headless HPC)
%   - <output_dir>/grand_average_erp_<trial_types>.mat

if nargin < 2 || isempty(trial_type_values)
    error('trial_type_values is required (CSV, string array, or cell array).');
end
if nargin < 3 || isempty(channels)
    channels = 21;
end
if nargin < 4 || isempty(output_dir)
    output_dir = fullfile(input_folder, 'erp_plots');
end
if nargin < 5 || isempty(figure_title)
    figure_title = 'Grand Average ERP by trial_type';
end

trial_type_values = parse_trial_type_values(trial_type_values);
channels = parse_channel_list(channels);

fprintf('\n==============================================\n');
fprintf('  HPC EPOCH -> ERP GRAND AVERAGE PLOT\n');
fprintf('==============================================\n');
fprintf('Input folder:   %s\n', input_folder);
fprintf('Trial types:    %s\n', strjoin(trial_type_values, ', '));
fprintf('Channels:       %s\n', num2str(channels));
fprintf('Output folder:  %s\n', output_dir);
fprintf('==============================================\n\n');

if ~exist(input_folder, 'dir')
    error('Input folder does not exist: %s', input_folder);
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

set(0, 'DefaultFigureVisible', 'off');
java.lang.System.setProperty('java.awt.headless', 'true');

addpath('/home/devon7y/scratch/devon7y/eeglab2022.1');
eeglab nogui;

epoch_files = dir(fullfile(input_folder, '*_epoch.set'));
if isempty(epoch_files)
    error('No *_epoch.set files found in: %s', input_folder);
end

n_subjects = numel(epoch_files);
n_conditions = numel(trial_type_values);

subject_ids = cell(n_subjects, 1);
trial_counts = zeros(n_subjects, n_conditions);

all_data_by_condition = cell(n_conditions, 1);
time_vector = [];
xlabel_text = 'Time (ms)';
channels_in_use = [];
n_timepoints = [];

fprintf('Processing %d epoched dataset(s)\n', n_subjects);

for i = 1:n_subjects
    fname = epoch_files(i).name;
    fpath = epoch_files(i).folder;
    subject_ids{i} = parse_subject_label(fname);

    fprintf('\n----------------------------------------------\n');
    fprintf('Subject %s (%d/%d): %s\n', subject_ids{i}, i, n_subjects, fname);

    EEG = pop_loadset('filename', fname, 'filepath', fpath);

    if i == 1
        channels_in_use = validate_channels(channels, EEG.nbchan);
        n_timepoints = size(EEG.data, 2);
        [time_vector, xlabel_text] = infer_time_vector(EEG, n_timepoints);

        for c = 1:n_conditions
            all_data_by_condition{c} = nan(numel(channels_in_use), n_timepoints, n_subjects);
        end
    else
        this_n_time = size(EEG.data, 2);
        if this_n_time ~= n_timepoints
            error(['Time-point mismatch for %s: expected %d points, found %d. ' ...
                   'All epoched files must share epoch length.'], ...
                fname, n_timepoints, this_n_time);
        end
        if any(channels_in_use > EEG.nbchan)
            error(['Channel index out of range for %s: requested %s but dataset has %d channels.'], ...
                fname, num2str(channels_in_use), EEG.nbchan);
        end
    end

    epoch_trial_type = extract_epoch_trial_type_labels(EEG);

    for c = 1:n_conditions
        label = trial_type_values{c};
        matching_epochs = find(strcmp(epoch_trial_type, label));
        trial_counts(i, c) = numel(matching_epochs);

        fprintf('  %-30s trials: %d\n', label, trial_counts(i, c));

        if ~isempty(matching_epochs)
            participant_mean = mean(EEG.data(channels_in_use, :, matching_epochs), 3);
            all_data_by_condition{c}(:, :, i) = participant_mean;
        end
    end
end

grand_average_by_condition = cell(n_conditions, 1);
participants_per_condition = zeros(1, n_conditions);
for c = 1:n_conditions
    participants_per_condition(c) = sum(trial_counts(:, c) > 0);
    grand_average_by_condition{c} = mean(all_data_by_condition{c}, 3, 'omitnan');
end

fig = figure('Color', 'w', 'Position', [100 100 1200 700], 'Visible', 'off');
hold on;
colors = lines(n_conditions);
legend_entries = {};
plotted_any = false;

for c = 1:n_conditions
    if participants_per_condition(c) < 1
        fprintf('WARNING: No matching trials found for "%s" across all subjects.\n', trial_type_values{c});
        continue;
    end

    curve_data = grand_average_by_condition{c};
    if size(curve_data, 1) > 1
        curve_data = mean(curve_data, 1, 'omitnan');
    end

    plot(time_vector, curve_data, 'Color', colors(c, :), 'LineWidth', 2.2);
    legend_entries{end+1} = sprintf('%s (n=%d)', trial_type_values{c}, participants_per_condition(c)); %#ok<AGROW>
    plotted_any = true;
end

if ~plotted_any
    close(fig);
    error('None of the requested trial_type labels were found in the epoch files.');
end

if strcmp(xlabel_text, 'Time (ms)')
    xline(0, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.0);
end
yline(0, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);
grid on;
box off;

xlabel(xlabel_text);
ylabel('Voltage (\muV)');
title(figure_title, 'Interpreter', 'none');
legend(legend_entries, 'Location', 'best', 'Interpreter', 'none');
hold off;

trial_type_tag = build_trial_type_tag(trial_type_values);
output_png = fullfile(output_dir, sprintf('grand_average_erp_%s.png', trial_type_tag));
output_svg = fullfile(output_dir, sprintf('grand_average_erp_%s.svg', trial_type_tag));
output_mat = fullfile(output_dir, sprintf('grand_average_erp_%s.mat', trial_type_tag));

exportgraphics(fig, output_png, 'Resolution', 300);
svg_saved = save_svg_with_fallback(fig, output_svg);
close(fig);

save(output_mat, 'trial_type_values', 'channels_in_use', 'time_vector', ...
    'grand_average_by_condition', 'trial_counts', 'participants_per_condition', ...
    'subject_ids');

fprintf('\n==============================================\n');
fprintf('  ERP PLOT COMPLETE\n');
fprintf('==============================================\n');
for c = 1:n_conditions
    fprintf('  %-30s participants with data: %d\n', ...
        trial_type_values{c}, participants_per_condition(c));
end
fprintf('Saved PNG: %s\n', output_png);
if svg_saved
    fprintf('Saved SVG: %s\n', output_svg);
else
    fprintf('SVG export unavailable in this MATLAB/headless config; skipped SVG.\n');
end
fprintf('Saved MAT: %s\n', output_mat);
fprintf('==============================================\n\n');

if nargout > 0
    results = struct();
    results.trial_type_values = trial_type_values;
    results.channels = channels_in_use;
    results.time_vector = time_vector;
    results.grand_average_by_condition = grand_average_by_condition;
    results.trial_counts = trial_counts;
    results.participants_per_condition = participants_per_condition;
    results.subject_ids = subject_ids;
    results.output_png = output_png;
    if svg_saved
        results.output_svg = output_svg;
    else
        results.output_svg = '';
    end
    results.output_mat = output_mat;
end

end

function labels = parse_trial_type_values(value)
if ischar(value) || (isstring(value) && isscalar(value))
    raw = strtrim(char(value));
    if contains(raw, ';')
        parts = strsplit(raw, ';');
    else
        parts = strsplit(raw, ',');
    end
    labels = cellfun(@strtrim, parts, 'UniformOutput', false);
elseif isstring(value)
    labels = arrayfun(@(x) strtrim(char(x)), value(:)', 'UniformOutput', false);
elseif iscell(value)
    labels = cell(1, numel(value));
    for i = 1:numel(value)
        labels{i} = normalize_label_value(value{i});
    end
else
    error('trial_type_values must be CSV, string array, or cell array.');
end

labels = labels(~cellfun(@isempty, labels));
labels = unique(labels, 'stable');
if isempty(labels)
    error('No valid trial_type labels provided.');
end
end

function ok = save_svg_with_fallback(fig, output_svg)
ok = false;

try
    exportgraphics(fig, output_svg, 'ContentType', 'vector');
    ok = true;
    return;
catch ME
    warning('exportgraphics SVG failed: %s', ME.message);
end

try
    print(fig, output_svg, '-dsvg');
    ok = true;
    return;
catch ME
    warning('print -dsvg failed: %s', ME.message);
end

try
    saveas(fig, output_svg);
    ok = exist(output_svg, 'file') == 2;
catch ME
    warning('saveas SVG failed: %s', ME.message);
end
end

function out = parse_channel_list(value)
if isnumeric(value)
    out = value(:)';
elseif ischar(value) || (isstring(value) && isscalar(value))
    parts = strsplit(char(value), ',');
    out = nan(1, numel(parts));
    for i = 1:numel(parts)
        out(i) = str2double(strtrim(parts{i}));
    end
elseif isstring(value)
    out = nan(1, numel(value));
    for i = 1:numel(value)
        out(i) = str2double(strtrim(char(value(i))));
    end
elseif iscell(value)
    out = nan(1, numel(value));
    for i = 1:numel(value)
        out(i) = str2double(strtrim(char(value{i})));
    end
else
    error('channels must be numeric vector or CSV string.');
end

if isempty(out)
    error('No channel indices provided.');
end
if any(~isfinite(out))
    error('channels contains non-numeric values.');
end
if any(mod(out, 1) ~= 0) || any(out < 1)
    error('channels must contain positive integer indices.');
end

out = unique(double(out), 'stable');
end

function channels_out = validate_channels(channels_in, nbchan)
channels_out = channels_in(channels_in >= 1 & channels_in <= nbchan);
if isempty(channels_out)
    error('No requested channels are valid for this dataset (nbchan=%d).', nbchan);
end
if numel(channels_out) < numel(channels_in)
    invalid = setdiff(channels_in, channels_out);
    error('Invalid channel indices for nbchan=%d: %s', nbchan, num2str(invalid));
end
end

function [time_vector, xlabel_text] = infer_time_vector(EEG, n_timepoints)
xlabel_text = 'Time (ms)';
if isfield(EEG, 'times') && numel(EEG.times) == n_timepoints
    time_vector = double(EEG.times(:))';
else
    time_vector = 1:n_timepoints;
    xlabel_text = 'Time Points';
end
end

function labels = extract_epoch_trial_type_labels(EEG)
labels = repmat({''}, 1, EEG.trials);
best_latency = inf(1, EEG.trials);

if isfield(EEG, 'event') && ~isempty(EEG.event)
    for ev = 1:numel(EEG.event)
        if ~isfield(EEG.event(ev), 'epoch')
            continue;
        end

        epoch_idx = parse_epoch_index(EEG.event(ev).epoch);
        if ~isfinite(epoch_idx) || epoch_idx < 1 || epoch_idx > EEG.trials
            continue;
        end

        trial_label = '';
        if isfield(EEG.event(ev), 'trial_type')
            trial_label = normalize_label_value(EEG.event(ev).trial_type);
        end
        if isempty(trial_label)
            continue;
        end

        event_latency = inf;
        if isfield(EEG.event(ev), 'latency')
            event_latency = parse_latency_value(EEG.event(ev).latency);
        end

        if isempty(labels{epoch_idx}) || event_latency < best_latency(epoch_idx)
            labels{epoch_idx} = trial_label;
            best_latency(epoch_idx) = event_latency;
        end
    end
end

if isfield(EEG, 'epoch') && ~isempty(EEG.epoch)
    n_epoch_struct = min(EEG.trials, numel(EEG.epoch));
    for ep = 1:n_epoch_struct
        if isempty(labels{ep}) && isfield(EEG.epoch(ep), 'eventtrial_type')
            labels{ep} = normalize_label_value(EEG.epoch(ep).eventtrial_type);
        end
    end
end
end

function out = parse_epoch_index(value)
if isnumeric(value)
    if isempty(value)
        out = nan;
    else
        out = double(value(1));
    end
elseif iscell(value)
    if isempty(value)
        out = nan;
    else
        out = parse_epoch_index(value{1});
    end
else
    out = str2double(strtrim(char(value)));
end

if ~isfinite(out)
    out = nan;
else
    out = round(out);
end
end

function out = parse_latency_value(value)
if isnumeric(value)
    if isempty(value)
        out = inf;
    else
        out = double(value(1));
    end
elseif iscell(value)
    if isempty(value)
        out = inf;
    else
        out = parse_latency_value(value{1});
    end
else
    out = str2double(strtrim(char(value)));
    if ~isfinite(out)
        out = inf;
    end
end
end

function out = normalize_label_value(value)
if iscell(value)
    out = '';
    for i = 1:numel(value)
        out = normalize_label_value(value{i});
        if ~isempty(out)
            return;
        end
    end
elseif isstring(value)
    if isempty(value)
        out = '';
    else
        out = strtrim(char(value(1)));
    end
elseif ischar(value)
    out = strtrim(value);
elseif isnumeric(value)
    if isempty(value)
        out = '';
    else
        out = strtrim(num2str(value(1)));
    end
else
    out = '';
end
end

function subject_label = parse_subject_label(filename)
[~, basename, ~] = fileparts(filename);
tokens = regexp(basename, '^(\d+)_epoch$', 'tokens');
if ~isempty(tokens)
    subject_label = tokens{1}{1};
else
    subject_label = basename;
end
end

function out = build_trial_type_tag(labels)
parts = cell(1, numel(labels));
for i = 1:numel(labels)
    token = lower(labels{i});
    token = regexprep(token, '[^a-z0-9]+', '_');
    token = regexprep(token, '^_+|_+$', '');
    token = regexprep(token, '_+', '_');
    if isempty(token)
        token = sprintf('cond%d', i);
    end
    parts{i} = token;
end

out = strjoin(parts, '_vs_');
if isempty(out)
    out = 'selected_conditions';
end
if numel(out) > 120
    out = out(1:120);
end
end
