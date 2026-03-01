function results = hpc_epoch_to_erp_plot_single_subject(input_folder, trial_type_values, channels, output_dir, figure_title, time_window_ms, condition_order_source)
% HPC_EPOCH_TO_ERP_PLOT_SINGLE_SUBJECT - Plot grand-average ERP curves from
% LIMO single-subject parameter files.
%
% Usage:
%   hpc_epoch_to_erp_plot_single_subject(input_folder, trial_type_values, channels)
%   hpc_epoch_to_erp_plot_single_subject(input_folder, trial_type_values, channels, output_dir)
%   hpc_epoch_to_erp_plot_single_subject(input_folder, trial_type_values, channels, output_dir, figure_title)
%   hpc_epoch_to_erp_plot_single_subject(input_folder, trial_type_values, channels, output_dir, figure_title, time_window_ms)
%   hpc_epoch_to_erp_plot_single_subject(..., condition_order_source)
%
% Inputs:
%   input_folder      - Folder containing parameter_<N>_single_subjects_Mean.mat files.
%   trial_type_values - Trial type labels to plot. Accepts CSV/semicolon list, string array,
%                       or cell array (e.g., {'Study_hits','Study_misses'}).
%   channels          - Required channel indices.
%                       Accepts numeric vector or CSV string (e.g., '21,22').
%   output_dir        - Optional output folder (default: <input_folder>/erp_plots_single_subject).
%   figure_title      - Optional custom figure title.
%   time_window_ms    - Optional time windows (default: none).
%                       Accepted formats:
%                       - single window: [300 500] or '300-500'
%                       - multiple windows: [300 500; 600 800] or '300-500;600-800'
%                       If provided and exactly 2 conditions are requested, paired t-tests
%                       are run on subject-level average voltage per window.
%   condition_order_source - Optional condition-order path:
%                       - file: condition_order.mat / condition_order.txt
%                       - folder containing condition_order files
%                       If omitted, the newest condition_order file under input_folder is used.
%
% Outputs:
%   - <output_dir>/grand_average_erp_single_subject_<trial_types>.png
%   - <output_dir>/grand_average_erp_single_subject_<trial_types>.svg (best effort on headless HPC)
%   - <output_dir>/grand_average_erp_single_subject_<trial_types>.mat
%   - <output_dir>/grand_average_erp_single_subject_<trial_types>_stats_<timestamp>.txt

if nargin < 2 || isempty(trial_type_values)
    error('trial_type_values is required (CSV, string array, or cell array).');
end
if nargin < 3 || isempty(channels)
    error('channels is required (numeric vector or CSV string).');
end
if nargin < 4 || isempty(output_dir)
    output_dir = fullfile(input_folder, 'erp_plots_single_subject');
end
if nargin < 5 || isempty(figure_title)
    figure_title = 'Grand Average ERP by trial_type (single-subject parameter files)';
end
if nargin < 6
    time_window_ms = [];
end
if nargin < 7
    condition_order_source = '';
end

trial_type_values = parse_trial_type_values(trial_type_values);
channels = parse_channel_list(channels);
time_windows_ms = parse_time_windows(time_window_ms);

fprintf('\n==============================================\n');
fprintf('  HPC SINGLE-SUBJECT PARAMETER -> ERP PLOT\n');
fprintf('==============================================\n');
fprintf('Input folder:   %s\n', input_folder);
fprintf('Trial types:    %s\n', strjoin(trial_type_values, ', '));
fprintf('Channels:       %s\n', num2str(channels));
fprintf('Output folder:  %s\n', output_dir);
if ~isempty(condition_order_source)
    fprintf('Condition map:  %s\n', condition_order_source);
else
    fprintf('Condition map:  (auto-detect)\n');
end
if ~isempty(time_windows_ms)
    fprintf('Time windows:\n');
    for w = 1:size(time_windows_ms, 1)
        fprintf('  - %.3f to %.3f\n', time_windows_ms(w, 1), time_windows_ms(w, 2));
    end
else
    fprintf('Time windows:   (none)\n');
end
fprintf('==============================================\n\n');

if ~exist(input_folder, 'dir')
    error('Input folder does not exist: %s', input_folder);
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

set(0, 'DefaultFigureVisible', 'off');
java.lang.System.setProperty('java.awt.headless', 'true');

n_conditions = numel(trial_type_values);
all_data_by_condition = cell(n_conditions, 1);
time_vector = [];
xlabel_text = 'Time (ms)';
channels_in_use = [];
n_timepoints = [];
n_subjects = [];
subject_ids = {};
trial_counts = [];
parameter_indices = nan(1, n_conditions);
parameter_files = cell(1, n_conditions);

[condition_order, condition_order_path] = resolve_condition_order(input_folder, condition_order_source);
fprintf('Resolved condition order from: %s\n', condition_order_path);
for i = 1:numel(condition_order)
    fprintf('  %d -> %s\n', i, condition_order{i});
end

for c = 1:n_conditions
    label = trial_type_values{c};
    param_idx = find_condition_index(condition_order, label);
    if isempty(param_idx)
        error(['Requested label "%s" was not found in condition_order (%s). ' ...
               'Available labels: %s'], ...
            label, condition_order_path, strjoin(condition_order, ', '));
    end

    param_file = fullfile(input_folder, sprintf('parameter_%d_single_subjects_Mean.mat', param_idx));
    if ~exist(param_file, 'file')
        error(['Condition "%s" maps to parameter %d by condition_order, but file is missing: %s'], ...
            label, param_idx, param_file);
    end

    parameter_indices(c) = param_idx;
    parameter_files{c} = param_file;
    fprintf('Condition %-30s -> parameter_%d_single_subjects_Mean.mat\n', label, param_idx);

    loaded = load(param_file);
    [condition_data, cond_time_vector, cond_xlabel] = extract_single_subject_condition_data(loaded, param_file);

    if isempty(channels_in_use)
        n_subjects = size(condition_data, 3);
        n_timepoints = size(condition_data, 2);
        channels_in_use = validate_channels(channels, size(condition_data, 1));
        subject_ids = build_subject_ids(n_subjects);
        trial_counts = ones(n_subjects, n_conditions);

        for k = 1:n_conditions
            all_data_by_condition{k} = nan(numel(channels_in_use), n_timepoints, n_subjects);
        end

        if ~isempty(cond_time_vector) && numel(cond_time_vector) == n_timepoints
            time_vector = cond_time_vector;
            xlabel_text = cond_xlabel;
        else
            time_vector = 1:n_timepoints;
            xlabel_text = 'Time Points';
        end
        fprintf('Detected %d subject(s) x %d timepoint(s)\n', n_subjects, n_timepoints);
    else
        if size(condition_data, 2) ~= n_timepoints
            error(['Time-point mismatch for condition "%s": expected %d points, found %d. ' ...
                   'All parameter files must share the same timeline.'], ...
                label, n_timepoints, size(condition_data, 2));
        end
        if size(condition_data, 3) ~= n_subjects
            error(['Subject-count mismatch for condition "%s": expected %d subjects, found %d. ' ...
                   'All parameter files must share the same subject count.'], ...
                label, n_subjects, size(condition_data, 3));
        end
        if any(channels_in_use > size(condition_data, 1))
            error(['Channel index out of range for condition "%s": requested %s but file has %d channels.'], ...
                label, num2str(channels_in_use), size(condition_data, 1));
        end
    end

    all_data_by_condition{c} = condition_data(channels_in_use, :, :);
end

grand_average_by_condition = cell(n_conditions, 1);
participants_per_condition = zeros(1, n_conditions);
for c = 1:n_conditions
    valid_subject_mask = squeeze(any(any(isfinite(all_data_by_condition{c}), 1), 2));
    if isrow(valid_subject_mask)
        valid_subject_mask = valid_subject_mask';
    end
    participants_per_condition(c) = sum(valid_subject_mask);
    trial_counts(:, c) = double(valid_subject_mask);
    grand_average_by_condition{c} = mean(all_data_by_condition{c}, 3, 'omitnan');
end

n_windows = size(time_windows_ms, 1);
window_indices = cell(n_windows, 1);
resolved_time_windows = nan(n_windows, 2);
for w = 1:n_windows
    [idx, resolved_window] = resolve_time_window_indices(time_vector, time_windows_ms(w, :));
    if isempty(idx)
        error('Requested time window %.3f-%.3f has no overlap with the data timeline.', ...
            time_windows_ms(w, 1), time_windows_ms(w, 2));
    end
    window_indices{w} = idx;
    resolved_time_windows(w, :) = resolved_window;
end

subject_window_means = nan(n_subjects, n_conditions, n_windows);
for w = 1:n_windows
    idx = window_indices{w};
    for c = 1:n_conditions
        subject_window_means(:, c, w) = compute_subject_window_means(all_data_by_condition{c}, idx);
    end
end

stats_by_window = repmat(init_stats_struct(), n_windows, 1);
for w = 1:n_windows
    stats_by_window(w) = run_window_stats_for_window(subject_window_means(:, :, w), subject_ids, ...
        trial_type_values, resolved_time_windows(w, :));
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
    error('No requested conditions could be plotted from the single-subject parameter files.');
end

if strcmp(xlabel_text, 'Time (ms)')
    xline(0, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.0);
end
yline(0, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);
grid on;
box off;

if n_windows > 0
    y_limits = ylim;
    for w = 1:n_windows
        idx = window_indices{w};
        x1 = time_vector(idx(1));
        x2 = time_vector(idx(end));
        highlight = patch([x1 x2 x2 x1], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
            [0.93 0.89 0.65], 'FaceAlpha', 0.16, 'EdgeColor', 'none');
        uistack(highlight, 'bottom');
    end
end

xlabel(xlabel_text);
ylabel('Voltage (\muV)');
title(figure_title, 'Interpreter', 'none');
legend(legend_entries, 'Location', 'best', 'Interpreter', 'none');

if n_windows > 0
    stat_text = build_stats_annotation(stats_by_window);
    [stat_color, stat_bg] = stats_annotation_colors(stats_by_window);
    x_limits = xlim;
    y_limits = ylim;
    text_x = x_limits(1) + 0.02 * (x_limits(2) - x_limits(1));
    text_y = y_limits(2) - 0.04 * (y_limits(2) - y_limits(1));
    text(text_x, text_y, stat_text, ...
        'Color', stat_color, ...
        'BackgroundColor', stat_bg, ...
        'EdgeColor', [0.8 0.8 0.8], ...
        'Margin', 5, ...
        'VerticalAlignment', 'top', ...
        'Interpreter', 'none', ...
        'FontSize', 10, ...
        'FontWeight', 'bold');
end

hold off;

trial_type_tag = build_trial_type_tag(trial_type_values);
output_png = fullfile(output_dir, sprintf('grand_average_erp_single_subject_%s.png', trial_type_tag));
output_svg = fullfile(output_dir, sprintf('grand_average_erp_single_subject_%s.svg', trial_type_tag));
output_mat = fullfile(output_dir, sprintf('grand_average_erp_single_subject_%s.mat', trial_type_tag));
timestamp_tag = datestr(now, 'yyyymmdd_HHMMSS');
output_stats_txt = fullfile(output_dir, sprintf('grand_average_erp_single_subject_%s_stats_%s.txt', ...
    trial_type_tag, timestamp_tag));

exportgraphics(fig, output_png, 'Resolution', 300);
svg_saved = save_svg_with_fallback(fig, output_svg);
close(fig);

save(output_mat, 'trial_type_values', 'channels_in_use', 'time_vector', ...
    'grand_average_by_condition', 'trial_counts', 'participants_per_condition', ...
    'subject_ids', 'time_windows_ms', 'resolved_time_windows', 'window_indices', ...
    'subject_window_means', 'stats_by_window', 'condition_order', ...
    'condition_order_path', 'parameter_indices', 'parameter_files');

write_stats_report(output_stats_txt, input_folder, trial_type_values, channels_in_use, ...
    participants_per_condition, subject_ids, time_windows_ms, resolved_time_windows, ...
    window_indices, subject_window_means, stats_by_window, output_png, output_svg, ...
    output_mat, svg_saved);

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
fprintf('Saved STATS TXT: %s\n', output_stats_txt);
fprintf('==============================================\n\n');

if nargout > 0
    results = struct();
    results.trial_type_values = trial_type_values;
    results.condition_order = condition_order;
    results.condition_order_path = condition_order_path;
    results.parameter_indices = parameter_indices;
    results.parameter_files = parameter_files;
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
    results.output_stats_txt = output_stats_txt;
    results.time_windows_ms = time_windows_ms;
    results.resolved_time_windows = resolved_time_windows;
    results.window_sample_indices = window_indices;
    results.subject_window_means = subject_window_means;
    results.stats_by_window = stats_by_window;
end

end

function [condition_order, resolved_source] = resolve_condition_order(input_folder, condition_order_source)
condition_order = {};
resolved_source = '';

if nargin >= 2 && ~isempty(condition_order_source)
    [condition_order, resolved_source] = load_condition_order_any(condition_order_source);
    if isempty(condition_order)
        error('Could not read condition_order from: %s', condition_order_source);
    end
    return;
end

candidate_struct = struct('path', {}, 'datenum', {});

direct_candidates = { ...
    fullfile(input_folder, 'condition_order.mat'), ...
    fullfile(input_folder, 'condition_order.txt'), ...
    fullfile(input_folder, 'limo_first_level', 'condition_order.mat'), ...
    fullfile(input_folder, 'limo_first_level', 'condition_order.txt'), ...
    fullfile(input_folder, 'epoch', 'limo_first_level', 'condition_order.mat'), ...
    fullfile(input_folder, 'epoch', 'limo_first_level', 'condition_order.txt') ...
};

for i = 1:numel(direct_candidates)
    f = direct_candidates{i};
    if exist(f, 'file')
        info = dir(f);
        candidate_struct(end+1).path = f; %#ok<AGROW>
        candidate_struct(end).datenum = info.datenum;
    end
end

recursive_mat = dir(fullfile(input_folder, '**', 'condition_order.mat'));
for i = 1:numel(recursive_mat)
    f = fullfile(recursive_mat(i).folder, recursive_mat(i).name);
    candidate_struct(end+1).path = f; %#ok<AGROW>
    candidate_struct(end).datenum = recursive_mat(i).datenum;
end

recursive_txt = dir(fullfile(input_folder, '**', 'condition_order.txt'));
for i = 1:numel(recursive_txt)
    f = fullfile(recursive_txt(i).folder, recursive_txt(i).name);
    candidate_struct(end+1).path = f; %#ok<AGROW>
    candidate_struct(end).datenum = recursive_txt(i).datenum;
end

if isempty(candidate_struct)
    error(['No condition_order.mat/txt found under %s. ' ...
           'Provide condition_order_source explicitly.'], input_folder);
end

[~, ord] = sort([candidate_struct.datenum], 'descend');
candidate_struct = candidate_struct(ord);
candidate_paths = unique({candidate_struct.path}, 'stable');

for i = 1:numel(candidate_paths)
    [condition_order, resolved_source] = load_condition_order_any(candidate_paths{i});
    if ~isempty(condition_order)
        return;
    end
end

error(['Found condition_order candidate files, but none could be parsed. ' ...
       'Please provide condition_order_source explicitly.']);
end

function [condition_order, resolved_source] = load_condition_order_any(path_value)
condition_order = {};
resolved_source = '';

if exist(path_value, 'file')
    condition_order = load_condition_order_file(path_value);
    if ~isempty(condition_order)
        resolved_source = path_value;
    end
    return;
end

if exist(path_value, 'dir')
    [condition_order, resolved_source] = load_condition_order_from_folder(path_value);
    return;
end

error('condition_order_source does not exist: %s', path_value);
end

function [condition_order, resolved_source] = load_condition_order_from_folder(folder_path)
condition_order = {};
resolved_source = '';

candidate_files = { ...
    fullfile(folder_path, 'condition_order.mat'), ...
    fullfile(folder_path, 'condition_order.txt'), ...
    fullfile(folder_path, 'limo_first_level', 'condition_order.mat'), ...
    fullfile(folder_path, 'limo_first_level', 'condition_order.txt'), ...
    fullfile(folder_path, 'epoch', 'limo_first_level', 'condition_order.mat'), ...
    fullfile(folder_path, 'epoch', 'limo_first_level', 'condition_order.txt') ...
};

for i = 1:numel(candidate_files)
    f = candidate_files{i};
    if exist(f, 'file')
        condition_order = load_condition_order_file(f);
        if ~isempty(condition_order)
            resolved_source = f;
            return;
        end
    end
end

recursive_hits = [ ...
    dir(fullfile(folder_path, '**', 'condition_order.mat')); ...
    dir(fullfile(folder_path, '**', 'condition_order.txt')) ...
];

if isempty(recursive_hits)
    return;
end

[~, ord] = sort([recursive_hits.datenum], 'descend');
recursive_hits = recursive_hits(ord);

for i = 1:numel(recursive_hits)
    f = fullfile(recursive_hits(i).folder, recursive_hits(i).name);
    condition_order = load_condition_order_file(f);
    if ~isempty(condition_order)
        resolved_source = f;
        return;
    end
end
end

function condition_order = load_condition_order_file(file_path)
condition_order = {};
[~, ~, ext] = fileparts(file_path);
ext = lower(ext);

switch ext
    case '.mat'
        data = load(file_path);
        if isfield(data, 'condition_order')
            condition_order = normalize_condition_order_values(data.condition_order);
        end
    case '.txt'
        condition_order = parse_condition_order_text_file(file_path);
    otherwise
        condition_order = {};
end
end

function labels = parse_condition_order_text_file(file_path)
labels = {};
fid = fopen(file_path, 'r');
if fid == -1
    return;
end
cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line) || startsWith(line, '#')
        continue;
    end

    parts = strsplit(line, sprintf('\t'));
    if numel(parts) >= 2
        label = strtrim(strjoin(parts(2:end), ' '));
    else
        tokens = regexp(line, '^\s*\d+\s+(.+)$', 'tokens', 'once');
        if ~isempty(tokens)
            label = strtrim(tokens{1});
        else
            label = strtrim(line);
        end
    end

    label = normalize_label_value(label);
    if ~isempty(label)
        labels{end+1} = label; %#ok<AGROW>
    end
end

labels = unique(labels, 'stable');
end

function labels = normalize_condition_order_values(raw_values)
labels = {};
if isempty(raw_values)
    return;
end

if ischar(raw_values) || (isstring(raw_values) && isscalar(raw_values))
    tmp = strrep(char(raw_values), sprintf('\n'), ',');
    parts = strsplit(tmp, ',');
    for i = 1:numel(parts)
        label = normalize_label_value(parts{i});
        if ~isempty(label)
            labels{end+1} = label; %#ok<AGROW>
        end
    end
elseif isstring(raw_values)
    for i = 1:numel(raw_values)
        label = normalize_label_value(raw_values(i));
        if ~isempty(label)
            labels{end+1} = label; %#ok<AGROW>
        end
    end
elseif iscell(raw_values)
    for i = 1:numel(raw_values)
        label = normalize_label_value(raw_values{i});
        if ~isempty(label)
            labels{end+1} = label; %#ok<AGROW>
        end
    end
end

labels = unique(labels, 'stable');
end

function idx = find_condition_index(condition_order, label)
idx = find(strcmp(condition_order, label), 1, 'first');
if isempty(idx)
    idx_case = find(strcmpi(condition_order, label), 1, 'first');
    if ~isempty(idx_case)
        idx = idx_case;
    end
end
end

function [condition_data, time_vector, xlabel_text] = extract_single_subject_condition_data(loaded, source_file)
if ~isfield(loaded, 'Data') || ~isstruct(loaded.Data) || ~isfield(loaded.Data, 'data')
    error('Missing Data.data in file: %s', source_file);
end

raw_data = loaded.Data.data;
if ~isnumeric(raw_data)
    error('Data.data is not numeric in file: %s', source_file);
end

raw_data = double(raw_data);
if ndims(raw_data) == 4
    if size(raw_data, 3) ~= 1
        error(['Expected single-subject parameter data with singleton 3rd dimension in %s. ' ...
               'Found size: %s'], source_file, mat2str(size(raw_data)));
    end
    raw_data = squeeze(raw_data);
end

if ndims(raw_data) == 2
    raw_data = reshape(raw_data, size(raw_data, 1), size(raw_data, 2), 1);
elseif ndims(raw_data) ~= 3
    error('Unexpected Data.data shape in %s: %s', source_file, mat2str(size(raw_data)));
end

condition_data = raw_data;
[time_vector, xlabel_text] = infer_time_vector_from_loaded_data(loaded, size(condition_data, 2));
end

function [time_vector, xlabel_text] = infer_time_vector_from_loaded_data(loaded, n_timepoints)
time_vector = [];
xlabel_text = 'Time (ms)';

if isfield(loaded, 'Data') && isstruct(loaded.Data) && ...
        isfield(loaded.Data, 'limo') && isstruct(loaded.Data.limo) && ...
        isfield(loaded.Data.limo, 'data') && isstruct(loaded.Data.limo.data)

    limo_data = loaded.Data.limo.data;

    if isfield(limo_data, 'timevect') && isnumeric(limo_data.timevect) && ...
            numel(limo_data.timevect) == n_timepoints
        time_vector = double(limo_data.timevect(:))';
        return;
    end

    if isfield(limo_data, 'start') && isfield(limo_data, 'sampling_rate') && ...
            isnumeric(limo_data.start) && isnumeric(limo_data.sampling_rate) && ...
            isfinite(double(limo_data.sampling_rate)) && double(limo_data.sampling_rate) > 0
        step_ms = 1000 / double(limo_data.sampling_rate);
        time_vector = double(limo_data.start) + (0:(n_timepoints - 1)) * step_ms;
        return;
    end

    if isfield(limo_data, 'start') && isfield(limo_data, 'end') && ...
            isnumeric(limo_data.start) && isnumeric(limo_data.end) && ...
            isfinite(double(limo_data.start)) && isfinite(double(limo_data.end))
        time_vector = linspace(double(limo_data.start), double(limo_data.end), n_timepoints);
        return;
    end
end

time_vector = 1:n_timepoints;
xlabel_text = 'Time Points';
end

function subject_ids = build_subject_ids(n_subjects)
subject_ids = cell(n_subjects, 1);
for i = 1:n_subjects
    subject_ids{i} = sprintf('S%03d', i);
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

function windows = parse_time_windows(value)
if isempty(value)
    windows = [];
    return;
end

if isnumeric(value)
    windows = parse_time_windows_numeric(value);
elseif ischar(value) || (isstring(value) && isscalar(value))
    windows = parse_time_windows_string(char(value));
elseif isstring(value)
    windows = [];
    for i = 1:numel(value)
        windows = [windows; parse_time_windows_string(char(value(i)))]; %#ok<AGROW>
    end
elseif iscell(value)
    windows = [];
    for i = 1:numel(value)
        windows = [windows; parse_time_windows(value{i})]; %#ok<AGROW>
    end
else
    error('time_window_ms must be empty, numeric, or a string like "300-500;600-800".');
end

if isempty(windows)
    return;
end
if size(windows, 2) ~= 2
    error('time_window_ms must resolve to an N-by-2 list of [start end] windows.');
end
if any(~isfinite(windows(:)))
    error('time_window_ms contains non-numeric values.');
end

windows = sort(double(windows), 2);
if any(windows(:, 1) == windows(:, 2))
    error('time_window_ms windows must have different start/end values.');
end
end

function windows = parse_time_windows_numeric(value)
vals = double(value);
if isempty(vals)
    windows = [];
    return;
end

if isvector(vals)
    vals = vals(:)';
    if mod(numel(vals), 2) ~= 0
        error('Numeric time_window_ms vector must contain an even number of values.');
    end
    windows = reshape(vals, 2, [])';
    return;
end

if size(vals, 2) == 2
    windows = vals;
    return;
end
if size(vals, 1) == 2
    windows = vals';
    return;
end

error('Numeric time_window_ms must be [start end], N-by-2, or a vector of paired values.');
end

function windows = parse_time_windows_string(raw)
raw = strtrim(raw);
if isempty(raw)
    windows = [];
    return;
end

raw = strrep(raw, sprintf('\n'), ';');
segments = strsplit(raw, ';');
windows = [];

for s = 1:numel(segments)
    segment = strtrim(segments{s});
    if isempty(segment)
        continue;
    end

    segment = regexprep(segment, '(?i)\bms\b', '');
    segment = regexprep(segment, '(?i)\bto\b', '-');

    pair_tokens = regexp(segment, '^\s*([+-]?\d*\.?\d+)\s*[-,:]\s*([+-]?\d*\.?\d+)\s*$', ...
        'tokens', 'once');
    if isempty(pair_tokens)
        pair_tokens = regexp(segment, '^\s*([+-]?\d*\.?\d+)\s+([+-]?\d*\.?\d+)\s*$', ...
            'tokens', 'once');
    end

    if ~isempty(pair_tokens)
        w = [str2double(pair_tokens{1}) str2double(pair_tokens{2})];
        if any(~isfinite(w))
            error('Could not parse time_window_ms segment: "%s".', segment);
        end
        windows(end+1, :) = w; %#ok<AGROW>
        continue;
    end

    number_tokens = regexp(segment, '[+-]?\d*\.?\d+', 'match');
    if isempty(number_tokens) || mod(numel(number_tokens), 2) ~= 0
        error(['Could not parse time_window_ms segment: "%s". ' ...
               'Use "start-end;start-end" (example "300-500;600-800").'], segment);
    end

    numbers = str2double(number_tokens);
    if any(~isfinite(numbers))
        error('Could not parse numeric values in time_window_ms segment: "%s".', segment);
    end
    windows = [windows; reshape(numbers, 2, [])']; %#ok<AGROW>
end
end

function [idx, resolved_window] = resolve_time_window_indices(time_vector, window)
idx = [];
resolved_window = [];

if isempty(window)
    return;
end

idx = find(time_vector >= window(1) & time_vector <= window(2));
if isempty(idx)
    return;
end

resolved_window = [time_vector(idx(1)) time_vector(idx(end))];
end

function means = compute_subject_window_means(condition_data, window_idx)
n_subjects_local = size(condition_data, 3);
means = nan(n_subjects_local, 1);
for s = 1:n_subjects_local
    subject_slice = condition_data(:, window_idx, s);
    means(s) = mean(subject_slice(:), 'omitnan');
end
end

function stats_info = init_stats_struct()
stats_info = struct();
stats_info.requested = false;
stats_info.performed = false;
stats_info.significant = false;
stats_info.alpha = 0.05;
stats_info.reason = '';
stats_info.p_value = nan;
stats_info.t_stat = nan;
stats_info.df = nan;
stats_info.mean_condition_1 = nan;
stats_info.mean_condition_2 = nan;
stats_info.mean_diff = nan;
stats_info.std_diff = nan;
stats_info.cohens_dz = nan;
stats_info.n_pairs = 0;
stats_info.condition_1 = '';
stats_info.condition_2 = '';
stats_info.included_subject_ids = {};
stats_info.excluded_subject_ids = {};
stats_info.requested_time_window = [nan nan];
stats_info.resolved_time_window = [nan nan];
end

function stats_info = run_window_stats_for_window(subject_window_means, subject_ids, ...
    trial_type_values, resolved_window)
stats_info = init_stats_struct();
stats_info.requested = true;
stats_info.requested_time_window = resolved_window;
stats_info.resolved_time_window = resolved_window;

n_conditions_local = size(subject_window_means, 2);
if n_conditions_local ~= 2
    stats_info.reason = sprintf('t-test requires exactly 2 conditions (got %d).', n_conditions_local);
    return;
end

stats_info.condition_1 = trial_type_values{1};
stats_info.condition_2 = trial_type_values{2};

v1 = subject_window_means(:, 1);
v2 = subject_window_means(:, 2);
valid_mask = isfinite(v1) & isfinite(v2);

stats_info.n_pairs = sum(valid_mask);
stats_info.included_subject_ids = subject_ids(valid_mask);
stats_info.excluded_subject_ids = subject_ids(~valid_mask);

if stats_info.n_pairs < 2
    stats_info.reason = sprintf('Need at least 2 paired subjects for t-test (found %d).', stats_info.n_pairs);
    return;
end

paired_1 = v1(valid_mask);
paired_2 = v2(valid_mask);

stats_info.mean_condition_1 = mean(paired_1, 'omitnan');
stats_info.mean_condition_2 = mean(paired_2, 'omitnan');
diffs = paired_1 - paired_2;
stats_info.mean_diff = mean(diffs, 'omitnan');
stats_info.std_diff = std(diffs, 0, 'omitnan');
if stats_info.std_diff > 0
    stats_info.cohens_dz = stats_info.mean_diff / stats_info.std_diff;
end

try
    [h, p, ~, stats] = ttest(paired_1, paired_2, 'Alpha', stats_info.alpha);
    stats_info.performed = true;
    stats_info.significant = logical(h);
    stats_info.p_value = p;
    stats_info.t_stat = stats.tstat;
    stats_info.df = stats.df;
    stats_info.reason = 'Paired t-test completed.';
catch ME
    stats_info.reason = sprintf('ttest failed: %s', ME.message);
end
end

function txt = build_stats_annotation(stats_by_window)
if isempty(stats_by_window)
    txt = 'No time window requested; no statistics.';
    return;
end

lines = cell(numel(stats_by_window), 1);
for w = 1:numel(stats_by_window)
    stats_info = stats_by_window(w);

    if all(isfinite(stats_info.resolved_time_window))
        window_label = sprintf('%.1f-%.1f ms', ...
            stats_info.resolved_time_window(1), stats_info.resolved_time_window(2));
    else
        window_label = sprintf('Window %d', w);
    end

    if ~stats_info.performed
        lines{w} = sprintf('%s: %s', window_label, stats_info.reason);
        continue;
    end

    if stats_info.significant
        sig_word = 'sig';
    else
        sig_word = 'n.s.';
    end

    lines{w} = sprintf('%s: t(%d)=%.3f, p=%.4g, n=%d (%s)', ...
        window_label, round(stats_info.df), stats_info.t_stat, ...
        stats_info.p_value, stats_info.n_pairs, sig_word);
end

txt = strjoin(lines, sprintf('\n'));
end

function [txt_color, bg_color] = stats_annotation_colors(stats_by_window)
if isempty(stats_by_window)
    txt_color = [0.22 0.22 0.22];
    bg_color = [0.95 0.95 0.95];
    return;
end

performed_mask = arrayfun(@(s) s.performed, stats_by_window);
significant_mask = arrayfun(@(s) s.performed && s.significant, stats_by_window);

if any(significant_mask)
    txt_color = [0.05 0.35 0.10];
    bg_color = [0.91 0.98 0.91];
elseif any(performed_mask)
    txt_color = [0.45 0.22 0.06];
    bg_color = [0.99 0.95 0.90];
else
    txt_color = [0.22 0.22 0.22];
    bg_color = [0.95 0.95 0.95];
end
end

function write_stats_report(output_path, input_folder, trial_type_values, channels_in_use, ...
    participants_per_condition, subject_ids, time_windows_ms, resolved_time_windows, ...
    window_indices, subject_window_means, stats_by_window, output_png, output_svg, ...
    output_mat, svg_saved)
fid = fopen(output_path, 'w');
if fid == -1
    error('Could not create stats report file: %s', output_path);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'ERP Window Statistics Report\n');
fprintf(fid, '============================\n');
fprintf(fid, 'Generated at: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, 'Input folder: %s\n', input_folder);
fprintf(fid, 'Conditions: %s\n', strjoin(trial_type_values, ', '));
fprintf(fid, 'Channels: %s\n', num2str(channels_in_use));
fprintf(fid, 'Total subjects considered: %d\n', numel(subject_ids));
fprintf(fid, '\nParticipants with data by condition:\n');
for c = 1:numel(trial_type_values)
    fprintf(fid, '  - %s: %d\n', trial_type_values{c}, participants_per_condition(c));
end

if isempty(time_windows_ms)
    fprintf(fid, '\nTime windows: not provided\n');
else
    fprintf(fid, '\nRequested time windows (ms):\n');
    for w = 1:size(time_windows_ms, 1)
        fprintf(fid, '  - [%.4f, %.4f]\n', time_windows_ms(w, 1), time_windows_ms(w, 2));
    end
end

n_windows = size(subject_window_means, 3);
if n_windows < numel(window_indices)
    n_windows = numel(window_indices);
end
if n_windows < numel(stats_by_window)
    n_windows = numel(stats_by_window);
end

if n_windows == 0
    fprintf(fid, '\nNo per-window statistics were requested.\n');
else
    for w = 1:n_windows
        fprintf(fid, '\n----------------------------------------\n');
        fprintf(fid, 'Window %d\n', w);

        if size(time_windows_ms, 1) >= w
            fprintf(fid, '  requested: %.4f to %.4f ms\n', ...
                time_windows_ms(w, 1), time_windows_ms(w, 2));
        end

        if size(resolved_time_windows, 1) >= w
            fprintf(fid, '  resolved: %.4f to %.4f ms\n', ...
                resolved_time_windows(w, 1), resolved_time_windows(w, 2));
        end

        if numel(window_indices) >= w && ~isempty(window_indices{w})
            fprintf(fid, '  samples in window: %d\n', numel(window_indices{w}));
        else
            fprintf(fid, '  samples in window: 0\n');
        end

        if size(subject_window_means, 3) >= w
            fprintf(fid, '  subject-level means (uV):\n');
            for c = 1:numel(trial_type_values)
                vals = subject_window_means(:, c, w);
                valid_vals = vals(isfinite(vals));
                fprintf(fid, '    - %s: n=%d, mean=%.6f, sd=%.6f\n', ...
                    trial_type_values{c}, numel(valid_vals), ...
                    mean(valid_vals, 'omitnan'), std(valid_vals, 0, 'omitnan'));
            end
        end

        if numel(stats_by_window) >= w
            stats_info = stats_by_window(w);
        else
            stats_info = init_stats_struct();
            stats_info.reason = 'No stats struct available for this window.';
        end

        fprintf(fid, '  statistical test:\n');
        fprintf(fid, '    requested: %s\n', logical_to_yesno(stats_info.requested));
        fprintf(fid, '    performed: %s\n', logical_to_yesno(stats_info.performed));
        fprintf(fid, '    reason: %s\n', stats_info.reason);

        if stats_info.performed
            fprintf(fid, '    test: paired t-test\n');
            fprintf(fid, '    alpha: %.3f\n', stats_info.alpha);
            fprintf(fid, '    condition_1: %s\n', stats_info.condition_1);
            fprintf(fid, '    condition_2: %s\n', stats_info.condition_2);
            fprintf(fid, '    n_pairs: %d\n', stats_info.n_pairs);
            fprintf(fid, '    t(%d)=%.8f\n', round(stats_info.df), stats_info.t_stat);
            fprintf(fid, '    p_value=%.10g\n', stats_info.p_value);
            fprintf(fid, '    significant: %s\n', logical_to_yesno(stats_info.significant));
            fprintf(fid, '    mean_%s=%.8f\n', stats_info.condition_1, stats_info.mean_condition_1);
            fprintf(fid, '    mean_%s=%.8f\n', stats_info.condition_2, stats_info.mean_condition_2);
            fprintf(fid, '    mean_diff_(%s-%s)=%.8f\n', ...
                stats_info.condition_1, stats_info.condition_2, stats_info.mean_diff);
            fprintf(fid, '    sd_diff=%.8f\n', stats_info.std_diff);
            fprintf(fid, '    cohens_dz=%.8f\n', stats_info.cohens_dz);
            fprintf(fid, '    included_subject_ids: %s\n', join_labels(stats_info.included_subject_ids));
            fprintf(fid, '    excluded_subject_ids: %s\n', join_labels(stats_info.excluded_subject_ids));
        end
    end
end

fprintf(fid, '\nOutput files:\n');
fprintf(fid, '  PNG: %s\n', output_png);
if svg_saved
    fprintf(fid, '  SVG: %s\n', output_svg);
else
    fprintf(fid, '  SVG: (not saved)\n');
end
fprintf(fid, '  MAT: %s\n', output_mat);
fprintf(fid, '  STATS_TXT: %s\n', output_path);
end

function out = logical_to_yesno(value)
if value
    out = 'yes';
else
    out = 'no';
end
end

function out = join_labels(values)
if isempty(values)
    out = '(none)';
    return;
end
tmp = values(:)';
out = strjoin(tmp, ', ');
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
    out = parse_channel_tokens(char(value));
elseif isstring(value)
    out = [];
    for i = 1:numel(value)
        out = [out parse_channel_tokens(char(value(i)))]; %#ok<AGROW>
    end
elseif iscell(value)
    out = [];
    for i = 1:numel(value)
        if isnumeric(value{i})
            out = [out double(value{i}(:)')]; %#ok<AGROW>
        else
            out = [out parse_channel_tokens(char(value{i}))]; %#ok<AGROW>
        end
    end
else
    error('channels must be numeric vector or a comma/semicolon-separated string.');
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

function out = parse_channel_tokens(raw)
parts = regexp(strtrim(raw), '[,;\s]+', 'split');
parts = parts(~cellfun(@isempty, parts));

if isempty(parts)
    out = [];
    return;
end

out = nan(1, numel(parts));
for i = 1:numel(parts)
    out(i) = str2double(strtrim(parts{i}));
end
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
