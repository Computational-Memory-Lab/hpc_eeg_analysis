function hpc_limo_channel_time_plots(input_folder, output_dir, title_base, channel_time_title, topoplot_title, lr_title, topoplot_layout_type)
% HPC_LIMO_CHANNEL_TIME_PLOTS - Generate channel-time and topoplot figures
%
% Usage:
%   hpc_limo_channel_time_plots(input_folder, output_dir)
%   hpc_limo_channel_time_plots(input_folder, output_dir, title_base)
%   hpc_limo_channel_time_plots(input_folder, output_dir, title_base, channel_time_title, topoplot_title)
%   hpc_limo_channel_time_plots(input_folder, output_dir, title_base, channel_time_title, topoplot_title, lr_title)
%   hpc_limo_channel_time_plots(input_folder, output_dir, title_base, channel_time_title, topoplot_title, lr_title, topoplot_layout_type)
%
% Inputs:
%   input_folder - Path to a limo_second_level_<output_tag> folder containing
%                  LIMO.mat and paired_samples_ttest_parameter_*.mat
%   output_dir   - Directory to save the output .png files
%   title_base   - Optional custom base contrast title used to build default
%                  titles, e.g. 'Test_hits vs Correct_rejection'
%   channel_time_title - Optional full custom channel-time figure title.
%   topoplot_title     - Optional full custom topoplot figure title.
%   lr_title           - Optional full custom likelihood-ratio figure title.
%   topoplot_layout_type - Optional topoplot layout:
%                        'grid'   (default): near-square grid layout
%                        'zigzag' : staggered left-to-right zigzag line
%                        'line'   : single horizontal row
%
% Outputs:
%   - <output_dir>/<test_name>_channel_time_plot.png      Signed -log10(p) channel-time plot
%   - <output_dir>/<test_name>_channel_time_plot.svg      Signed -log10(p) channel-time plot (SVG)
%   - <output_dir>/<test_name>_channel_time_plot_LR.png   Likelihood ratio plot (if data exists)
%   - <output_dir>/<test_name>_topoplots.png              Topoplot figure at multiple time points
%   - <output_dir>/<test_name>_topoplots.svg              Topoplot figure at multiple time points (SVG)
%
% The test_name is derived from the input_folder name.
% Example: limo_second_level_test_hits_vs_test_misses

if nargin < 3 || isempty(title_base)
    title_base = '';
end
if nargin < 4 || isempty(channel_time_title)
    channel_time_title = '';
end
if nargin < 5 || isempty(topoplot_title)
    topoplot_title = '';
end
if nargin < 6 || isempty(lr_title)
    lr_title = '';
end
if nargin < 7 || isempty(topoplot_layout_type)
    topoplot_layout_type = 'grid';
end

title_base = strtrim(char(title_base));
channel_time_title = strtrim(char(channel_time_title));
topoplot_title = strtrim(char(topoplot_title));
lr_title = strtrim(char(lr_title));
topoplot_layout_type = normalize_topoplot_layout_type(topoplot_layout_type);

% ==================== PARSE INPUT FOLDER ====================
[~, folder_name] = fileparts(input_folder);

% Find the paired_samples_ttest_parameter_*.mat file in the input folder
mat_files = dir(fullfile(input_folder, 'paired_samples_ttest_parameter_*.mat'));
mat_files = mat_files(~contains({mat_files.name}, '_likelihood'));
if isempty(mat_files)
    error('No paired_samples_ttest_parameter_*.mat found in: %s', input_folder);
end
tok = regexp(mat_files(1).name, 'paired_samples_ttest_parameter_(\d+)\.mat', 'tokens');
parameter_num = str2double(tok{1}{1});
test_name = folder_name;
display_contrast_name = derive_display_contrast_name(input_folder, folder_name);
if ~isempty(title_base)
    display_contrast_name = normalize_title_label(title_base);
end

default_channel_time_title = sprintf('%s Channel-Time Plot (TFCE-corrected)', display_contrast_name);
default_topoplot_title = sprintf('%s Topoplots (TFCE-corrected)', display_contrast_name);
default_lr_title = sprintf('%s Likelihood Ratios', display_contrast_name);

fprintf('Default Step 6 titles:\n');
fprintf('  Channel-time: %s\n', default_channel_time_title);
fprintf('  Topoplots:    %s\n', default_topoplot_title);
fprintf('  LR:           %s\n', default_lr_title);
fprintf('  Topoplot layout: %s\n', topoplot_layout_type);

% ==================== SETUP ====================
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('Starting EEGLAB in headless mode...\n');
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1');
eeglab nogui;

% ==================== LOAD DATA ====================
limo_path       = fullfile(input_folder, 'LIMO.mat');
base_data_path  = fullfile(input_folder, sprintf('paired_samples_ttest_parameter_%d.mat', parameter_num));
likelihood_data_path = fullfile(input_folder, sprintf('paired_samples_ttest_parameter_%d_likelihood.mat', parameter_num));

if ~exist(limo_path, 'file')
    error('LIMO.mat not found in: %s', input_folder);
end
if ~exist(base_data_path, 'file')
    error('paired_samples_ttest_parameter_%d.mat not found in: %s', parameter_num, input_folder);
end

[~, FileName, ext] = fileparts(base_data_path);
FileName = [FileName ext];

test_limo  = load(limo_path);
base_data  = load(base_data_path);
LIMO       = test_limo.LIMO;

% Load likelihood data if it exists
likelihood_data     = [];
likelihood_values   = [];
log10_likelihood_values = [];
has_likelihood_data = exist(likelihood_data_path, 'file');
if has_likelihood_data
    likelihood_data = load(likelihood_data_path);
end

% Get TFCE-corrected p-values
[tfce_p_values, ~, ~] = limo_stat_values(FileName, 0.05, 3, LIMO);

% Load TFCE scores to report thresholds
tfce_file    = fullfile(LIMO.dir, 'tfce', ['tfce_' FileName]);
H0_tfce_file = fullfile(LIMO.dir, 'H0', ['tfce_H0_' FileName]);

if exist(tfce_file, 'file') && exist(H0_tfce_file, 'file')
    tfce_data     = load(tfce_file);
    tfce_score    = tfce_data.tfce_score;
    H0_tfce_data  = load(H0_tfce_file);
    H0_tfce_score = H0_tfce_data.tfce_H0_score;

    [~, ~, bootstrap_threshold] = limo_max_correction(tfce_score, H0_tfce_score, 0.05);
    max_observed_tfce = max(tfce_score(:));

    if max_observed_tfce < bootstrap_threshold
        fprintf('%s: Max observed TFCE = %.3f, Bootstrap threshold = %.3f [NOT SIGNIFICANT]\n', ...
            test_name, max_observed_tfce, bootstrap_threshold);
    else
        fprintf('%s: Max observed TFCE = %.3f, Bootstrap threshold = %.3f [SIGNIFICANT]\n', ...
            test_name, max_observed_tfce, bootstrap_threshold);
    end
end

% Extract test statistics
base_array  = base_data.paired_samples;
df_values   = squeeze(base_array(:,:,3));
t_values    = squeeze(base_array(:,:,4));
p_values    = tfce_p_values;

% Extract likelihood data
if has_likelihood_data
    likelihood_array        = likelihood_data.paired_samples;
    likelihood_values       = real(squeeze(likelihood_array(:,:,6)));
    log10_likelihood_values = real(squeeze(likelihood_array(:,:,7)));
end

% ==================== CHANNEL-TIME PLOT ====================
output_file_png = fullfile(output_dir, sprintf('%s_channel_time_plot.png', test_name));
output_file_svg = fullfile(output_dir, sprintf('%s_channel_time_plot.svg', test_name));
alpha_threshold = 0.05;
alpha_nonsig    = 0.4;

signed_logp_data = -log10(p_values) .* sign(t_values);
alpha_mask_significance = ones(size(p_values));
alpha_mask_significance(p_values >= alpha_threshold) = alpha_nonsig;

abs_max_val = max(abs(signed_logp_data(:)));
if isempty(abs_max_val) || abs_max_val == 0
    abs_max_val = 1;
end
cc = limo_color_images([-abs_max_val, abs_max_val]);

times = linspace(LIMO.data.start, LIMO.data.end, size(t_values,2));
channel_time_raster_scale = 4;

fig_exp = figure('Visible', 'off', 'Position', [100 100 1500 700], 'Color', 'w', 'InvertHardcopy', 'off');
colormap(fig_exp, cc);

imgax_exp = axes('Position', [0.08 0.12 0.78 0.68]);
set(imgax_exp, 'Color', [0.9 0.9 0.9], 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'LineWidth', 1);

signed_logp_data = real(signed_logp_data);
num_channels = size(signed_logp_data, 1);
channels = 1:num_channels;

times_hr = linspace(times(1), times(end), max(2, size(signed_logp_data, 2) * channel_time_raster_scale));
channels_hr = linspace(1, num_channels, max(2, num_channels * channel_time_raster_scale));
[time_grid, channel_grid] = meshgrid(times, channels);
[time_grid_hr, channel_grid_hr] = meshgrid(times_hr, channels_hr);

signed_logp_display = interp2(time_grid, channel_grid, signed_logp_data, ...
    time_grid_hr, channel_grid_hr, 'linear');
alpha_mask_display = interp2(time_grid, channel_grid, alpha_mask_significance, ...
    time_grid_hr, channel_grid_hr, 'nearest');
alpha_mask_display = min(max(alpha_mask_display, 0), 1);

h_img_exp = imagesc(times_hr, channels_hr, signed_logp_display);
xlim([times(1), times(end)]);

set(h_img_exp, 'AlphaData', alpha_mask_display);
set(imgax_exp, 'Color', [0.9 0.9 0.9]);
caxis([-abs_max_val, abs_max_val]);

xlabel('Time (ms)', 'Color', 'k');
ylabel('Channel', 'Color', 'k');

set(imgax_exp, 'YDir', 'reverse');
ylim([1 num_channels]);
grid(imgax_exp, 'off');
set(imgax_exp, 'XColor', 'k', 'YColor', 'k');
hold on;

% Add vertical line at time=0
ylim_vals = get(imgax_exp, 'YLim');
plot(imgax_exp, [0 0], ylim_vals, 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'LineWidth', 1);

colormap(imgax_exp, cc);

h = colorbar(imgax_exp, 'Position', [0.88 0.12 0.04 0.68]);
title(h, 'signed -log10(p)', 'FontSize', 10, 'Color', 'k');
set(h, 'Color', 'k', 'Box', 'on', 'LineWidth', 1, 'TickLength', 0);

% Add significance threshold lines to colorbar
sig_threshold = -log10(alpha_threshold);
cbar_pos    = get(h, 'Position');
cbar_limits = get(h, 'Limits');

y_pos = NaN;
y_neg = NaN;

if sig_threshold >= cbar_limits(1) && sig_threshold <= cbar_limits(2)
    norm_pos_pos = (sig_threshold - cbar_limits(1)) / (cbar_limits(2) - cbar_limits(1));
    y_pos = cbar_pos(2) + norm_pos_pos * cbar_pos(4);
end
if -sig_threshold >= cbar_limits(1) && -sig_threshold <= cbar_limits(2)
    norm_pos_neg = (-sig_threshold - cbar_limits(1)) / (cbar_limits(2) - cbar_limits(1));
    y_neg = cbar_pos(2) + norm_pos_neg * cbar_pos(4);
end

if isnan(y_neg) && isnan(y_pos)
    y_neg = cbar_pos(2);
    y_pos = cbar_pos(2) + cbar_pos(4);
end

if ~isnan(y_neg) && ~isnan(y_pos)
    border_inset   = 0.0006;
    overlay_x      = cbar_pos(1) + border_inset;
    overlay_y      = y_neg;
    overlay_width  = cbar_pos(3) - 2*border_inset;
    overlay_height = y_pos - y_neg;
    annotation('rectangle', [overlay_x, overlay_y, overlay_width, overlay_height], ...
        'FaceColor', [0.9137, 0.9020, 0.9020], ...
        'FaceAlpha', 1 - alpha_nonsig, ...
        'EdgeColor', 'none');
end
if ~isnan(y_pos)
    annotation('line', [cbar_pos(1), cbar_pos(1) + cbar_pos(3)], [y_pos, y_pos], ...
        'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);
end
if ~isnan(y_neg)
    annotation('line', [cbar_pos(1), cbar_pos(1) + cbar_pos(3)], [y_neg, y_neg], ...
        'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);
end

% Formatting
set(imgax_exp, 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'LineWidth', 1, 'TickLength', [0 0]);

chan_labels = {};
does_not_begin_with_E = false;
if isfield(LIMO.data, 'chanlocs') && ~isempty(LIMO.data.chanlocs)
    chan_labels = {LIMO.data.chanlocs.labels};
    does_not_begin_with_E = all(cellfun(@(x) isempty(regexp(x, '^E', 'once')), chan_labels));
end

if does_not_begin_with_E
    ytick_values = 1:num_channels;
else
    ytick_values = round(linspace(1, num_channels, 10));
end
set(imgax_exp, 'YTick', ytick_values);

if ~isempty(chan_labels) && length(chan_labels) == num_channels
    set(imgax_exp, 'YTickLabel', chan_labels(ytick_values));
end

tick_interval = 200;
min_tick = ceil(times(1) / tick_interval) * tick_interval;
max_tick = floor(times(end) / tick_interval) * tick_interval;
regular_ticks = min_tick:tick_interval:max_tick;
all_ticks = unique([regular_ticks, times(end)]);
set(imgax_exp, 'XTick', all_ticks);

set(imgax_exp, 'FontSize', 7);
xlabel(imgax_exp, 'Time (ms)', 'Color', 'k', 'FontSize', 8);
ylabel(imgax_exp, 'Channel', 'Color', 'k', 'FontSize', 8);

h_xlabel = get(imgax_exp, 'XLabel');
h_ylabel = get(imgax_exp, 'YLabel');
xlabel_pos = get(h_xlabel, 'Position');
ylabel_pos = get(h_ylabel, 'Position');
set(h_xlabel, 'Position', [xlabel_pos(1), xlabel_pos(2) + 3, xlabel_pos(3)]);
set(h_ylabel, 'Position', [ylabel_pos(1) - 25, ylabel_pos(2), ylabel_pos(3)]);

title(h, 'signed -log10(p)', 'FontSize', 8, 'Color', 'k');
set(h, 'FontSize', 7);

plot_title = channel_time_title;
if isempty(plot_title)
    plot_title = default_channel_time_title;
end
if isempty(channel_time_title) && ~any(p_values(:) < alpha_threshold)
    plot_title = sprintf('%s - No Significance', plot_title);
end

annotation('textbox', [0, 0.78, 1, 0.10], ...
           'String', plot_title, ...
           'EdgeColor', 'none', ...
           'HorizontalAlignment', 'center', ...
           'VerticalAlignment', 'middle', ...
           'FontSize', 14, ...
           'Color', 'k', ...
           'Interpreter', 'none');

set(fig_exp, 'PaperPositionMode', 'auto');
set(fig_exp, 'PaperUnits', 'inches');
set(fig_exp, 'InvertHardCopy', 'off');

% Keep PNG at high-resolution raster (upsampled display), but build the
% channel-time SVG from the native channel x time grid to avoid explosive
% vector primitive counts and SVG renderer memory pressure.
save_channel_time_png_and_svg(fig_exp, output_file_png, output_file_svg, ...
    imgax_exp, h_img_exp, times, channels, signed_logp_data, alpha_mask_significance);
close(fig_exp);

fprintf('Saved channel-time plot PNG: %s\n', output_file_png);

% ==================== LIKELIHOOD RATIO PLOT ====================
if has_likelihood_data
    output_file_lr = fullfile(output_dir, sprintf('%s_channel_time_plot_LR.png', test_name));

    LR_strong_threshold = 10;
    alpha_nonsig_lr = 0.4;

    LR_data = log10_likelihood_values;
    LR_data(isinf(LR_data) & LR_data > 0) = 10;
    LR_data(isinf(LR_data) & LR_data < 0) = -10;

    alpha_mask_evidence = ones(size(likelihood_values));
    alpha_mask_evidence(likelihood_values < LR_strong_threshold & likelihood_values > 1/LR_strong_threshold) = alpha_nonsig_lr;

    abs_max_lr = max(abs(LR_data(:)));
    if isempty(abs_max_lr) || abs_max_lr == 0
        abs_max_lr = 1;
    end
    cc_lr = limo_color_images([-abs_max_lr, abs_max_lr]);

    fig_bf = figure('Visible', 'off', 'Position', [100 100 1000 700], 'Color', 'w', 'InvertHardcopy', 'off');
    colormap(fig_bf, cc_lr);

    imgax_bf = axes('Position', [0.08 0.12 0.78 0.68]);
    set(imgax_bf, 'Color', [0.9 0.9 0.9], 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'LineWidth', 1);

    h_img_bf = imagesc(times, 1:size(LR_data,1), LR_data);
    xlim([times(1), times(end)]);
    set(h_img_bf, 'AlphaData', alpha_mask_evidence);
    set(imgax_bf, 'Color', [0.9 0.9 0.9]);
    caxis([-abs_max_lr, abs_max_lr]);

    xlabel('Time (ms)', 'Color', 'k');
    ylabel('Channel', 'Color', 'k');
    set(imgax_bf, 'YDir', 'reverse');
    ylim([1 size(LR_data,1)]);
    grid(imgax_bf, 'off');
    set(imgax_bf, 'XColor', 'k', 'YColor', 'k');
    hold on;

    ylim_vals_bf = get(imgax_bf, 'YLim');
    plot(imgax_bf, [0 0], ylim_vals_bf, 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'LineWidth', 1);

    colormap(imgax_bf, cc_lr);

    h_bf = colorbar(imgax_bf, 'Position', [0.88 0.12 0.04 0.68]);
    title(h_bf, 'log10(LR)', 'FontSize', 10, 'Color', 'k');
    set(h_bf, 'Color', 'k', 'Box', 'on', 'LineWidth', 1, 'TickLength', 0);

    log10_LR_threshold_strong_H1 = log10(LR_strong_threshold);
    log10_LR_threshold_strong_H0 = log10(1/LR_strong_threshold);

    cbar_pos_bf    = get(h_bf, 'Position');
    cbar_limits_bf = get(h_bf, 'Limits');

    y_strong_H1 = NaN;
    y_strong_H0 = NaN;

    if log10_LR_threshold_strong_H1 >= cbar_limits_bf(1) && log10_LR_threshold_strong_H1 <= cbar_limits_bf(2)
        norm_pos = (log10_LR_threshold_strong_H1 - cbar_limits_bf(1)) / (cbar_limits_bf(2) - cbar_limits_bf(1));
        y_strong_H1 = cbar_pos_bf(2) + norm_pos * cbar_pos_bf(4);
    end
    if log10_LR_threshold_strong_H0 >= cbar_limits_bf(1) && log10_LR_threshold_strong_H0 <= cbar_limits_bf(2)
        norm_pos = (log10_LR_threshold_strong_H0 - cbar_limits_bf(1)) / (cbar_limits_bf(2) - cbar_limits_bf(1));
        y_strong_H0 = cbar_pos_bf(2) + norm_pos * cbar_pos_bf(4);
    end

    if isnan(y_strong_H0) && isnan(y_strong_H1)
        y_strong_H0 = cbar_pos_bf(2);
        y_strong_H1 = cbar_pos_bf(2) + cbar_pos_bf(4);
    end

    if ~isnan(y_strong_H0) && ~isnan(y_strong_H1)
        border_inset   = 0.0006;
        overlay_x      = cbar_pos_bf(1) + border_inset;
        overlay_y      = y_strong_H0;
        overlay_width  = cbar_pos_bf(3) - 2*border_inset;
        overlay_height = y_strong_H1 - y_strong_H0;
        annotation('rectangle', [overlay_x, overlay_y, overlay_width, overlay_height], ...
            'FaceColor', [0.9137, 0.9020, 0.9020], ...
            'FaceAlpha', 1 - alpha_nonsig_lr, ...
            'EdgeColor', 'none');
    end
    if ~isnan(y_strong_H1)
        annotation('line', [cbar_pos_bf(1), cbar_pos_bf(1) + cbar_pos_bf(3)], [y_strong_H1, y_strong_H1], ...
            'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);
    end
    if ~isnan(y_strong_H0)
        annotation('line', [cbar_pos_bf(1), cbar_pos_bf(1) + cbar_pos_bf(3)], [y_strong_H0, y_strong_H0], ...
            'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);
    end

    set(imgax_bf, 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'LineWidth', 1, 'TickLength', [0 0]);

    num_channels_bf = size(LR_data, 1);
    chan_labels_bf = {};
    does_not_begin_with_E_bf = false;
    if isfield(LIMO.data, 'chanlocs') && ~isempty(LIMO.data.chanlocs)
        chan_labels_bf = {LIMO.data.chanlocs.labels};
        does_not_begin_with_E_bf = all(cellfun(@(x) isempty(regexp(x, '^E', 'once')), chan_labels_bf));
    end

    if does_not_begin_with_E_bf
        ytick_values_bf = 1:num_channels_bf;
    else
        ytick_values_bf = round(linspace(1, num_channels_bf, 10));
    end
    set(imgax_bf, 'YTick', ytick_values_bf);

    if ~isempty(chan_labels_bf) && length(chan_labels_bf) == num_channels_bf
        set(imgax_bf, 'YTickLabel', chan_labels_bf(ytick_values_bf));
    end

    tick_interval_bf = 250;
    min_tick_bf = ceil(times(1) / tick_interval_bf) * tick_interval_bf;
    max_tick_bf = floor(times(end) / tick_interval_bf) * tick_interval_bf;
    all_ticks_bf = unique([min_tick_bf:tick_interval_bf:max_tick_bf, times(end)]);
    set(imgax_bf, 'XTick', all_ticks_bf);

    set(imgax_bf, 'FontSize', 7);
    xlabel(imgax_bf, 'Time (ms)', 'Color', 'k', 'FontSize', 8);
    ylabel(imgax_bf, 'Channel', 'Color', 'k', 'FontSize', 8);

    h_xlabel_bf = get(imgax_bf, 'XLabel');
    h_ylabel_bf = get(imgax_bf, 'YLabel');
    xlabel_pos_bf = get(h_xlabel_bf, 'Position');
    ylabel_pos_bf = get(h_ylabel_bf, 'Position');
    set(h_xlabel_bf, 'Position', [xlabel_pos_bf(1), xlabel_pos_bf(2) + 3, xlabel_pos_bf(3)]);
    set(h_ylabel_bf, 'Position', [ylabel_pos_bf(1) - 25, ylabel_pos_bf(2), ylabel_pos_bf(3)]);

    title(h_bf, 'log10(LR)', 'FontSize', 8, 'Color', 'k');
    set(h_bf, 'FontSize', 7);

    plot_title_lr = lr_title;
    if isempty(plot_title_lr)
        plot_title_lr = default_lr_title;
    end
    if isempty(lr_title) && ~any(likelihood_values(:) > LR_strong_threshold) && ~any(likelihood_values(:) < 1/LR_strong_threshold)
        plot_title_lr = sprintf('%s - No Strong Evidence', plot_title_lr);
    end

    annotation('textbox', [0, 0.78, 1, 0.10], ...
               'String', plot_title_lr, ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'middle', ...
               'FontSize', 14, ...
               'Color', 'k', ...
               'Interpreter', 'none');

    set(fig_bf, 'PaperPositionMode', 'auto');
    set(fig_bf, 'PaperUnits', 'inches');
    set(fig_bf, 'InvertHardCopy', 'off');

    print(fig_bf, output_file_lr, '-dpng', '-r600');
    close(fig_bf);

    fprintf('Saved LR plot: %s\n', output_file_lr);
else
    fprintf('Skipping LR plot for %s (no likelihood data file)\n', test_name);
end

% ==================== TOPOPLOT GRID ====================
output_file_topo_png = fullfile(output_dir, sprintf('%s_topoplots.png', test_name));
output_file_topo_svg = fullfile(output_dir, sprintf('%s_topoplots.svg', test_name));
max_latency = min(times(end), 2500);
topoplot_step_ms = 100;
latencies = topoplot_step_ms:topoplot_step_ms:max_latency;

if times(1) > 0
    latencies = latencies(latencies >= times(1));
elseif times(1) < 0
    latencies = [times(1):topoplot_step_ms:0, topoplot_step_ms:topoplot_step_ms:max_latency];
    latencies = latencies(latencies >= times(1) & latencies <= times(end));
end

if isempty(latencies)
    latencies = linspace(times(1), times(end), 12);
end

time_indices = arrayfun(@(x) find(abs(times - x) == min(abs(times - x)), 1), latencies);

abs_max_topo = max(abs(signed_logp_data(:)));
if isempty(abs_max_topo) || abs_max_topo == 0
    abs_max_topo = 1;
end
cc_topo = limo_color_images([-abs_max_topo, abs_max_topo]);
topoplot_bg_rgb = [232 230 232] / 255;
cc_topo = set_colormap_mid_color(cc_topo, topoplot_bg_rgb, 2);
colorbar_limits_topo = [-abs_max_topo, abs_max_topo];

num_latencies = length(latencies);
[topo_fig_position, topo_ax_positions, topo_cbar_position] = ...
    compute_topoplot_layout(num_latencies, topoplot_layout_type);

fig_topo = figure('Visible', 'off', 'Position', topo_fig_position, 'Color', 'w', 'InvertHardcopy', 'off');
ax_handles_topo = gobjects(1, num_latencies);

for i = 1:num_latencies
    ax_seed = axes('Parent', fig_topo, 'Position', topo_ax_positions(i, :));
    axes(ax_seed);

    data_timepoint = signed_logp_data(:, time_indices(i));
    p_timepoint    = p_values(:, time_indices(i));

    data_masked = data_timepoint;
    data_masked(p_timepoint >= alpha_threshold) = 0;

    topoplot(data_masked, LIMO.data.chanlocs, ...
        'maplimits', colorbar_limits_topo, ...
        'electrodes', 'off', ...
        'style', 'both', ...
        'shading', 'interp', ...
        'numcontour', 6);
    ax_topo = gca;
    if isgraphics(ax_seed) && ax_seed ~= ax_topo
        delete(ax_seed);
    end
    set(ax_topo, 'Units', 'normalized', ...
        'Position', topo_ax_positions(i, :), ...
        'ActivePositionProperty', 'position', ...
        'PositionConstraint', 'innerposition', ...
        'Color', topoplot_bg_rgb);
    ax_handles_topo(i) = ax_topo;

    text(ax_topo, 0.5, -0.12, sprintf('%d ms', round(latencies(i))), ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 8, ...
        'FontWeight', 'bold', ...
        'Color', 'k', ...
        'Clipping', 'off');
end

colormap(fig_topo, cc_topo);

for i = 1:length(ax_handles_topo)
    set(ax_handles_topo(i), 'CLim', colorbar_limits_topo);
end

ax_ref_topo = axes('Position', [0.001 0.001 0.001 0.001], 'Visible', 'off');
imagesc(ax_ref_topo, linspace(colorbar_limits_topo(1), colorbar_limits_topo(2), 100)');
set(ax_ref_topo, 'CLim', colorbar_limits_topo);
colormap(ax_ref_topo, cc_topo);

h_topo = colorbar(ax_ref_topo, 'Position', topo_cbar_position);
tick_values_topo = [colorbar_limits_topo(1), colorbar_limits_topo(1)/2, 0, colorbar_limits_topo(2)/2, colorbar_limits_topo(2)];
tick_labels_topo = arrayfun(@(x) sprintf('%.1f', x), tick_values_topo, 'UniformOutput', false);
set(h_topo, 'Ticks', tick_values_topo, 'TickLabels', tick_labels_topo, 'Color', 'k', 'Box', 'on');
title(h_topo, 'signed -log10(p)', 'FontSize', 10, 'Color', 'k');
set(ax_ref_topo, 'XTick', [], 'YTick', [], 'Box', 'off');

plot_title_topo = topoplot_title;
if isempty(plot_title_topo)
    plot_title_topo = default_topoplot_title;
end
if isempty(topoplot_title) && ~any(p_values(:) < 0.05)
    plot_title_topo = sprintf('%s - No Significance', plot_title_topo);
end

sgtitle(plot_title_topo, 'FontSize', 14, 'Color', 'k', 'Interpreter', 'none');

set(fig_topo, 'PaperPositionMode', 'auto');
set(fig_topo, 'PaperUnits', 'inches');
set(fig_topo, 'InvertHardCopy', 'off');

save_png_and_svg(fig_topo, output_file_topo_png, output_file_topo_svg);
close(fig_topo);

fprintf('Saved topoplot PNG: %s\n', output_file_topo_png);

fprintf('\nAll plots completed! Saved to: %s\n', output_dir);

end

function display_name = derive_display_contrast_name(input_folder, folder_name)
[label1, label2] = read_contrast_labels_from_metadata(input_folder);
if ~isempty(label1) && ~isempty(label2)
    display_name = sprintf('%s vs %s', ...
        normalize_title_label(label1), normalize_title_label(label2));
    return;
end

prefix = 'limo_second_level_';
if startsWith(folder_name, prefix)
    contrast_key = folder_name((length(prefix) + 1):end);
else
    contrast_key = folder_name;
end

parts = strsplit(contrast_key, '_vs_');
if numel(parts) == 2 && ~isempty(parts{1}) && ~isempty(parts{2})
    display_name = sprintf('%s vs %s', ...
        normalize_title_label(parts{1}), normalize_title_label(parts{2}));
else
    display_name = normalize_title_label(contrast_key);
end
end

function [label1, label2] = read_contrast_labels_from_metadata(input_folder)
label1 = '';
label2 = '';
meta_file = fullfile(input_folder, 'contrast_metadata.txt');
if ~exist(meta_file, 'file')
    return;
end

fid = fopen(meta_file, 'r');
if fid == -1
    return;
end

cleanup_obj = onCleanup(@() fclose(fid));
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line)
        continue;
    end
    line = strtrim(line);
    if isempty(line)
        continue;
    end

    parts = strsplit(line, sprintf('\t'));
    if numel(parts) < 2
        continue;
    end

    key = lower(strtrim(parts{1}));
    value = strtrim(strjoin(parts(2:end), sprintf('\t')));
    if strcmp(key, 'label1')
        label1 = value;
    elseif strcmp(key, 'label2')
        label2 = value;
    end
end
end

function out = normalize_title_label(in)
out = strtrim(char(in));
if isempty(out)
    return;
end

if strcmpi(out, 'correct_rejections')
    out = 'Correct_rejection';
    return;
end

out(1) = upper(out(1));
end

function layout_type = normalize_topoplot_layout_type(layout_value)
layout_type = lower(strtrim(char(layout_value)));
if isempty(layout_type)
    layout_type = 'grid';
    return;
end

switch layout_type
    case {'grid', 'square', 'square_grid', 'square-grid'}
        layout_type = 'grid';
    case {'zigzag', 'zigzag_line', 'zigzag-line', 'staggered', 'staggered_line', 'staggered-line'}
        layout_type = 'zigzag';
    case {'line', 'horizontal', 'row', 'single_row', 'single-row', 'one_row', 'one-row'}
        layout_type = 'line';
    otherwise
        warning('Unknown topoplot layout type "%s". Falling back to "grid".', layout_type);
        layout_type = 'grid';
end
end

function [fig_position, ax_positions, colorbar_position] = compute_topoplot_layout(num_plots, layout_type)
if num_plots < 1
    num_plots = 1;
end

plot_region = [0.04 0.12 0.84 0.74];

switch layout_type
    case 'zigzag'
        % Compute in pixels so each topoplot is truly square on-screen.
        % This avoids tiny circles on wide canvases and preserves 90-degree
        % zigzag turns in rendered output (dx == dy in pixels).
        fig_width = max(2600, 900 + 95 * num_plots);
        fig_height = 1000;

        left_px   = plot_region(1) * fig_width;
        bottom_px = plot_region(2) * fig_height;
        region_w_px = plot_region(3) * fig_width;
        region_h_px = plot_region(4) * fig_height;

        ax_positions = zeros(num_plots, 4);
        if num_plots == 1
            d_px = min(region_w_px, region_h_px) * 0.75;
            x_center_px = left_px + region_w_px / 2;
            y_center_px = bottom_px + region_h_px / 2;
            ax_positions(1, :) = [ ...
                (x_center_px - d_px/2) / fig_width, ...
                (y_center_px - d_px/2) / fig_height, ...
                d_px / fig_width, ...
                d_px / fig_height];
        else
            near_touch_ratio = 0.98;
            diameter_per_step = near_touch_ratio * sqrt(2);
            step_px_by_width = region_w_px / ((num_plots - 1) + diameter_per_step);
            step_px_by_height = region_h_px / (1 + diameter_per_step);
            step_px = min(step_px_by_width, step_px_by_height);
            d_px = diameter_per_step * step_px;
            y_step_px = step_px; % true 90-degree zigzag

            x_first_px = left_px + d_px / 2;
            y_mid_px = bottom_px + region_h_px / 2;
            y_high_px = y_mid_px + y_step_px / 2;
            y_low_px  = y_mid_px - y_step_px / 2;

            for idx = 1:num_plots
                x_center_px = x_first_px + (idx - 1) * step_px;
                if mod(idx, 2) == 1
                    y_center_px = y_high_px;
                else
                    y_center_px = y_low_px;
                end
                ax_positions(idx, :) = [ ...
                    (x_center_px - d_px/2) / fig_width, ...
                    (y_center_px - d_px/2) / fig_height, ...
                    d_px / fig_width, ...
                    d_px / fig_height];
            end
        end
        colorbar_position = [0.92 0.16 0.03 0.68];

    case 'line'
        % Single horizontal row with even spacing.
        fig_width = max(2600, 900 + 130 * num_plots);
        fig_height = 700;

        left_px   = plot_region(1) * fig_width;
        bottom_px = plot_region(2) * fig_height;
        region_w_px = plot_region(3) * fig_width;
        region_h_px = plot_region(4) * fig_height;

        ax_positions = zeros(num_plots, 4);
        if num_plots == 1
            d_px = min(region_w_px, region_h_px) * 0.72;
            x_center_px = left_px + region_w_px / 2;
            y_center_px = bottom_px + region_h_px * 0.55;
            ax_positions(1, :) = [ ...
                (x_center_px - d_px/2) / fig_width, ...
                (y_center_px - d_px/2) / fig_height, ...
                d_px / fig_width, ...
                d_px / fig_height];
        else
            gap_ratio = 0.08;
            d_px_by_width = region_w_px / (num_plots + (num_plots - 1) * gap_ratio);
            d_px_by_height = region_h_px * 0.78;
            d_px = min(d_px_by_width, d_px_by_height);
            gap_px = gap_ratio * d_px;
            total_w_px = num_plots * d_px + (num_plots - 1) * gap_px;

            x_first_px = left_px + (region_w_px - total_w_px) / 2 + d_px / 2;
            y_center_px = bottom_px + region_h_px * 0.55;

            for idx = 1:num_plots
                x_center_px = x_first_px + (idx - 1) * (d_px + gap_px);
                ax_positions(idx, :) = [ ...
                    (x_center_px - d_px/2) / fig_width, ...
                    (y_center_px - d_px/2) / fig_height, ...
                    d_px / fig_width, ...
                    d_px / fig_height];
            end
        end
        colorbar_position = [0.92 0.16 0.03 0.68];

    otherwise
        % Near-square grid (ncols ~= nrows) for compact overview layout.
        ncols = ceil(sqrt(num_plots));
        nrows = ceil(num_plots / ncols);
        cell_pixels = 220;
        fig_width = max(1200, 200 + ncols * cell_pixels);
        fig_height = max(700, 200 + nrows * cell_pixels);

        cell_width = plot_region(3) / ncols;
        cell_height = plot_region(4) / nrows;
        plot_size = min(cell_width, cell_height) * 0.90;

        ax_positions = zeros(num_plots, 4);
        for idx = 1:num_plots
            row = ceil(idx / ncols);
            col = mod(idx - 1, ncols) + 1;
            x_center = plot_region(1) + (col - 0.5) * cell_width;
            y_center = plot_region(2) + plot_region(4) - (row - 0.5) * cell_height;
            ax_positions(idx, :) = [x_center - plot_size/2, y_center - plot_size/2, plot_size, plot_size];
        end
        colorbar_position = [0.92 0.1 0.03 0.8];
end

fig_position = [100 50 fig_width fig_height];
end

function cmap_out = set_colormap_mid_color(cmap_in, mid_rgb, half_width)
cmap_out = cmap_in;
if isempty(cmap_in)
    return;
end

n = size(cmap_in, 1);
mid = round((n + 1) / 2);
idx_start = max(1, mid - half_width);
idx_end = min(n, mid + half_width);
cmap_out(idx_start:idx_end, :) = repmat(mid_rgb, idx_end - idx_start + 1, 1);
end

function save_png_and_svg(fig_handle, png_path, svg_path)
print(fig_handle, png_path, '-dpng', '-r600');

try
    set(fig_handle, 'Renderer', 'painters');
    print(fig_handle, svg_path, '-dsvg');
    fprintf('Saved SVG: %s\n', svg_path);
catch ME
    warning('SVG export failed for %s: %s', svg_path, ME.message);
end
end

function save_channel_time_png_and_svg(fig_handle, png_path, svg_path, ax_handle, img_handle, x_vals, y_vals, c_data, alpha_data)
% Save PNG from the original raster image, then save SVG with only the
% channel-time image converted to vector surface primitives.
print(fig_handle, png_path, '-dpng', '-r600');

try
    x_edges = image_grid_edges(x_vals);
    y_edges = image_grid_edges(y_vals);
    [x_edge_grid, y_edge_grid] = meshgrid(x_edges, y_edges);

    c_surface = pad_replicate_post(c_data);
    a_surface = pad_replicate_post(alpha_data);

    hold_state = ishold(ax_handle);
    hold(ax_handle, 'on');
    if isgraphics(img_handle)
        delete(img_handle);
    end

    h_surface = surface(ax_handle, ...
        x_edge_grid, y_edge_grid, zeros(size(x_edge_grid)), c_surface, ...
        'EdgeColor', 'none', ...
        'FaceColor', 'flat', ...
        'FaceAlpha', 'flat', ...
        'AlphaData', a_surface, ...
        'AlphaDataMapping', 'none');
    uistack(h_surface, 'bottom');
    if ~hold_state
        hold(ax_handle, 'off');
    end

    set(fig_handle, 'Renderer', 'painters');
    print(fig_handle, svg_path, '-dsvg');
    fprintf('Saved SVG: %s\n', svg_path);
catch ME
    warning('Channel-time SVG vector export failed for %s: %s', svg_path, ME.message);
end
end

function edges = image_grid_edges(vals)
vals = vals(:)';
if numel(vals) < 2
    if isempty(vals)
        edges = [0 1];
    else
        edges = [vals(1) - 0.5, vals(1) + 0.5];
    end
    return;
end

edges = zeros(1, numel(vals) + 1);
edges(2:end-1) = (vals(1:end-1) + vals(2:end)) / 2;
edges(1) = vals(1) - (vals(2) - vals(1)) / 2;
edges(end) = vals(end) + (vals(end) - vals(end-1)) / 2;
end

function out = pad_replicate_post(in)
[nr, nc] = size(in);
out = zeros(nr + 1, nc + 1, class(in));
out(1:nr, 1:nc) = in;
out(nr + 1, 1:nc) = in(nr, :);
out(1:nr, nc + 1) = in(:, nc);
out(nr + 1, nc + 1) = in(nr, nc);
end
