function hpc_limo_channel_time_plots(input_folder, output_dir)
% HPC_LIMO_CHANNEL_TIME_PLOTS - Generate channel-time and topoplot figures
%
% Usage:
%   hpc_limo_channel_time_plots(input_folder, output_dir)
%
% Inputs:
%   input_folder - Path to a limo_second_level_<output_tag> folder containing
%                  LIMO.mat and paired_samples_ttest_parameter_*.mat
%   output_dir   - Directory to save the output .png files
%
% Outputs:
%   - <output_dir>/<test_name>_channel_time_plot.png      Signed -log10(p) channel-time plot
%   - <output_dir>/<test_name>_channel_time_plot_LR.png   Likelihood ratio plot (if data exists)
%   - <output_dir>/<test_name>_topoplots.png              Topoplot grid at multiple time points
%
% The test_name is derived from the input_folder name.
% Example: limo_second_level_test_hits_vs_test_misses

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
output_file  = fullfile(output_dir, sprintf('%s_channel_time_plot.png', test_name));
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

fig_exp = figure('Visible', 'off', 'Position', [100 100 1000 700], 'Color', 'w', 'InvertHardcopy', 'off');
colormap(fig_exp, cc);

imgax_exp = axes('Position', [0.08 0.12 0.78 0.68]);
set(imgax_exp, 'Color', [0.9 0.9 0.9], 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'LineWidth', 1);

signed_logp_data = real(signed_logp_data);
h_img_exp = imagesc(times, 1:size(signed_logp_data,1), signed_logp_data);
xlim([times(1), times(end)]);

set(h_img_exp, 'AlphaData', alpha_mask_significance);
set(imgax_exp, 'Color', [0.9 0.9 0.9]);
caxis([-abs_max_val, abs_max_val]);

xlabel('Time (ms)', 'Color', 'k');
ylabel('Channel', 'Color', 'k');

set(imgax_exp, 'YDir', 'reverse');
ylim([1 size(signed_logp_data,1)]);
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

num_channels = size(signed_logp_data, 1);
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

tick_interval = 250;
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

plot_title = sprintf('%s With TFCE Correction', test_name);
if ~any(p_values(:) < alpha_threshold)
    plot_title = sprintf('%s - No Significance', plot_title);
end

annotation('textbox', [0, 0.78, 1, 0.10], ...
           'String', plot_title, ...
           'EdgeColor', 'none', ...
           'HorizontalAlignment', 'center', ...
           'VerticalAlignment', 'middle', ...
           'FontSize', 14, ...
           'Color', 'k');

set(fig_exp, 'PaperPositionMode', 'auto');
set(fig_exp, 'PaperUnits', 'inches');
set(fig_exp, 'InvertHardCopy', 'off');

print(fig_exp, output_file, '-dpng', '-r600');
close(fig_exp);

fprintf('Saved channel-time plot: %s\n', output_file);

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

    plot_title_lr = sprintf('%s Likelihood Ratios', test_name);
    if ~any(likelihood_values(:) > LR_strong_threshold) && ~any(likelihood_values(:) < 1/LR_strong_threshold)
        plot_title_lr = sprintf('%s - No Strong Evidence', plot_title_lr);
    end

    annotation('textbox', [0, 0.78, 1, 0.10], ...
               'String', plot_title_lr, ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'middle', ...
               'FontSize', 14, ...
               'Color', 'k');

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
output_file_topo = fullfile(output_dir, sprintf('%s_topoplots.png', test_name));
max_latency = min(times(end), 2500);
latencies = 200:200:max_latency;

if times(1) > 0
    latencies = latencies(latencies >= times(1));
elseif times(1) < 0
    latencies = [times(1):200:0, 200:200:max_latency];
    latencies = latencies(latencies >= times(1) & latencies <= times(end));
end

if isempty(latencies)
    latencies = linspace(times(1), times(end), 12);
end

time_indices = arrayfun(@(x) find(abs(times - x) == min(abs(times - x)), 1), latencies);

abs_max_topo = max(abs(t_values(:)));
if isempty(abs_max_topo) || abs_max_topo == 0
    abs_max_topo = 1;
end
cc_topo = limo_color_images([-abs_max_topo, abs_max_topo]);
colorbar_limits_topo = [-abs_max_topo, abs_max_topo];

num_latencies = length(latencies);
ncols = 4;
nrows = ceil(num_latencies / ncols);

fig_topo = figure('Visible', 'off', 'Position', [100 50 1200 300*nrows], 'Color', 'w', 'InvertHardcopy', 'off');
ax_handles_topo = [];

for i = 1:num_latencies
    ax_topo = subplot(nrows, ncols, i);
    ax_handles_topo(end+1) = ax_topo;

    data_timepoint = t_values(:, time_indices(i));
    p_timepoint    = p_values(:, time_indices(i));

    data_masked = data_timepoint;
    data_masked(p_timepoint >= 0.05) = 0;

    topoplot(data_masked, LIMO.data.chanlocs, ...
        'maplimits', colorbar_limits_topo, ...
        'electrodes', 'off', ...
        'style', 'both', ...
        'shading', 'interp', ...
        'numcontour', 6);

    title(sprintf('%d ms', round(latencies(i))), 'FontSize', 8, 'Color', 'k');

    pos = get(ax_topo, 'Position');
    center_x = pos(1) + pos(3)/2;
    center_y = pos(2) + pos(4)/2;
    new_width  = pos(3) * 1.1;
    new_height = pos(4) * 1.1;
    new_x = center_x - new_width/2;
    new_y = center_y - new_height/2;
    set(ax_topo, 'Position', [new_x, new_y, new_width, new_height]);
end

colormap(fig_topo, cc_topo);

for i = 1:length(ax_handles_topo)
    set(ax_handles_topo(i), 'CLim', colorbar_limits_topo);
end

ax_ref_topo = axes('Position', [0.001 0.001 0.001 0.001], 'Visible', 'off');
imagesc(ax_ref_topo, linspace(colorbar_limits_topo(1), colorbar_limits_topo(2), 100)');
set(ax_ref_topo, 'CLim', colorbar_limits_topo);
colormap(ax_ref_topo, cc_topo);

h_topo = colorbar(ax_ref_topo, 'Position', [0.92 0.1 0.03 0.8]);
tick_values_topo = [colorbar_limits_topo(1), colorbar_limits_topo(1)/2, 0, colorbar_limits_topo(2)/2, colorbar_limits_topo(2)];
tick_labels_topo = arrayfun(@(x) sprintf('%.1f', x), tick_values_topo, 'UniformOutput', false);
set(h_topo, 'Ticks', tick_values_topo, 'TickLabels', tick_labels_topo, 'Color', 'k', 'Box', 'on');
title(h_topo, 't-values', 'FontSize', 10, 'Color', 'k');
set(ax_ref_topo, 'XTick', [], 'YTick', [], 'Box', 'off');

plot_title_topo = sprintf('%s Topoplots (TFCE-corrected)', test_name);
if ~any(p_values(:) < 0.05)
    plot_title_topo = sprintf('%s - No Significance', plot_title_topo);
end

sgtitle(plot_title_topo, 'FontSize', 14, 'Color', 'k');

for i = 1:length(ax_handles_topo)
    pos = get(ax_handles_topo(i), 'Position');
    set(ax_handles_topo(i), 'Position', [pos(1), pos(2)*0.925, pos(3), pos(4)]);
end

set(fig_topo, 'PaperPositionMode', 'auto');
set(fig_topo, 'PaperUnits', 'inches');
set(fig_topo, 'InvertHardCopy', 'off');

print(fig_topo, output_file_topo, '-dpng', '-r600');
close(fig_topo);

fprintf('Saved topoplot grid: %s\n', output_file_topo);

fprintf('\nAll plots completed! Saved to: %s\n', output_dir);

end
