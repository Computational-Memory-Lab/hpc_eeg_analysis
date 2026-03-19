function summary = summarize_session_log(session_log_path)
%SUMMARIZE_SESSION_LOG Print behavioral summary from a ####_session.log file.
%   SUMMARY = SUMMARIZE_SESSION_LOG(SESSION_LOG_PATH) reads the tab-delimited
%   session log and prints timing and performance metrics to the console.
%
%   Metrics are computed from rows where include_in_analysis == 1 for:
%   - list-level durations
%   - ITEM/ASSOC RT, accuracy, and SDT-style rates
%
%   Total session time is computed from the full log (including practice rows).

if ~(ischar(session_log_path) || (isstring(session_log_path) && isscalar(session_log_path)))
    error('session_log_path must be a character vector or string scalar.');
end
session_log_path = char(session_log_path);

if ~isfile(session_log_path)
    error('File not found: %s', session_log_path);
end

T = readtable( ...
    session_log_path, ...
    'FileType', 'text', ...
    'Delimiter', '\t', ...
    'VariableNamingRule', 'preserve');

required_cols = {'timestamp', 'list_id', 'phase', 'accuracy', 'rt', 'event_label', 'include_in_analysis'};
missing_cols = required_cols(~ismember(required_cols, T.Properties.VariableNames));
if ~isempty(missing_cols)
    error('Missing required columns: %s', strjoin(missing_cols, ', '));
end

timestamp = to_numeric(T.timestamp);
rt = to_numeric(T.rt);
accuracy = to_numeric(T.accuracy);
include_in_analysis = to_numeric(T.include_in_analysis) == 1;
phase = upper(strtrim(string(T.phase)));
phase(ismissing(phase)) = "";
list_id = strtrim(string(T.list_id));
list_id(ismissing(list_id)) = "";
event_label = lower(strtrim(string(T.event_label)));
event_label(ismissing(event_label)) = "";

% Session timing across the full file.
valid_start = ~isnan(timestamp);
if ~any(valid_start)
    error('No valid numeric timestamps found in %s', session_log_path);
end
session_start = min(timestamp(valid_start));
valid_end = ~isnan(timestamp) & ~isnan(rt);
if any(valid_end)
    session_end = max(timestamp(valid_end) + rt(valid_end));
else
    session_end = max(timestamp(valid_start));
end
total_session_ms = max(0, session_end - session_start);

% List durations from analysis-included task rows only.
analysis_rows = include_in_analysis & ismember(phase, ["STUDY", "ITEM", "ASSOC"]) & (list_id ~= "");
analysis_lists = unique(list_id(analysis_rows), 'stable');
n_lists = numel(analysis_lists);
list_duration_ms = nan(n_lists, 1);
list_type = strings(n_lists, 1);

for i = 1:n_lists
    this_list = analysis_lists(i);
    mask = analysis_rows & (list_id == this_list);

    ts = timestamp(mask);
    te = timestamp(mask) + rt(mask);
    ts = ts(~isnan(ts));
    te = te(~isnan(te));
    if ~isempty(ts) && ~isempty(te)
        list_duration_ms(i) = max(te) - min(ts);
    end

    ph = phase(mask);
    if any(ph == "ITEM")
        list_type(i) = "ITEM";
    elseif any(ph == "ASSOC")
        list_type(i) = "ASSOC";
    else
        list_type(i) = "UNKNOWN";
    end
end

avg_all_list_ms = mean(list_duration_ms, 'omitnan');
item_list_mask = (list_type == "ITEM") & ~isnan(list_duration_ms);
assoc_list_mask = (list_type == "ASSOC") & ~isnan(list_duration_ms);
n_item_lists = sum(item_list_mask);
n_assoc_lists = sum(assoc_list_mask);
avg_item_list_ms = mean(list_duration_ms(item_list_mask), 'omitnan');
avg_assoc_list_ms = mean(list_duration_ms(assoc_list_mask), 'omitnan');

item_stats = compute_phase_stats("ITEM", include_in_analysis, phase, rt, accuracy, event_label);
assoc_stats = compute_phase_stats("ASSOC", include_in_analysis, phase, rt, accuracy, event_label);

participant_id = extract_participant_id(session_log_path);

summary = struct();
summary.participant_id = participant_id;
summary.total_session_ms = total_session_ms;
summary.n_lists = n_lists;
summary.avg_all_list_ms = avg_all_list_ms;
summary.n_item_lists = n_item_lists;
summary.n_assoc_lists = n_assoc_lists;
summary.avg_item_list_ms = avg_item_list_ms;
summary.avg_assoc_list_ms = avg_assoc_list_ms;
summary.item = item_stats;
summary.assoc = assoc_stats;

fprintf('Participant %s:\n', participant_id);
fprintf('Total session time: %s\n', format_hours_minutes(total_session_ms));
fprintf('Average time per all lists (n=%d): %s\n\n', n_lists, format_minutes_seconds(avg_all_list_ms));

fprintf('Item:\n');
fprintf('Average time per list (n=%d): %s\n', n_item_lists, format_minutes_seconds(avg_item_list_ms));
fprintf('Average RT = %.2f ms\n', item_stats.avg_rt_ms);
fprintf('Accuracy = %.2f%%\n', item_stats.accuracy_pct);
fprintf('Hit rate: %.2f%%\n', item_stats.hit_rate_pct);
fprintf('Miss rate: %.2f%%\n', item_stats.miss_rate_pct);
fprintf('FA rate: %.2f%%\n', item_stats.fa_rate_pct);
fprintf('CR rate: %.2f%%\n', item_stats.cr_rate_pct);
fprintf('d'' = %.3f, c = %.3f\n\n', item_stats.d_prime, item_stats.criterion_c);

fprintf('Assoc:\n');
fprintf('Average time per list (n=%d): %s\n', n_assoc_lists, format_minutes_seconds(avg_assoc_list_ms));
fprintf('Average RT = %.2f ms\n', assoc_stats.avg_rt_ms);
fprintf('Accuracy = %.2f%%\n', assoc_stats.accuracy_pct);
fprintf('Hit rate: %.2f%%\n', assoc_stats.hit_rate_pct);
fprintf('Miss rate: %.2f%%\n', assoc_stats.miss_rate_pct);
fprintf('FA rate: %.2f%%\n', assoc_stats.fa_rate_pct);
fprintf('CR rate: %.2f%%\n', assoc_stats.cr_rate_pct);
fprintf('d'' = %.3f, c = %.3f\n', assoc_stats.d_prime, assoc_stats.criterion_c);

end

function stats = compute_phase_stats(phase_name, include_in_analysis, phase, rt, accuracy, event_label)
mask = include_in_analysis & (phase == phase_name);
n_trials = sum(mask);

stats = struct();
stats.n_trials = n_trials;
stats.avg_rt_ms = mean(rt(mask), 'omitnan');
stats.accuracy_pct = 100 * mean(accuracy(mask), 'omitnan');

stats.hits_n = sum(mask & contains(event_label, 'test_hits'));
stats.misses_n = sum(mask & contains(event_label, 'test_misses'));
stats.fas_n = sum(mask & contains(event_label, 'false_alarms'));
stats.crs_n = sum(mask & contains(event_label, 'correct_rejections'));

stats.n_signal = stats.hits_n + stats.misses_n;
stats.n_noise = stats.fas_n + stats.crs_n;

if stats.n_signal > 0
    hit_rate = stats.hits_n / stats.n_signal;
    miss_rate = stats.misses_n / stats.n_signal;
else
    hit_rate = NaN;
    miss_rate = NaN;
end

if stats.n_noise > 0
    fa_rate = stats.fas_n / stats.n_noise;
    cr_rate = stats.crs_n / stats.n_noise;
else
    fa_rate = NaN;
    cr_rate = NaN;
end

stats.hit_rate = hit_rate;
stats.miss_rate = miss_rate;
stats.fa_rate = fa_rate;
stats.cr_rate = cr_rate;
stats.hit_rate_pct = 100 * hit_rate;
stats.miss_rate_pct = 100 * miss_rate;
stats.fa_rate_pct = 100 * fa_rate;
stats.cr_rate_pct = 100 * cr_rate;

% Keep overall-trial proportions for backward compatibility.
if n_trials > 0
    stats.hits_pct = 100 * stats.hits_n / n_trials;
    stats.misses_pct = 100 * stats.misses_n / n_trials;
    stats.fas_pct = 100 * stats.fas_n / n_trials;
    stats.crs_pct = 100 * stats.crs_n / n_trials;
else
    stats.hits_pct = NaN;
    stats.misses_pct = NaN;
    stats.fas_pct = NaN;
    stats.crs_pct = NaN;
end

if ~isnan(hit_rate) && ~isnan(fa_rate)
    adj_hit = adjust_rate_for_sdt(hit_rate, stats.n_signal);
    adj_fa = adjust_rate_for_sdt(fa_rate, stats.n_noise);
    z_hit = norminv_local(adj_hit);
    z_fa = norminv_local(adj_fa);
    stats.d_prime = z_hit - z_fa;
    stats.criterion_c = -0.5 * (z_hit + z_fa);
else
    stats.d_prime = NaN;
    stats.criterion_c = NaN;
end
end

function out = format_hours_minutes(ms)
if isnan(ms)
    out = 'N/A';
    return;
end
total_seconds = max(0, round(ms / 1000));
hours = floor(total_seconds / 3600);
minutes = floor(mod(total_seconds, 3600) / 60);
if hours > 0
    out = sprintf('%dh %02dm', hours, minutes);
else
    out = sprintf('%dm', minutes);
end
end

function out = format_minutes_seconds(ms)
if isnan(ms)
    out = 'N/A';
    return;
end
total_seconds = max(0, round(ms / 1000));
minutes = floor(total_seconds / 60);
seconds = mod(total_seconds, 60);
out = sprintf('%dm %02ds', minutes, seconds);
end

function participant_id = extract_participant_id(session_log_path)
[~, base, ~] = fileparts(session_log_path);

token = regexp(base, '^(\d+)_session$', 'tokens', 'once');
if isempty(token)
    token = regexp(base, '(\d+)', 'tokens', 'once');
end

if isempty(token)
    participant_id = base;
else
    participant_id = token{1};
end
end

function out = to_numeric(v)
if isnumeric(v)
    out = double(v);
elseif islogical(v)
    out = double(v);
elseif iscell(v)
    out = str2double(string(v));
else
    out = str2double(string(v));
end
end

function p_adj = adjust_rate_for_sdt(p, n)
if isnan(p) || isnan(n) || n <= 0
    p_adj = NaN;
    return;
end
lower = 0.5 / n;
upper = 1 - lower;
p_adj = min(max(p, lower), upper);
end

function z = norminv_local(p)
z = -sqrt(2) * erfcinv(2 * p);
end
