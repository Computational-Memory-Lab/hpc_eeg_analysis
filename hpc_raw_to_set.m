function hpc_raw_to_set(input_folder, subject_filter)
% HPC_RAW_TO_SET - Convert EGI .raw files to .set and import behavioral events
%
% Usage:
%   hpc_raw_to_set(input_folder)
%   hpc_raw_to_set(input_folder, subject_filter)
%
% Inputs:
%   input_folder - Folder containing subject .raw files plus required
%                  session .log files and EEG .eeglog files.
%   subject_filter - Optional numeric subject ID. When provided, process
%                    only that subject.
%
% Expected session log format:
%   - Contains a tab-delimited header row with the required columns.
%   - Required columns: timestamp, trigger_code, event_label
%   - Optional columns: block, trial, word1_id/w1id, word2_id/w2id, target, response,
%     accuracy/acc, rt/rt_ms, include_in_analysis
%
% Outputs (created under the parent directory of input_folder):
%   - initial_set/<ID>.set
%   - behavioral_set/processed_<ID>.set
%   - behavioral_set/behavioral_<ID>.log
%   - behavioral_set/EEGevents_<ID>.txt
%   - behavioral_set/alignment_parameters_S<ID>.mat

if nargin < 1 || nargin > 2 || ~(ischar(input_folder) || (isstring(input_folder) && isscalar(input_folder)))
    error('Usage: hpc_raw_to_set(input_folder [, subject_filter])');
end

if isstring(input_folder)
    input_folder = char(input_folder);
end
input_folder = normalize_input_folder(input_folder);

if nargin < 2
    subject_filter = [];
else
    subject_filter = parse_optional_subject_filter(subject_filter);
end

if ~exist(input_folder, 'dir')
    error('Input folder not found: %s', input_folder);
end

if ~exist('pop_readegi', 'file')
    error('EEGLAB EGI reader (pop_readegi) not found on MATLAB path.');
end
if ~exist('pop_loadset', 'file')
    error('EEGLAB (pop_loadset) not found on MATLAB path.');
end
if ~exist('pop_importevent', 'file')
    error('EEGLAB (pop_importevent) not found on MATLAB path.');
end

[parent_folder, ~, ~] = fileparts(input_folder);
if isempty(parent_folder)
    parent_folder = pwd;
end

initial_set_folder = fullfile(parent_folder, 'initial_set');
behavioral_set_folder = fullfile(parent_folder, 'behavioral_set');

if ~exist(initial_set_folder, 'dir')
    mkdir(initial_set_folder);
end
if ~exist(behavioral_set_folder, 'dir')
    mkdir(behavioral_set_folder);
end

raw_files = dir(fullfile(input_folder, '*.raw'));
if isempty(raw_files)
    error('No .raw files found in: %s', input_folder);
end

fprintf('\n==============================================\n');
fprintf('  HPC RAW -> SET PIPELINE\n');
fprintf('==============================================\n');
fprintf('Input folder:         %s\n', input_folder);
fprintf('Initial set folder:   %s\n', initial_set_folder);
fprintf('Behavioral set folder:%s\n', behavioral_set_folder);
fprintf('Raw files found:      %d\n', numel(raw_files));
fprintf('==============================================\n\n');

success_count = 0;
failure_messages = {};
matched_subject_count = 0;

for i = 1:numel(raw_files)
    raw_name = raw_files(i).name;
    raw_path = fullfile(raw_files(i).folder, raw_name);
    [~, raw_stem, ~] = fileparts(raw_name);

    token = regexp(raw_stem, '^(\d+)$', 'tokens', 'once');
    if isempty(token)
        token = regexp(raw_stem, '(\d+)', 'tokens', 'once');
    end

    if isempty(token)
        warning('Skipping %s: could not extract subject ID from filename.', raw_name);
        continue;
    end

    subject_id = str2double(token{1});

    if ~isempty(subject_filter) && subject_id ~= subject_filter
        continue;
    end
    matched_subject_count = matched_subject_count + 1;

    fprintf('----------------------------------------------\n');
    fprintf('Subject %d (%d/%d)\n', subject_id, i, numel(raw_files));
    fprintf('Raw file: %s\n', raw_path);

    try
        initial_set_file = fullfile(initial_set_folder, sprintf('%d.set', subject_id));

        if exist(initial_set_file, 'file')
            fprintf('Initial set already exists, reusing: %s\n', initial_set_file);
        else
            fprintf('Converting raw -> initial set...\n');
            EEG = pop_readegi(raw_path, [], [], 'auto');
            EEG.setname = sprintf('%draw', subject_id);
            EEG = eeg_checkset(EEG);
            EEG = pop_saveset(EEG, 'filename', sprintf('%d.set', subject_id), 'filepath', initial_set_folder); %#ok<NASGU>
            fprintf('Saved initial set: %s\n', initial_set_file);
        end

        session_file = find_subject_session_log(input_folder, subject_id);
        eeglog_file = find_subject_eeglog(input_folder, subject_id);

        output_set_file = fullfile(behavioral_set_folder, sprintf('processed_%d.set', subject_id));
        if exist(output_set_file, 'file')
            fprintf('Processed set already exists, skipping behavioral import: %s\n', output_set_file);
        else
            process_subject_behavioral_alignment(subject_id, session_file, initial_set_file, eeglog_file, behavioral_set_folder);
        end

        success_count = success_count + 1;
        fprintf('Subject %d complete.\n\n', subject_id);

    catch ME
        msg = sprintf('Subject %d failed: %s', subject_id, ME.message);
        warning('%s', msg);
        failure_messages{end+1} = msg; %#ok<AGROW>
        fprintf('\n');
    end
end

if ~isempty(subject_filter) && matched_subject_count == 0
    error('Subject %d not found in raw files under: %s', subject_filter, input_folder);
end

fprintf('==============================================\n');
fprintf('  RAW -> SET PIPELINE SUMMARY\n');
fprintf('==============================================\n');
fprintf('Successful subjects: %d\n', success_count);
fprintf('Failed subjects:     %d\n', numel(failure_messages));
fprintf('Initial set folder:  %s\n', initial_set_folder);
fprintf('Behavioral set folder:%s\n', behavioral_set_folder);

if ~isempty(failure_messages)
    for i = 1:numel(failure_messages)
        fprintf('  - %s\n', failure_messages{i});
    end
    error('hpc_raw_to_set completed with failures. See summary above.');
else
    fprintf('All subjects processed successfully.\n');
end
fprintf('==============================================\n\n');

end

function session_file = find_subject_session_log(input_folder, subject_id)
all_logs = dir(fullfile(input_folder, '*.log'));
if ~isempty(all_logs)
    is_eeglog = endsWith({all_logs.name}, '.eeglog');
    all_logs = all_logs(~is_eeglog);
end

candidates = filter_files_by_subject_id(all_logs, subject_id);

if isempty(candidates) && numel(all_logs) == 1
    % Backward-compatible fallback for single-subject folders.
    candidates = all_logs;
end

session_file = require_single_match(candidates, input_folder, subject_id, 'session .log');
end

function eeglog_file = find_subject_eeglog(input_folder, subject_id)
candidates_all = dir(fullfile(input_folder, '*.eeglog'));
candidates = filter_files_by_subject_id(candidates_all, subject_id);

if isempty(candidates) && numel(candidates_all) == 1
    % Backward-compatible fallback for single-subject folders.
    candidates = candidates_all;
end

eeglog_file = require_single_match(candidates, input_folder, subject_id, '.eeglog');
end

function filtered = filter_files_by_subject_id(files, subject_id)
filtered = files([]);
if isempty(files)
    return;
end

sid = num2str(subject_id);
sid_pattern = ['(^|[^0-9])' regexptranslate('escape', sid) '([^0-9]|$)'];

for i = 1:numel(files)
    [~, stem, ~] = fileparts(files(i).name);
    if ~isempty(regexp(stem, sid_pattern, 'once'))
        filtered(end + 1) = files(i); %#ok<AGROW>
    end
end
end

function matched_file = require_single_match(candidates, input_folder, subject_id, file_type)
if isempty(candidates)
    error('No %s file found for subject %d in %s', file_type, subject_id, input_folder);
end

if numel(candidates) > 1
    names = strjoin({candidates.name}, ', ');
    error('Multiple %s files matched subject %d: %s', file_type, subject_id, names);
end

if isfield(candidates, 'folder') && ~isempty(candidates(1).folder)
    matched_file = fullfile(candidates(1).folder, candidates(1).name);
else
    matched_file = fullfile(input_folder, candidates(1).name);
end
end

function process_subject_behavioral_alignment(subject_id, session_file, eeg_file, eeglog_file, output_folder)
if ~exist(session_file, 'file')
    error('Session file not found: %s', session_file);
end
if ~exist(eeg_file, 'file')
    error('EEG set file not found: %s', eeg_file);
end
if ~exist(eeglog_file, 'file')
    error('EEGlog file not found: %s', eeglog_file);
end

[~, session_filename, session_ext] = fileparts(session_file);
session_filename = [session_filename session_ext];
[eeg_path, eeg_filename, eeg_ext] = fileparts(eeg_file);
eeg_filename = [eeg_filename eeg_ext];

fprintf('\n==============================================\n');
fprintf('  PROCESSING SUBJECT %d\n', subject_id);
fprintf('==============================================\n');
fprintf('Session file: %s\n', session_file);
fprintf('EEG set:      %s\n', eeg_file);
fprintf('EEGlog file:  %s\n\n', eeglog_file);

%% Part 1: Read behavioral events directly from session log
fprintf('STEP 1: Loading behavioral events from session log...\n');
events = read_session_events(session_file);
fprintf('  Events retained for import: %d\n', numel(events));

% Event label summary
labels = {events.event_label};
unique_labels = unique(labels);
for i = 1:numel(unique_labels)
    c = sum(strcmp(labels, unique_labels{i}));
    fprintf('    %-30s %d\n', unique_labels{i}, c);
end

behavioral_file = fullfile(output_folder, sprintf('behavioral_%d.log', subject_id));
fid = fopen(behavioral_file, 'w');
if fid == -1
    error('Could not open output file for writing: %s', behavioral_file);
end
fprintf(fid, '# Behavioral events for subject %d\n', subject_id);
fprintf(fid, '# Source session file: %s\n', session_filename);
fprintf(fid, 'timestamp\tblock\ttrial\tword1_id\tword2_id\ttarget\tresponse\taccuracy\trt\ttrigger_code\tevent_label\n');
for i = 1:numel(events)
    fprintf(fid, '%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%s\n', ...
        events(i).timestamp, events(i).block, events(i).trial, events(i).word1_id, ...
        events(i).word2_id, events(i).target, events(i).response, events(i).accuracy, ...
        events(i).rt, events(i).trigger_code, events(i).event_label);
end
fclose(fid);
fprintf('  Saved: %s\n\n', behavioral_file);

%% Part 2: EEG Time Alignment
fprintf('STEP 2: Aligning to EEG time...\n');
EEG = pop_loadset('filename', eeg_filename, 'filepath', eeg_path);

nstime = [];
for ii = 1:length(EEG.event)
    nstime = [nstime; EEG.event(1, ii).latency]; %#ok<AGROW>
end

diffnstime_sample = diff(nstime);
diffnstime_ms = diff(nstime) * 2;

dbtime = find(diffnstime_sample < 8);
dbpulses = length(dbtime);
fprintf('  Duplicate pulses found: %d\n', dbpulses);

correctednstime_sample = nstime(setdiff(1:length(nstime), dbtime + 1));
correcteddiffnstime_sample = diff(correctednstime_sample);
correcteddiffnstime_ms = diffnstime_ms(setdiff(1:length(diffnstime_ms), dbtime)); %#ok<NASGU>

fid = fopen(eeglog_file, 'r');
if fid == -1
    error('Could not open eeglog file: %s', eeglog_file);
end
fs = repmat('%s', 1, 4);
data = textscan(fid, fs, 'delimiter', '\t');
fclose(fid);

eeglogtime = str2double(data{1});
class = data{3};

upin = strcmp('UP', class);
tupin = strcmp('TRAIN_UP', class);
upindex = (upin + tupin);
uptime = eeglogtime(find(upindex == 1));
diffeeglogtime_sample = round(diff(uptime) / 2);

fprintf('  Aligning pulse sequences...\n');
fprintf('  PyEPL pulses: %d, EEG pulses: %d\n', length(diffeeglogtime_sample), length(correcteddiffnstime_sample));

max_dev = [0 0];
alignindex = 1;

if length(diffeeglogtime_sample) > length(correcteddiffnstime_sample)
    check = 0;
    ii = 1;
    while check == 0 && ii < length(diffeeglogtime_sample)
        end_idx = length(correcteddiffnstime_sample) + ii - 2;
        if end_idx > length(diffeeglogtime_sample)
            break;
        end
        aligndiff = diffeeglogtime_sample(ii:end_idx) - correcteddiffnstime_sample(2:end);
        maxdiff = max(aligndiff);
        mindiff = min(aligndiff);

        if maxdiff > 16 || mindiff < -16
            ii = ii + 1;
        else
            check = 1;
            max_dev = [maxdiff mindiff];
            alignindex = ii;
        end
    end

    end_idx = length(correctednstime_sample) + ii - 2;
    if end_idx > length(uptime)
        end_idx = length(uptime);
    end
    trimmed_uptime = uptime(ii:end_idx);
    trimmed_nstime = correctednstime_sample(2:(end_idx - ii + 2)) * 2;

elseif length(diffeeglogtime_sample) < length(correcteddiffnstime_sample)
    check = 0;
    ii = 1;
    while check == 0 && ii < length(diffeeglogtime_sample)
        end_idx = length(diffeeglogtime_sample) - ii + 2;
        if end_idx > length(correcteddiffnstime_sample)
            break;
        end
        aligndiff = diffeeglogtime_sample(ii:end) - correcteddiffnstime_sample(2:end_idx);
        maxdiff = max(aligndiff);
        mindiff = min(aligndiff);

        if maxdiff > 16 || mindiff < -16
            ii = ii + 1;
        else
            check = 1;
            max_dev = [maxdiff mindiff];
            alignindex = ii;
        end
    end

    trimmed_uptime = uptime(ii:end);
    trimmed_nstime = correctednstime_sample(2:(length(uptime(ii:end)) + 1)) * 2;

else
    ii = 1;
    aligndiff = diffeeglogtime_sample - correcteddiffnstime_sample;
    maxdiff = max(aligndiff);
    mindiff = min(aligndiff);
    if maxdiff < 16 && mindiff > -16
        max_dev = [maxdiff mindiff];
        alignindex = ii;
        trimmed_uptime = uptime(ii:end);
        trimmed_nstime = correctednstime_sample(2:(length(uptime(ii:end)) + 1)) * 2;
    end
end

if isequal(max_dev, [0 0])
    error(['Alignment failed! Could not find matching pulse pattern between PyEPL and EEG.\n' ...
           'PyEPL pulses: %d, EEG pulses: %d\n' ...
           'This usually means:\n' ...
           '  1. Wrong eeg.eeglog file selected\n' ...
           '  2. EEG and PyEPL recordings are from different sessions\n' ...
           '  3. Pulse timing tolerance (+/-16 samples) is too strict'], ...
           length(diffeeglogtime_sample), length(correcteddiffnstime_sample));
end

P = polyfit(trimmed_uptime - uptime(ii), trimmed_nstime, 1);
fprintf('  Slope: %.16f, Intercept: %.16f\n', P);
fprintf('  Max deviation: [%d, %d] samples\n', max_dev);
fprintf('  Alignment index: %d\n', alignindex);

alignment_time = now;
alignment_time_string = datestr(alignment_time); %#ok<NASGU>
alignment_file = fullfile(output_folder, sprintf('alignment_parameters_S%d.mat', subject_id));
save(alignment_file, 'P', 'max_dev', 'alignment_time*', 'dbpulses');
fprintf('  Saved: %s\n', alignment_file);

if P(1) <= 0.999
    error('Alignment quality poor (slope = %.6f). Check eeglog!', P(1));
end

for i = 1:numel(events)
    events(i).latency = (events(i).timestamp - uptime(ii)) * P(1) + P(2);
end

eeg_events_file = fullfile(output_folder, sprintf('EEGevents_%d.txt', subject_id));
fid = fopen(eeg_events_file, 'w');
if fid == -1
    error('Could not open output file for writing: %s', eeg_events_file);
end
for i = 1:numel(events)
    fprintf(fid, '%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%s\t%.10f\n', ...
        events(i).timestamp, events(i).block, events(i).trial, events(i).word1_id, ...
        events(i).word2_id, events(i).target, events(i).response, events(i).accuracy, ...
        events(i).rt, events(i).trigger_code, events(i).event_label, events(i).latency);
end
fclose(fid);
fprintf('  Saved: %s\n\n', eeg_events_file);

%% Part 3: Import events into EEGLAB
fprintf('STEP 3: Importing events into EEGLAB...\n');

EEG = pop_importevent(EEG, 'append', 'no', 'event', eeg_events_file, ...
    'fields', {'unixtime', 'block', 'trial', 'W1id', 'W2id', ...
               'target', 'response', 'accuracy', 'RT', 'type', 'event_label', 'latency'}, ...
    'timeunit', 0.001);

fprintf('  Events imported: %d\n', length(EEG.event));

output_eeg_file = sprintf('processed_%d.set', subject_id);
EEG = pop_saveset(EEG, 'filename', output_eeg_file, 'filepath', output_folder); %#ok<NASGU>
fprintf('  Saved: %s\n', fullfile(output_folder, output_eeg_file));

fprintf('\n==============================================\n');
fprintf('  SUBJECT %d COMPLETE\n', subject_id);
fprintf('==============================================\n');
fprintf('  Behavioral log:   %s\n', behavioral_file);
fprintf('  EEG events file:  %s\n', eeg_events_file);
fprintf('  EEG dataset:      %s\n', fullfile(output_folder, output_eeg_file));
fprintf('  Alignment params: %s\n', alignment_file);
fprintf('==============================================\n\n');

end

function events = read_session_events(session_file)
lines = read_text_lines(session_file);
if isempty(lines)
    error('Session log is empty: %s', session_file);
end

[header_idx, header_parts, header_norm] = detect_header_line(lines, session_file);

idx_timestamp = resolve_required_column(header_norm, {'timestamp', 'unixtime', 'unix_time'}, 'timestamp');
idx_block = resolve_optional_column(header_norm, {'block'});
idx_trial = resolve_optional_column(header_norm, {'trial'});
idx_trigger = resolve_required_column(header_norm, {'trigger_code', 'trigger', 'type', 'event_code'}, 'trigger_code');
idx_event_label = resolve_required_column(header_norm, {'event_label', 'trial_type', 'condition'}, 'event_label');

idx_word1 = resolve_optional_column(header_norm, {'word1_id', 'w1id', 'word1', 'word1id'});
idx_word2 = resolve_optional_column(header_norm, {'word2_id', 'w2id', 'word2', 'word2id'});
idx_target = resolve_optional_column(header_norm, {'target'});
idx_response = resolve_optional_column(header_norm, {'response'});
idx_accuracy = resolve_optional_column(header_norm, {'accuracy', 'acc'});
idx_rt = resolve_optional_column(header_norm, {'rt', 'rt_ms'});
idx_include = resolve_optional_column(header_norm, {'include_in_analysis', 'include'});

events = struct('timestamp', {}, 'block', {}, 'trial', {}, 'word1_id', {}, ...
                'word2_id', {}, 'target', {}, 'response', {}, 'accuracy', {}, ...
                'rt', {}, 'trigger_code', {}, 'event_label', {}, 'latency', {});

for line_idx = (header_idx + 1):numel(lines)
    raw_line = lines{line_idx};
    if isempty(strtrim(raw_line))
        continue;
    end

    parts = regexp(raw_line, '\t', 'split');
    if numel(parts) < numel(header_parts)
        parts(end+1:numel(header_parts)) = {''};
    elseif numel(parts) > numel(header_parts)
        parts = parts(1:numel(header_parts));
    end

    event_label = strtrim(parts{idx_event_label});
    if isempty(event_label)
        continue;
    end

    if ~isempty(idx_include)
        include_value = parse_include_flag(parts{idx_include});
        if ~include_value
            continue;
        end
    end

    timestamp = parse_numeric(parts{idx_timestamp});
    trigger_code = parse_numeric(parts{idx_trigger});

    if isnan(timestamp)
        warning('Skipping line %d in %s: invalid timestamp.', line_idx, session_file);
        continue;
    end
    if isnan(trigger_code)
        warning('Skipping line %d in %s: invalid trigger_code.', line_idx, session_file);
        continue;
    end

    e = struct();
    e.timestamp = timestamp;
    e.block = parse_optional_numeric(parts, idx_block, -1);
    e.trial = parse_optional_numeric(parts, idx_trial, -1);
    e.word1_id = parse_optional_numeric(parts, idx_word1, -1);
    e.word2_id = parse_optional_numeric(parts, idx_word2, -1);
    e.target = parse_optional_numeric(parts, idx_target, -1);
    e.response = parse_optional_numeric(parts, idx_response, -1);
    e.accuracy = parse_optional_numeric(parts, idx_accuracy, -1);
    e.rt = parse_optional_numeric(parts, idx_rt, -1);
    e.trigger_code = trigger_code;
    e.event_label = sanitize_event_label(event_label);
    e.latency = NaN;

    events(end+1) = e; %#ok<AGROW>
end

if isempty(events)
    error('No analyzable events were found in session log: %s', session_file);
end

[~, sort_idx] = sort([events.timestamp]);
events = events(sort_idx);
end

function [header_idx, header_parts, header_norm] = detect_header_line(lines, session_file)
header_idx = [];
header_parts = {};
header_norm = {};

for i = 1:numel(lines)
    line = strtrim(lines{i});
    if isempty(line)
        continue;
    end

    parts = regexp(line, '\t', 'split');
    norm_parts = cellfun(@normalize_col_name, parts, 'UniformOutput', false);

    has_timestamp = has_any_alias(norm_parts, {'timestamp', 'unixtime', 'unix_time'});
    has_trigger = has_any_alias(norm_parts, {'trigger_code', 'trigger', 'type', 'event_code'});
    has_event_label = has_any_alias(norm_parts, {'event_label', 'trial_type', 'condition'});

    if has_timestamp && has_trigger && has_event_label
        header_idx = i;
        header_parts = parts;
        header_norm = norm_parts;
        return;
    end
end

error(['Could not locate a valid session-log header in %s. ' ...
       'Expected a tab-delimited header containing timestamp, trigger_code, and event_label (or aliases).'], ...
      session_file);
end

function tf = has_any_alias(header_norm, aliases)
tf = false;
for i = 1:numel(aliases)
    candidate = normalize_col_name(aliases{i});
    if any(strcmp(header_norm, candidate))
        tf = true;
        return;
    end
end
end

function lines = read_text_lines(file_path)
fid = fopen(file_path, 'r');
if fid == -1
    error('Could not open file: %s', file_path);
end
lines = {};
while ~feof(fid)
    tline = fgetl(fid);
    if ischar(tline)
        lines{end+1} = tline; %#ok<AGROW>
    end
end
fclose(fid);
end

function idx = resolve_required_column(header_norm, aliases, display_name)
idx = resolve_optional_column(header_norm, aliases);
if isempty(idx)
    error('Required column "%s" not found in session log header.', display_name);
end
end

function idx = resolve_optional_column(header_norm, aliases)
idx = [];
for i = 1:numel(aliases)
    candidate = normalize_col_name(aliases{i});
    hit = find(strcmp(header_norm, candidate), 1, 'first');
    if ~isempty(hit)
        idx = hit;
        return;
    end
end
end

function out = normalize_col_name(in)
out = lower(strtrim(in));
out = strrep(out, char(65279), '');
out = strrep(out, ' ', '_');
out = strrep(out, '-', '_');
end

function n = parse_numeric(raw)
if isnumeric(raw)
    if isscalar(raw)
        n = raw;
    else
        n = NaN;
    end
    return;
end

if isstring(raw)
    raw = char(raw);
elseif ~ischar(raw)
    n = NaN;
    return;
end

raw = strtrim(raw);
n = str2double(raw);
if ~isnan(n)
    return;
end

match = regexp(raw, '-?\d+(\.\d+)?', 'match', 'once');
if isempty(match)
    n = NaN;
else
    n = str2double(match);
end
end

function n = parse_numeric_or_default(raw, default_value)
n = parse_numeric(raw);
if isnan(n)
    n = default_value;
end
end

function n = parse_optional_numeric(parts, idx, default_value)
if isempty(idx)
    n = default_value;
else
    n = parse_numeric_or_default(parts{idx}, default_value);
end
end

function include_value = parse_include_flag(raw)
value = lower(strtrim(raw));
if isempty(value)
    include_value = true;
elseif any(strcmp(value, {'1', 'true', 'yes', 'y'}))
    include_value = true;
elseif any(strcmp(value, {'0', 'false', 'no', 'n'}))
    include_value = false;
else
    include_value = true;
end
end

function label = sanitize_event_label(raw)
label = strtrim(raw);
label = strrep(label, sprintf('\t'), ' ');
label = strrep(label, sprintf('\n'), ' ');
label = strrep(label, sprintf('\r'), ' ');
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

function out = normalize_input_folder(path_in)
out = strtrim(char(path_in));
if isempty(out)
    return;
end

while numel(out) > 1 && (out(end) == '/' || out(end) == '\')
    out = out(1:end-1);
end
end
