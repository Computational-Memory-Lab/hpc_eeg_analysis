function hpc_raw_to_set(input_folder)
% HPC_RAW_TO_SET - Convert EGI .raw files to .set and import behavioral events
%
% Usage:
%   hpc_raw_to_set(input_folder)
%
% Inputs:
%   input_folder - Folder containing subject .raw files plus required
%                  session .log files and EEG .eeglog files.
%
% Outputs (created under the parent directory of input_folder):
%   - initial_set/<ID>.set                        Initial EEGLAB dataset converted from .raw
%   - initial_set/behavioral_set/processed_<ID>.set
%                                                 EEG dataset with imported behavioral events
%   - initial_set/behavioral_set/behavioral_<ID>.log
%   - initial_set/behavioral_set/EEGevents_<ID>.txt
%   - initial_set/behavioral_set/alignment_parameters_S<ID>.mat
%
% Notes:
%   - Subject IDs are parsed from raw filenames (e.g., 1007.raw -> 1007).
%   - This function is non-interactive (no file picker usage).

if nargin ~= 1 || ~(ischar(input_folder) || (isstring(input_folder) && isscalar(input_folder)))
    error('Usage: hpc_raw_to_set(input_folder)');
end

if isstring(input_folder)
    input_folder = char(input_folder);
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
behavioral_set_folder = fullfile(initial_set_folder, 'behavioral_set');

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
fprintf('Input folder:        %s\n', input_folder);
fprintf('Initial set folder:  %s\n', initial_set_folder);
fprintf('Behavioral set folder:%s\n', behavioral_set_folder);
fprintf('Raw files found:     %d\n', numel(raw_files));
fprintf('==============================================\n\n');

success_count = 0;
failure_messages = {};

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
pattern = fullfile(input_folder, sprintf('*%d*.log', subject_id));
candidates = dir(pattern);

if ~isempty(candidates)
    is_eeglog = endsWith({candidates.name}, '.eeglog');
    candidates = candidates(~is_eeglog);
end

if isempty(candidates)
    all_logs = dir(fullfile(input_folder, '*.log'));
    is_eeglog = endsWith({all_logs.name}, '.eeglog');
    candidates = all_logs(~is_eeglog);
end

session_file = require_single_match(candidates, input_folder, subject_id, 'session .log');
end

function eeglog_file = find_subject_eeglog(input_folder, subject_id)
pattern = fullfile(input_folder, sprintf('*%d*.eeglog', subject_id));
candidates = dir(pattern);

if isempty(candidates)
    candidates = dir(fullfile(input_folder, '*.eeglog'));
end

eeglog_file = require_single_match(candidates, input_folder, subject_id, '.eeglog');
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
% Process one subject: parse behavior, align timing, import events, and save .set

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

%% Part 1: Process Behavioral Data
fprintf('STEP 1: Processing behavioral data...\n');

fid = fopen(session_file);
if fid == -1
    error('Could not open session file: %s', session_file);
end
fgetl(fid);
fgetl(fid);

lines = {};
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line)
        lines{end+1} = line; %#ok<AGROW>
    end
end
fclose(fid);

study_trials = [];
assoc_trials = [];
item_trials = [];

num_study = 0;
num_assoc_intact = 0;
num_assoc_recomb = 0;
num_item_old = 0;
num_item_new = 0;
num_pilot_excluded = 0;

for i = 1:length(lines)
    parts = strsplit(lines{i}, '\t');
    if length(parts) < 3
        continue;
    end

    if strcmp(parts{3}, 'B') || strcmp(parts{3}, 'E')
        continue;
    elseif strcmp(parts{3}, 'LIST')
        continue;
    elseif startsWith(parts{3}, 'P')
        num_pilot_excluded = num_pilot_excluded + 1;
        continue;
    elseif length(parts) >= 13 && strcmp(parts{5}, 'ASSOC')
        timestamp = str2double(parts{1});
        block = str2double(parts{3});
        trial = str2double(parts{4});
        word1_id = str2double(parts{7});
        word2_id = str2double(parts{9});
        target = str2double(parts{10});
        response = str2double(parts{11});
        accuracy = str2double(parts{12});
        rt = str2double(parts{13});

        if target == 1
            type_code = 21;
            num_assoc_intact = num_assoc_intact + 1;
        else
            type_code = 22;
            num_assoc_recomb = num_assoc_recomb + 1;
        end

        assoc_trials(end+1,:) = [timestamp block trial word1_id word2_id target response accuracy rt type_code]; %#ok<AGROW>

    elseif length(parts) >= 11 && strcmp(parts{5}, 'ITEM')
        timestamp = str2double(parts{1});
        block = str2double(parts{3});
        trial = str2double(parts{4});
        word_id = str2double(parts{7});
        target = str2double(parts{8});
        response = str2double(parts{9});
        accuracy = str2double(parts{10});
        rt = str2double(parts{11});

        if target == 1
            type_code = 31;
            num_item_old = num_item_old + 1;
        else
            type_code = 32;
            num_item_new = num_item_new + 1;
        end

        item_trials(end+1,:) = [timestamp block trial word_id -1 target response accuracy rt type_code]; %#ok<AGROW>

    elseif length(parts) >= 9
        timestamp = str2double(parts{1});
        block = str2double(parts{3});
        trial = str2double(parts{4});
        word1_id = str2double(parts{6});
        word2_id = str2double(parts{8});
        ipi = str2double(parts{9});

        study_trials(end+1,:) = [timestamp block trial word1_id word2_id -1 -1 -1 ipi -2]; %#ok<AGROW>
        num_study = num_study + 1;
    end
end

if num_pilot_excluded > 0
    fprintf('  Pilot events excluded: %d\n', num_pilot_excluded);
end
fprintf('  Study trials: %d\n', num_study);
fprintf('  ASSOC trials: %d (Intact: %d, Recombined: %d)\n', num_assoc_intact + num_assoc_recomb, num_assoc_intact, num_assoc_recomb);
fprintf('  ITEM trials: %d (Old: %d, New: %d)\n\n', num_item_old + num_item_new, num_item_old, num_item_new);

for i = 1:size(assoc_trials, 1)
    if assoc_trials(i, 10) == 21
        w1_id = assoc_trials(i, 4);
        accuracy = assoc_trials(i, 8);

        study_idx = find(study_trials(:, 4) == w1_id);
        if ~isempty(study_idx)
            study_trials(study_idx(1), 8) = accuracy;
            study_trials(study_idx(1), 10) = 11;
        end
    end
end

unmatched = find(study_trials(:, 10) == -2);
study_trials(unmatched, 10) = 12;

alleeg = sortrows([study_trials; assoc_trials; item_trials], 1);

behavioral_file = fullfile(output_folder, sprintf('behavioral_%d.log', subject_id));
fid = fopen(behavioral_file, 'w');
if fid == -1
    error('Could not open output file for writing: %s', behavioral_file);
end

fprintf(fid, '# Behavioral events for subject %d\n', subject_id);
fprintf(fid, '# Source session file: %s\n', session_filename);
fprintf(fid, '# Columns: timestamp, block, trial, pad1, pad2, word1_id, word2_id, ipi/target, accuracy, event_type\n');

for i = 1:size(alleeg, 1)
    code = alleeg(i, 10);
    switch code
        case 11
            label = 'Study_Intact';
        case 12
            label = 'Study_Recombined';
        case 21
            label = 'Test_Intact';
        case 22
            label = 'Test_Recombined';
        case 31
            label = 'Test_Old';
        case 32
            label = 'Test_New';
        otherwise
            label = 'Unknown';
    end

    fprintf(fid, '%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%s\n', alleeg(i, 1:9), label);
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
            check = 0;
            ii = ii + 1;
        elseif maxdiff < 16 && mindiff > -16
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
            check = 0;
            ii = ii + 1;
        elseif maxdiff < 16 && mindiff > -16
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

newlat = (alleeg(:, 1) - uptime(ii)) * P(1) + P(2);
newevents = [alleeg newlat];

eeg_events_file = fullfile(output_folder, sprintf('EEGevents_%d.txt', subject_id));
dlmwrite(eeg_events_file, newevents, 'delimiter', '\t', 'precision', 10);
fprintf('  Saved: %s\n\n', eeg_events_file);

%% Part 3: Import Events into EEGLAB
fprintf('STEP 3: Importing events into EEGLAB...\n');

EEG = pop_importevent(EEG, 'append', 'no', 'event', eeg_events_file, ...
    'fields', {'unixtime', 'block', 'trial', 'W1id', 'W2id', ...
               'target', 'response', 'accuracy', 'RT', 'type', 'latency'}, ...
    'timeunit', 0.001);

fprintf('  Events imported: %d\n', length(EEG.event));

for i = 1:length(EEG.event)
    switch EEG.event(i).type
        case 11
            EEG.event(i).event_label = 'Study_Intact';
        case 12
            EEG.event(i).event_label = 'Study_Recombined';
        case 21
            EEG.event(i).event_label = 'Test_Intact';
        case 22
            EEG.event(i).event_label = 'Test_Recombined';
        case 31
            EEG.event(i).event_label = 'Test_Old';
        case 32
            EEG.event(i).event_label = 'Test_New';
        otherwise
            EEG.event(i).event_label = 'Unknown';
    end
end

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
