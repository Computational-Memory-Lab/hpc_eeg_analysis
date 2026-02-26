function hpc_limo_first_level(input_folder, condition_order)
% HPC_LIMO_FIRST_LEVEL - Run LIMO 1st level analysis on epoched EEG data
%
% Usage:
%   hpc_limo_first_level(input_folder)
%   hpc_limo_first_level(input_folder, condition_order)
%
% Inputs:
%   input_folder     - Path to the epoch folder containing *_epoch.set files
%   condition_order  - Optional condition order used by std_makedesign.
%                      Accepts cell array or CSV string.
%                      Example: 'Study_hits,Study_misses,Test_hits,Test_misses'
%
% Outputs:
%   - <parent_of_input_folder>/limo_first_level/
%   - condition_order.mat / condition_order.txt (label -> parameter mapping)

% Start timing
tic;
start_time = datetime('now');
fprintf('\n========================================\n');
fprintf('LIMO 1ST LEVEL ANALYSIS\n');
fprintf('Started at: %s\n', start_time);
fprintf('========================================\n\n');

% Parse optional condition order input
user_supplied_condition_order = false;
if nargin >= 2 && ~isempty(condition_order)
    user_supplied_condition_order = true;
    condition_order = parse_condition_order(condition_order);
else
    condition_order = {};
end
input_folder = normalize_input_folder(input_folder);

% Force headless MATLAB
set(0, 'DefaultFigureVisible', 'off');
java.lang.System.setProperty('java.awt.headless', 'true');

% Add EEGLAB to path
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1');
fprintf('Added EEGLAB to path\n');

% Define directories
pipeline_root = fileparts(input_folder);
if isempty(pipeline_root)
    pipeline_root = pwd;
end
output_dir = fullfile(pipeline_root, 'limo_first_level');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
cd(output_dir);
fprintf('Input folder:  %s\n', input_folder);
fprintf('Output folder: %s\n', output_dir);

% Start EEGLAB safely
try
    eeglab nogui;
catch ME
    disp('EEGLAB failed to start:');
    disp(getReport(ME));
    return;
end

% Disable EEGLAB online checks
setpref('Internet', 'EeglabServerCheck', 'off');
setpref('Internet', 'EeglabUpdateServer', 'off');
setpref('Internet', 'EeglabPluginlistURL', 'off');

disp('LIMO will manage parallel pool automatically...');

% Define STUDY file path
study_file = fullfile(output_dir, 'PairAsso_Epoched.study');

% Load or create STUDY
if exist(study_file, 'file')
    fprintf('Loading existing STUDY file: %s\n', study_file);
    [STUDY, ALLEEG] = pop_loadstudy('filename', 'PairAsso_Epoched.study', 'filepath', output_dir);
else
    fprintf('Creating new STUDY from epoched files in: %s\n', input_folder);

    files_in_dir = dir(fullfile(input_folder, '*_epoch.set'));

    epoched_files = [];
    for i = 1:length(files_in_dir)
        fname = files_in_dir(i).name;
        if ~startsWith(fname, 'sub-') && ~contains(fname, '_ses-') && ~contains(fname, '_design')
            epoched_files = [epoched_files; files_in_dir(i)]; %#ok<AGROW>
        end
    end

    fprintf('Found %d epoched files\n', length(epoched_files));
    if isempty(epoched_files)
        error('No epoched files found in: %s', input_folder);
    end

    ALLEEG = [];
    study_commands = cell(length(epoched_files), 1);

    inferred_condition_order = {};

    for i = 1:length(epoched_files)
        fprintf('Loading file %d/%d: %s\n', i, length(epoched_files), epoched_files(i).name);
        EEG = pop_loadset('filename', epoched_files(i).name, 'filepath', epoched_files(i).folder);

        [~, basename, ~] = fileparts(epoched_files(i).name);
        tokens = regexp(basename, '^(\d+)_epoch$', 'tokens');
        if ~isempty(tokens)
            unique_subject_id = tokens{1}{1};
        else
            error('Failed to parse filename: %s (expected format: <ID>_epoch.set)', epoched_files(i).name);
        end

        EEG.subject = unique_subject_id;
        EEG.condition = '';

        % Infer condition order from trial_type field when not provided
        if ~user_supplied_condition_order
            inferred_condition_order = append_trial_types_in_order(inferred_condition_order, EEG);
        end

        study_commands{i} = {'index', i, 'subject', unique_subject_id};
        fprintf('  Dataset %d: Subject %s\n', i, unique_subject_id);

        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, i); %#ok<ASGLU>
    end

    if ~user_supplied_condition_order
        condition_order = inferred_condition_order;
    end

    [STUDY, ALLEEG] = std_editset([], ALLEEG, 'name', 'PairAsso_Epoched', ...
        'task', 'Associative Memory', ...
        'filename', 'PairAsso_Epoched.study', ...
        'filepath', output_dir, ...
        'commands', study_commands, ...
        'updatedat', 'on');
end

% If condition order was not provided and we loaded existing study, try loading saved order.
if isempty(condition_order)
    condition_order = load_saved_condition_order(output_dir);
end

if isempty(condition_order)
    % Last resort: infer from ALLEEG currently in memory.
    inferred = {};
    for i = 1:length(ALLEEG)
        inferred = append_trial_types_in_order(inferred, ALLEEG(i));
    end
    condition_order = inferred;
end

if isempty(condition_order)
    error(['Could not determine condition order. Provide condition_order explicitly, ' ...
           'or ensure trial_type exists in epoched datasets.']);
end

fprintf('Condition order for std_makedesign:\n');
for i = 1:numel(condition_order)
    fprintf('  %d -> %s\n', i, condition_order{i});
end

% Define / refresh design using the selected condition order
fprintf('Setting up design with %d condition(s)...\n', numel(condition_order));
STUDY = std_makedesign(STUDY, ALLEEG, 1, 'name', 'All_Conditions', ...
    'delfiles', 'off', ...
    'defaultdesign', 'off', ...
    'variable1', 'trial_type', ...
    'values1', condition_order, ...
    'vartype1', 'categorical', ...
    'subjselect', {ALLEEG.subject});

% Save condition order mapping
save(fullfile(output_dir, 'condition_order.mat'), 'condition_order');
write_condition_order_text(fullfile(output_dir, 'condition_order.txt'), condition_order);

% Save STUDY
[STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, 'filename', 'PairAsso_Epoched.study', 'filepath', output_dir);
fprintf('STUDY saved to: %s\n', study_file);

% Precompute ERP measures
fprintf('\n========================================\n');
fprintf('Precomputing ERP Measures\n');
fprintf('========================================\n');
try
    [STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG, {}, ...
        'design', 1, ...
        'erp', 'on', ...
        'recompute', 'on');
    fprintf('ERP measures computed successfully\n');
catch ME
    fprintf('ERROR: ERP precomputation failed:\n');
    disp(getReport(ME));
    return;
end

% Run LIMO 1st level analysis
try
    fprintf('\n========================================\n');
    fprintf('Running LIMO 1st Level Analysis\n');
    fprintf('========================================\n');
    STUDY = pop_limo(STUDY, ALLEEG, 'method', 'WLS', 'measure', 'daterp', ...
        'erase', 'on', 'splitreg', 'off', 'interaction', 'off');
    fprintf('LIMO analysis completed successfully\n');
catch ME
    fprintf('ERROR: 1st level LIMO failed:\n');
    disp(getReport(ME));
end

disp('Analysis complete.');

% Cleanup parallel pool
p = gcp('nocreate');
if ~isempty(p)
    disp('Shutting down parallel pool...');
    delete(p);
end

% End timing
end_time = datetime('now');
elapsed = toc;
fprintf('\n========================================\n');
fprintf('TIMING SUMMARY\n');
fprintf('========================================\n');
fprintf('Start time:  %s\n', start_time);
fprintf('End time:    %s\n', end_time);
fprintf('Elapsed:     %.2f seconds (%.2f minutes)\n', elapsed, elapsed / 60);
fprintf('========================================\n');
disp('Done.');

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

function out = parse_condition_order(value)
if ischar(value) || (isstring(value) && isscalar(value))
    parts = strsplit(char(value), ',');
    out = cellfun(@strtrim, parts, 'UniformOutput', false);
elseif iscell(value)
    out = cell(size(value));
    for i = 1:numel(value)
        out{i} = strtrim(char(value{i}));
    end
else
    error('condition_order must be a CSV string or cell array of labels.');
end

out = out(~cellfun(@isempty, out));
out = unique(out, 'stable');
end

function order = append_trial_types_in_order(order, EEG)
if ~isfield(EEG, 'event') || isempty(EEG.event)
    return;
end

for i = 1:length(EEG.event)
    if isfield(EEG.event(i), 'trial_type')
        label = normalize_label_value(EEG.event(i).trial_type);
        if ~isempty(label) && ~any(strcmp(order, label))
            order{end+1} = label; %#ok<AGROW>
        end
    end
end
end

function label = normalize_label_value(value)
if iscell(value)
    if isempty(value)
        label = '';
    else
        label = strtrim(char(value{1}));
    end
elseif isstring(value)
    label = strtrim(char(value));
elseif ischar(value)
    label = strtrim(value);
else
    label = '';
end
end

function order = load_saved_condition_order(output_dir)
order = {};
mat_file = fullfile(output_dir, 'condition_order.mat');
if exist(mat_file, 'file')
    data = load(mat_file);
    if isfield(data, 'condition_order') && iscell(data.condition_order)
        order = data.condition_order;
        return;
    end
end

text_file = fullfile(output_dir, 'condition_order.txt');
if exist(text_file, 'file')
    order = read_condition_order_text_file(text_file);
end
end

function write_condition_order_text(file_path, condition_order)
fid = fopen(file_path, 'w');
if fid == -1
    warning('Could not write condition order file: %s', file_path);
    return;
end
for i = 1:numel(condition_order)
    label = char(condition_order{i});
    label = strrep(label, sprintf('\t'), ' ');
    label = strrep(label, sprintf('\n'), ' ');
    label = strrep(label, sprintf('\r'), ' ');
    fprintf(fid, '%d\t%s\n', i, label);
end
fclose(fid);
end

function order = read_condition_order_text_file(file_path)
order = {};
fid = fopen(file_path, 'r');
if fid == -1
    return;
end

cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>
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
    if isempty(parts)
        continue;
    end

    if numel(parts) >= 2
        label = strtrim(strjoin(parts(2:end), sprintf('\t')));
    else
        label = strtrim(parts{1});
        token = regexp(label, '^\d+\s+(.+)$', 'tokens', 'once');
        if ~isempty(token)
            label = strtrim(token{1});
        end
    end

    if ~isempty(label)
        order{end + 1} = label; %#ok<AGROW>
    end
end

order = unique(order, 'stable');
end
