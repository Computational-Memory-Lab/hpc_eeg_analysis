function hpc_limo_second_level(input_folder, contrast, output_tag)
% HPC_LIMO_SECOND_LEVEL - Run LIMO group-level paired t-test
%
% Usage:
%   hpc_limo_second_level(input_folder, [p1 p2])
%   hpc_limo_second_level(input_folder, {'CondA', 'CondB'})
%   hpc_limo_second_level(input_folder, contrast, output_tag)
%
% Inputs:
%   input_folder - Path to limo_first_level folder
%   contrast     - Either numeric [p1 p2] parameter indices, or
%                  condition-name pair {'CondA','CondB'}
%   output_tag   - Optional folder suffix (e.g., 'study_hits_vs_study_misses')
%
% Outputs:
%   - <parent_of_input_folder>/limo_second_level_<output_tag>/

if nargin < 2 || isempty(contrast)
    error('Provide a contrast: numeric [p1 p2] or labels {''A'',''B''}.');
end
if nargin < 3
    output_tag = '';
end
input_folder = normalize_input_folder(input_folder);

is_numeric_contrast = isnumeric(contrast);
if is_numeric_contrast
    if numel(contrast) ~= 2
        error('Numeric contrast must be [p1 p2].');
    end
    p1 = contrast(1);
    p2 = contrast(2);
    if any(~isfinite([p1 p2])) || any(mod([p1 p2], 1) ~= 0)
        error('Numeric contrast values must be finite integers: [p1 p2].');
    end
    contrast_labels = {};
else
    contrast_labels = parse_contrast_labels(contrast);
    p1 = NaN;
    p2 = NaN;
end

if isempty(output_tag)
    if is_numeric_contrast
        output_tag = sprintf('%d_%d', p1, p2);
    else
        output_tag = sprintf('%s_vs_%s', sanitize_tag(contrast_labels{1}), sanitize_tag(contrast_labels{2}));
    end
end

% Start timing
tic;
start_time = datetime('now');
fprintf('\n========================================\n');
fprintf('LIMO PAIRED T-TEST ANALYSIS\n');
fprintf('Started at: %s\n', start_time);
fprintf('Output tag: %s\n', output_tag);
fprintf('========================================\n\n');

% Force headless MATLAB
set(0, 'DefaultFigureVisible', 'off');
java.lang.System.setProperty('java.awt.headless', 'true');

% Add EEGLAB to path
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1');
addpath('/home/devon7y/scratch/devon7y/hpc_eeg_analysis');
fprintf('Added EEGLAB to path\n');

% Force EEGLAB into offline-safe mode before startup. This avoids
% intermittent plugin metadata fetch failures on Fir.
configure_eeglab_offline_mode();

% Define directories
limo_analysis_dir = input_folder;
pipeline_root = fileparts(input_folder);
if isempty(pipeline_root)
    pipeline_root = pwd;
end
output_dir = fullfile(pipeline_root, sprintf('limo_second_level_%s', output_tag));
if ~exist(output_dir, 'dir')
    fprintf('Creating output directory: %s\n', output_dir);
    [success, msg] = mkdir(output_dir);
    if ~success
        error('Failed to create output directory: %s. Error: %s', output_dir, msg);
    end
end
fprintf('Input (1st level) dir: %s\n', limo_analysis_dir);
fprintf('Output dir:            %s\n', output_dir);

% Start EEGLAB safely
try
    eeglab nogui;
catch ME
    fprintf('EEGLAB startup failed on first attempt; retrying in offline mode...\n');
    disp(getReport(ME));
    configure_eeglab_offline_mode();
    try
        eeglab nogui;
    catch ME2
        disp('EEGLAB failed to start:');
        disp(getReport(ME2));
        return;
    end
end

% Disable EEGLAB online checks
setpref('Internet', 'EeglabServerCheck', 'off');
setpref('Internet', 'EeglabUpdateServer', 'off');
setpref('Internet', 'EeglabPluginlistURL', 'off');

% Load STUDY
study_files = dir(fullfile(limo_analysis_dir, '*.study'));
if isempty(study_files)
    error('No .study file found in: %s', limo_analysis_dir);
end
study_filename = study_files(1).name;
fprintf('Loading STUDY file: %s\n', study_filename);

try
    [STUDY, ALLEEG] = pop_loadstudy('filename', study_filename, 'filepath', limo_analysis_dir); %#ok<ASGLU>
    fprintf('STUDY loaded successfully\n');
catch ME
    fprintf('ERROR: Failed to load STUDY:\n');
    disp(getReport(ME));
    return;
end

% Resolve condition labels to numeric indices if needed
condition_order = load_condition_order(limo_analysis_dir, STUDY);
if ~is_numeric_contrast
    [p1, p2] = labels_to_parameter_indices(contrast_labels, condition_order);
end

if any([p1 p2] < 1)
    error('Invalid parameter indices: [%d %d]', p1, p2);
end

if isempty(contrast_labels) && numel(condition_order) >= max([p1 p2])
    contrast_labels = {condition_order{p1}, condition_order{p2}};
end

fprintf('Resolved parameter contrast: [%d %d]\n', p1, p2);
if ~isempty(contrast_labels)
    fprintf('Resolved label contrast: %s vs %s\n', contrast_labels{1}, contrast_labels{2});
end

% Find LIMO beta files for all subjects
fprintf('Finding LIMO beta files...\n');
limo_files = {};
subject_count = 0;

derivatives_dir = fullfile(limo_analysis_dir, 'derivatives');
limo_dirs = dir(fullfile(derivatives_dir, 'sub-*_*'));
limo_dirs = limo_dirs([limo_dirs.isdir]);

if isempty(limo_dirs)
    limo_dirs = dir(fullfile(derivatives_dir, '*'));
    limo_dirs = limo_dirs([limo_dirs.isdir]);
    limo_dirs = limo_dirs(~ismember({limo_dirs.name}, {'.', '..'}));
end

excluded_subjects = {};
for i = 1:length(limo_dirs)
    beta_search = dir(fullfile(derivatives_dir, limo_dirs(i).name, '**', 'Betas.mat'));
    if isempty(beta_search)
        beta_search = dir(fullfile(derivatives_dir, limo_dirs(i).name, '*', 'Betas.mat'));
    end

    for b = 1:length(beta_search)
        beta_file = fullfile(beta_search(b).folder, beta_search(b).name);
        try
            temp_data = load(beta_file);
            sz = size(temp_data.Betas);
            if sz(3) >= max(p1, p2)
                limo_files{end+1} = beta_file; %#ok<AGROW>
                subject_count = subject_count + 1;
                fprintf('  Found valid Betas file %d: %s\n', subject_count, limo_dirs(i).name);
            else
                fprintf('  WARNING: Skipping %s - has %d parameters (need at least %d)\n', ...
                    limo_dirs(i).name, sz(3), max(p1, p2));
                excluded_subjects{end+1} = limo_dirs(i).name; %#ok<AGROW>
            end
        catch ME
            fprintf('  WARNING: Could not load %s: %s\n', limo_dirs(i).name, ME.message);
            excluded_subjects{end+1} = limo_dirs(i).name; %#ok<AGROW>
        end
    end
end

if isempty(limo_files)
    error('No valid LIMO Betas.mat files found.');
end

fprintf('Found %d valid subjects with LIMO results\n', subject_count);
if ~isempty(excluded_subjects)
    fprintf('Excluded %d subjects\n', length(excluded_subjects));
end

% Create beta list file
beta_list_file = fullfile(output_dir, 'beta_files_list.txt');
fid = fopen(beta_list_file, 'w');
if fid == -1
    error('Failed to create beta file list at: %s', beta_list_file);
end
for i = 1:length(limo_files)
    fprintf(fid, '%s\n', limo_files{i});
end
fclose(fid);

% Locate channel locations file
chanlocs_file = fullfile(limo_analysis_dir, 'derivatives', 'limo_gp_level_chanlocs.mat');
if ~exist(chanlocs_file, 'file')
    error('Channel locations file not found: %s', chanlocs_file);
end

cd(output_dir);

fprintf('\n========================================\n');
fprintf('Running LIMO Paired T-Test\n');
fprintf('Parameters: [%d, %d]\n', p1, p2);
fprintf('Bootstrap iterations: 1000\n');
fprintf('TFCE correction: on\n');
fprintf('========================================\n');

try
    limo_random_select('paired t-test', ...
        chanlocs_file, ...
        'LIMOfiles', ...
        beta_list_file, ...
        'analysis_type', 'Full scalp analysis', ...
        'parameter', [p1 p2], ...
        'type', 'Channels', 'nboot', 1000, 'tfce', 1);
    fprintf('LIMO paired t-test completed successfully\n');
catch ME
    fprintf('ERROR: LIMO paired t-test failed:\n');
    disp(getReport(ME));
end

% Save resolved contrast metadata
meta_file = fullfile(output_dir, 'contrast_metadata.txt');
fid = fopen(meta_file, 'w');
if fid ~= -1
    fprintf(fid, 'p1\t%d\n', p1);
    fprintf(fid, 'p2\t%d\n', p2);
    if ~isempty(contrast_labels)
        fprintf(fid, 'label1\t%s\n', contrast_labels{1});
        fprintf(fid, 'label2\t%s\n', contrast_labels{2});
    end
    fclose(fid);
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
fprintf('Output saved to: %s\n', output_dir);
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

function configure_eeglab_offline_mode()
% Keep EEGLAB from querying remote plugin/update metadata.
setenv('EEGLAB_NO_UPDATE', '1');
setenv('HTTP_PROXY', '');
setenv('HTTPS_PROXY', '');
setenv('http_proxy', '');
setenv('https_proxy', '');

setpref('Internet', 'EeglabServerCheck', 'off');
setpref('Internet', 'EeglabUpdateServer', 'off');
setpref('Internet', 'EeglabPluginlistURL', 'off');
end

function labels = parse_contrast_labels(value)
if iscell(value)
    if numel(value) ~= 2
        error('Cell label contrast must contain exactly 2 entries.');
    end
    labels = {strtrim(char(value{1})), strtrim(char(value{2}))};
elseif ischar(value) || (isstring(value) && isscalar(value))
    parts = strsplit(char(value), ',');
    parts = cellfun(@strtrim, parts, 'UniformOutput', false);
    parts = parts(~cellfun(@isempty, parts));
    if numel(parts) ~= 2
        error('String label contrast must be "CondA,CondB".');
    end
    labels = parts;
else
    error('Unsupported label contrast type.');
end

if isempty(labels{1}) || isempty(labels{2})
    error('Contrast labels must be non-empty.');
end
end

function tag = sanitize_tag(label)
tag = lower(strtrim(label));
tag = regexprep(tag, '[^a-zA-Z0-9]+', '_');
tag = regexprep(tag, '_+', '_');
tag = regexprep(tag, '^_|_$', '');
if isempty(tag)
    tag = 'contrast';
end
end

function condition_order = load_condition_order(limo_analysis_dir, STUDY)
condition_order = {};
mat_file = fullfile(limo_analysis_dir, 'condition_order.mat');
if exist(mat_file, 'file')
    d = load(mat_file);
    if isfield(d, 'condition_order') && iscell(d.condition_order)
        condition_order = d.condition_order;
    end
end

if isempty(condition_order)
    txt_file = fullfile(limo_analysis_dir, 'condition_order.txt');
    if exist(txt_file, 'file')
        condition_order = read_condition_order_text_file(txt_file);
    end
end

if isempty(condition_order)
    try
        if isfield(STUDY, 'design') && ~isempty(STUDY.design)
            if isfield(STUDY.design(1), 'variable') && ~isempty(STUDY.design(1).variable)
                values = STUDY.design(1).variable(1).value;
                if iscell(values)
                    condition_order = cellfun(@(x) strtrim(char(x)), values, 'UniformOutput', false);
                end
            end
        end
    catch
        % Ignore fallback parsing errors
    end
end

condition_order = condition_order(~cellfun(@isempty, condition_order));
condition_order = unique(condition_order, 'stable');
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

function [p1, p2] = labels_to_parameter_indices(labels, condition_order)
if isempty(condition_order)
    error('Condition order is unavailable; cannot map labels to LIMO parameters.');
end

idx1 = find(strcmpi(condition_order, labels{1}), 1, 'first');
idx2 = find(strcmpi(condition_order, labels{2}), 1, 'first');

if isempty(idx1)
    error('Condition label not found in condition order: %s', labels{1});
end
if isempty(idx2)
    error('Condition label not found in condition order: %s', labels{2});
end

p1 = idx1;
p2 = idx2;
end
