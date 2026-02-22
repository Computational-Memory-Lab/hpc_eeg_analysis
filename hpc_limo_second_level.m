function hpc_limo_second_level(input_folder, parameters)
% HPC_LIMO_SECOND_LEVEL - Run LIMO group-level paired t-test
%
% Usage:
%   hpc_limo_second_level(input_folder, parameters)
%
% Inputs:
%   input_folder - Path to the limo_first_level folder (contains the .study
%                  file and derivatives/ directory with subject Betas.mat files)
%   parameters   - 1x2 array of condition indices to compare, e.g. [3 4]
%                  Condition order: 1=Study_hits, 2=Study_misses,
%                  3=Test_hits, 4=Test_misses, 5=Correct_rejections, 6=False_alarms
%
% Outputs:
%   - <input_folder>/limo_second_level_<p1>_<p2>/  LIMO 2nd level results
%     Contains: paired_samples_ttest_parameter_<p1p2>.mat, LIMO.mat, etc.
%
% Analysis:
%   - Paired t-test between conditions specified by parameters
%   - TFCE cluster correction
%   - 1000 bootstrap iterations

if numel(parameters) ~= 2
    error('parameters must be a 1x2 array, e.g. [3 4]');
end

p1 = parameters(1);
p2 = parameters(2);

% Start timing
tic;
start_time = datetime('now');
fprintf('\n========================================\n');
fprintf('LIMO PAIRED T-TEST ANALYSIS\n');
fprintf('Parameters: [%d %d]\n', p1, p2);
fprintf('Started at: %s\n', start_time);
fprintf('========================================\n\n');

% Force headless MATLAB
set(0,'DefaultFigureVisible','off');
java.lang.System.setProperty('java.awt.headless','true');

% Add EEGLAB to path
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1');
fprintf('Added EEGLAB to path\n');

% Define directories
limo_analysis_dir = input_folder;
output_dir = fullfile(input_folder, sprintf('limo_second_level_%d_%d', p1, p2));

if ~exist(output_dir,'dir')
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
    disp('EEGLAB failed to start:');
    disp(getReport(ME));
    return
end

% Disable EEGLAB online checks
setpref('Internet','EeglabServerCheck','off');
setpref('Internet','EeglabUpdateServer','off');
setpref('Internet','EeglabPluginlistURL','off');

% Find and load the STUDY file
study_files = dir(fullfile(limo_analysis_dir, '*.study'));
if isempty(study_files)
    error('No .study file found in: %s', limo_analysis_dir);
end
study_filename = study_files(1).name;
fprintf('Loading STUDY file: %s\n', study_filename);

try
    [STUDY, ALLEEG] = pop_loadstudy('filename', study_filename, 'filepath', limo_analysis_dir);
    fprintf('STUDY loaded successfully\n');
catch ME
    fprintf('ERROR: Failed to load STUDY:\n');
    disp(getReport(ME));
    return;
end

% Find LIMO beta files for all subjects
fprintf('Finding LIMO beta files...\n');
limo_files = {};
subject_count = 0;

derivatives_dir = fullfile(limo_analysis_dir, 'derivatives');
limo_dirs = dir(fullfile(derivatives_dir, 'sub-*_*'));
limo_dirs = limo_dirs([limo_dirs.isdir]);

% Also handle non-prefixed subject directories
if isempty(limo_dirs)
    limo_dirs = dir(fullfile(derivatives_dir, '*'));
    limo_dirs = limo_dirs([limo_dirs.isdir]);
    limo_dirs = limo_dirs(~ismember({limo_dirs.name}, {'.', '..'}));
end

excluded_subjects = {};
for i = 1:length(limo_dirs)
    % Look for Betas.mat files in LIMO subject directories
    beta_search = dir(fullfile(derivatives_dir, limo_dirs(i).name, '**', 'Betas.mat'));
    if isempty(beta_search)
        % Try non-recursive search
        beta_search = dir(fullfile(derivatives_dir, limo_dirs(i).name, '*', 'Betas.mat'));
    end

    for b = 1:length(beta_search)
        beta_file = fullfile(beta_search(b).folder, beta_search(b).name);
        try
            temp_data = load(beta_file);
            sz = size(temp_data.Betas);
            % Validate that Betas has enough parameters for both conditions
            if sz(3) >= max(p1, p2)
                limo_files{end+1} = beta_file;
                subject_count = subject_count + 1;
                fprintf('  Found valid Betas file %d: %s\n', subject_count, limo_dirs(i).name);
            else
                fprintf('  WARNING: Skipping %s - has %d parameters (need at least %d)\n', ...
                    limo_dirs(i).name, sz(3), max(p1, p2));
                excluded_subjects{end+1} = limo_dirs(i).name;
            end
        catch ME
            fprintf('  WARNING: Could not load %s: %s\n', limo_dirs(i).name, ME.message);
            excluded_subjects{end+1} = limo_dirs(i).name;
        end
    end
end

if isempty(limo_files)
    error('No valid LIMO Betas.mat files found. Make sure 1st level analysis completed successfully.');
end

fprintf('Found %d valid subjects with LIMO results\n', subject_count);
if ~isempty(excluded_subjects)
    fprintf('Excluded %d subjects:\n', length(excluded_subjects));
    for i = 1:length(excluded_subjects)
        fprintf('  - %s\n', excluded_subjects{i});
    end
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
fprintf('Beta file list saved to: %s\n', beta_list_file);

% Locate channel locations file
chanlocs_file = fullfile(limo_analysis_dir, 'derivatives', 'limo_gp_level_chanlocs.mat');
if ~exist(chanlocs_file, 'file')
    error('Channel locations file not found: %s', chanlocs_file);
end
fprintf('Channel locations file: %s\n', chanlocs_file);

cd(output_dir);

% Run LIMO paired t-test
fprintf('\n========================================\n');
fprintf('Running LIMO Paired T-Test\n');
fprintf('Parameters: [%d, %d]\n', p1, p2);
fprintf('Bootstrap iterations: 1000\n');
fprintf('Using TFCE for cluster correction\n');
fprintf('========================================\n');

try
    limo_random_select('paired t-test', ...
        chanlocs_file, ...
        'LIMOfiles', ...
        beta_list_file, ...
        'analysis_type','Full scalp analysis', ...
        'parameter', parameters, ...
        'type','Channels','nboot',1000,'tfce',1);
    fprintf('LIMO paired t-test completed successfully\n');
catch ME
    fprintf('ERROR: LIMO paired t-test failed:\n');
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
fprintf('Elapsed:     %.2f seconds (%.2f minutes)\n', elapsed, elapsed/60);
fprintf('========================================\n');
fprintf('Output saved to: %s\n', output_dir);
disp('Done.');

end
