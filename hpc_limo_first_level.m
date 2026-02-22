function hpc_limo_first_level(input_folder)
% HPC_LIMO_FIRST_LEVEL - Run LIMO 1st level analysis on epoched EEG data
%
% Usage:
%   hpc_limo_first_level(input_folder)
%
% Inputs:
%   input_folder - Path to the epoch folder containing *_epoch.set files
%
% Outputs:
%   - <input_folder>/limo_first_level/  LIMO 1st level results directory
%     Contains: PairAsso_Epoched.study, LIMO derivatives per subject
%
% Analysis:
%   - Method: Weighted Least Squares (WLS)
%   - Measure: ERP data (daterp)
%   - Conditions: 6 trial types (Study_hits, Study_misses, Test_hits,
%                 Test_misses, Correct_rejections, False_alarms)

% Start timing
tic;
start_time = datetime('now');
fprintf('\n========================================\n');
fprintf('LIMO 1ST LEVEL ANALYSIS\n');
fprintf('Started at: %s\n', start_time);
fprintf('========================================\n\n');

% Force headless MATLAB
set(0,'DefaultFigureVisible','off');
java.lang.System.setProperty('java.awt.headless','true');

% Add EEGLAB to path
addpath('/home/devon7y/scratch/devon7y/eeglab2022.1');
fprintf('Added EEGLAB to path\n');

% Define directories
output_dir = fullfile(input_folder, 'limo_first_level');
if ~exist(output_dir,'dir')
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
    return
end

% Disable EEGLAB online checks
setpref('Internet','EeglabServerCheck','off');
setpref('Internet','EeglabUpdateServer','off');
setpref('Internet','EeglabPluginlistURL','off');

disp('LIMO will manage parallel pool automatically...');

% Define STUDY file path
study_file = fullfile(output_dir, 'PairAsso_Epoched.study');

% Load or create STUDY
if exist(study_file, 'file')
    fprintf('Loading existing STUDY file: %s\n', study_file);
    [STUDY, ALLEEG] = pop_loadstudy('filename', 'PairAsso_Epoched.study', 'filepath', output_dir);
else
    fprintf('Creating new STUDY from epoched files in: %s\n', input_folder);

    % Collect all epoched .set files
    files_in_dir = dir(fullfile(input_folder, '*_epoch.set'));

    % Filter out LIMO-generated files
    epoched_files = [];
    for i = 1:length(files_in_dir)
        fname = files_in_dir(i).name;
        if ~startsWith(fname, 'sub-') && ~contains(fname, '_ses-') && ~contains(fname, '_design')
            epoched_files = [epoched_files; files_in_dir(i)];
        end
    end

    fprintf('Found %d epoched files\n', length(epoched_files));

    if isempty(epoched_files)
        error('No epoched files found in: %s', input_folder);
    end

    % Load all datasets
    ALLEEG = [];
    study_commands = cell(length(epoched_files), 1);

    for i = 1:length(epoched_files)
        fprintf('Loading file %d/%d: %s\n', i, length(epoched_files), epoched_files(i).name);
        EEG = pop_loadset('filename', epoched_files(i).name, 'filepath', epoched_files(i).folder);

        % Parse subject ID from filename: <ID>_epoch.set
        [~, basename, ~] = fileparts(epoched_files(i).name);

        tokens = regexp(basename, '^(\d+)_epoch$', 'tokens');
        if ~isempty(tokens)
            unique_subject_id = tokens{1}{1};
        else
            error('Failed to parse filename: %s (expected format: <ID>_epoch.set)', epoched_files(i).name);
        end

        EEG.subject = unique_subject_id;
        EEG.condition = '';

        study_commands{i} = {'index', i, 'subject', unique_subject_id};
        fprintf('  Dataset %d: Subject %s\n', i, unique_subject_id);

        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, i);
    end

    % Create STUDY
    fprintf('Creating STUDY structure...\n');
    [STUDY, ALLEEG] = std_editset([], ALLEEG, 'name', 'PairAsso_Epoched', ...
        'task', 'Associative Memory', ...
        'filename', 'PairAsso_Epoched.study', ...
        'filepath', output_dir, ...
        'commands', study_commands, ...
        'updatedat', 'on');

    % Define design for 6 experimental conditions
    fprintf('Setting up design for 6 experimental conditions...\n');
    STUDY = std_makedesign(STUDY, ALLEEG, 1, 'name', 'All_Conditions', ...
        'delfiles','off', ...
        'defaultdesign','off', ...
        'variable1','trial_type', ...
        'values1', {'Study_hits', 'Study_misses', 'Test_hits', 'Test_misses', 'Correct_rejections', 'False_alarms'}, ...
        'vartype1','categorical', ...
        'subjselect', {ALLEEG.subject});

    % Save STUDY
    [STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, 'filename', 'PairAsso_Epoched.study', 'filepath', output_dir);
    fprintf('STUDY saved to: %s\n', study_file);
end

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
    STUDY = pop_limo(STUDY, ALLEEG, 'method','WLS','measure','daterp', ...
        'erase','on','splitreg','off','interaction','off');
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
fprintf('Elapsed:     %.2f seconds (%.2f minutes)\n', elapsed, elapsed/60);
fprintf('========================================\n');
disp('Done.');

end
