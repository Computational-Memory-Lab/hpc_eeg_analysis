function hpc_raw_to_interpol(input_folder)
% HPC_RAW_TO_INTERPOL - Apply filtering and channel rejection to EEG data
%
% Usage:
%   hpc_raw_to_interpol(input_folder)
%
% Inputs:
%   input_folder - Path to folder containing raw .set files to process
%
% Outputs:
%   - <input_folder>/interpol/preprocessed_full_<ID>.set  Filtered and cleaned EEG datasets
%   - <input_folder>/interpol/preprocessing_summary_<ID>.txt  Per-subject summary log
%
% Processing Steps:
%   1. Flatline removal: Remove broken electrodes
%   2. High-pass filter: 0.1 Hz cutoff
%   3. Low-pass filter: 50 Hz cutoff
%   4. Line noise removal: 60 Hz and 120 Hz (pop_cleanline)
%   5. Kurtosis-based channel rejection: Z-score threshold = 2 SD
%   6. RANSAC-based channel rejection: Correlation threshold = 0.8
%   7. ICA decomposition: Extended Infomax ICA (runica)
%   8. IC classification: ICLabel deep learning classifier
%   9. IC rejection: Remove artifact components (Eye>90%, Muscle>90%, Heart>90%)
%   10. ASR: Artifact Subspace Reconstruction, cutoff = 20 SD
%   11. Bad window removal: CURRENTLY DISABLED
%
% NOTE: EEGLAB should be loaded by the SLURM script before calling this function
%   (via: eeglab nogui; in the matlab -r command)

%% CHECK DEPENDENCIES
if ~exist('pop_loadset', 'file')
    error('EEGLAB not found in MATLAB path. Please add EEGLAB using: addpath(genpath(''/path/to/eeglab''))');
end

% Add all cleanline external dependencies (not auto-loaded by eeglab nogui)
cleanline_external = '/home/devon7y/scratch/devon7y/eeglab2025.0.0/plugins/Cleanline2.1/external';
if exist(cleanline_external, 'dir')
    addpath(genpath(cleanline_external));
    path(path);
    fprintf('Added all cleanline external dependencies to path (high priority)\n');
else
    warning('Cleanline external directory not found at: %s', cleanline_external);
end

% Verify correct dpss is being used
dpss_path = which('dpss');
if contains(dpss_path, 'chronux')
    fprintf('Using chronux dpss (correct for cleanline)\n');
elseif contains(dpss_path, 'Fieldtrip') || contains(dpss_path, 'fieldtrip')
    warning('FieldTrip dpss detected - may cause conflicts. Trying to fix...');
    rmpath(fileparts(dpss_path));
    addpath(genpath(cleanline_external));
end

if ~exist('pop_cleanline', 'file')
    warning('pop_cleanline not found. You may need to install the cleanline plugin.');
end
if ~exist('clean_flatlines', 'file')
    error('clean_flatlines not found. You need to install the clean_rawdata plugin.');
end
if ~exist('clean_channels', 'file')
    error('clean_channels not found. You need to install the clean_rawdata plugin.');
end
if ~exist('clean_asr', 'file')
    error('clean_asr not found. You need to install the clean_rawdata plugin.');
end
if ~exist('clean_windows', 'file')
    error('clean_windows not found. You need to install the clean_rawdata plugin.');
end
if ~exist('pop_iclabel', 'file')
    error('ICLabel plugin not found. Install ICLabel for IC classification.');
end

%% CREATE OUTPUT DIRECTORY
interpol_folder = fullfile(input_folder, 'interpol');
if ~exist(interpol_folder, 'dir')
    mkdir(interpol_folder);
end

%% FIND INPUT FILES
set_files = dir(fullfile(input_folder, '*.set'));
if isempty(set_files)
    error('No .set files found in: %s', input_folder);
end

files_to_process = cell(1, length(set_files));
for i = 1:length(set_files)
    files_to_process{i} = fullfile(input_folder, set_files(i).name);
end

fprintf('\n==============================================\n');
fprintf('  EEG Preprocessing Pipeline\n');
fprintf('==============================================\n');
fprintf('Input folder:     %s\n', input_folder);
fprintf('Output folder:    %s\n', interpol_folder);
fprintf('Files to process: %d\n', length(files_to_process));
fprintf('==============================================\n\n');

%% PROCESS EACH FILE
for file_idx = 1:length(files_to_process)
    processed_file = files_to_process{file_idx};
    [~, filename, ext] = fileparts(processed_file);

    % Extract subject ID from filename
    tokens = regexp(filename, '\d+', 'match');
    if ~isempty(tokens)
        subject_id = str2double(tokens{1});
    else
        error('Could not extract subject ID from filename "%s"', filename);
    end

    fprintf('\n==============================================\n');
    fprintf('  FILE %d of %d: SUBJECT %d\n', file_idx, length(files_to_process), subject_id);
    fprintf('==============================================\n');
    fprintf('Input file: %s\n\n', [filename ext]);

    %% LOAD DATASET
    fprintf('==============================================\n');
    fprintf('  PREPROCESSING SUBJECT %d\n', subject_id);
    fprintf('==============================================\n\n');

    fprintf('STEP 1: Loading dataset...\n');
    EEG = pop_loadset(processed_file);
    fprintf('  Channels: %d\n', EEG.nbchan);
    fprintf('  Samples: %d\n', EEG.pnts);
    fprintf('  Sampling rate: %d Hz\n', EEG.srate);
    fprintf('  Events: %d\n\n', length(EEG.event));

    %% REMOVE FLATLINED CHANNELS
    fprintf('STEP 2: Removing flatlined channels (broken electrodes)...\n');

    initial_channels = EEG.nbchan;
    flatline_threshold = 5 * EEG.srate; % 5 seconds in samples
    bad_channels = [];

    for ch = 1:EEG.nbchan
        channel_data = EEG.data(ch, :);
        diff_data = diff(channel_data);
        zero_diff = (abs(diff_data) < eps * 20);

        if any(zero_diff)
            d = diff([0; zero_diff(:); 0]);
            starts = find(d == 1);
            ends = find(d == -1) - 1;
            run_lengths = ends - starts + 1;
            max_run = max(run_lengths);
        else
            max_run = 0;
        end

        if max_run > flatline_threshold
            bad_channels = [bad_channels ch];
        end
    end

    removed_flat = length(bad_channels);
    if removed_flat > 0
        fprintf('  Detected %d flatlined channel(s): %s\n', removed_flat, mat2str(bad_channels));
        EEG = pop_select(EEG, 'nochannel', bad_channels);
        fprintf('  Removed %d flatlined channel(s)\n', removed_flat);
    else
        fprintf('  No flatlined channels detected\n');
    end
    fprintf('  Remaining channels: %d\n\n', EEG.nbchan);

    %% HIGH-PASS FILTER
    fprintf('STEP 3: Applying high-pass filter (0.1 Hz)...\n');
    EEG = pop_eegfiltnew(EEG, 'locutoff', 0.1, 'plotfreqz', 0);
    fprintf('  Complete\n\n');

    %% LOW-PASS FILTER
    fprintf('STEP 4: Applying low-pass filter (50 Hz)...\n');
    EEG = pop_eegfiltnew(EEG, 'hicutoff', 50, 'plotfreqz', 0);
    fprintf('  Complete\n\n');

    %% LINE NOISE REMOVAL
    fprintf('STEP 5: Removing line noise (60 Hz and 120 Hz)...\n');
    fprintf('  NOTE: Consider alternative methods for line noise removal:\n');
    fprintf('  https://eeglab.org/tutorials/05_Preprocess/Filtering.html#alternative-to-filtering-for-line-noise-removal\n');

    EEG = pop_cleanline(EEG, 'bandwidth', 2, 'chanlist', 1:EEG.nbchan, ...
        'computepower', 1, 'linefreqs', [60 120], 'newversion', 0, ...
        'normSpectrum', 0, 'p', 0.01, 'pad', 2, 'plotfigures', 0, ...
        'scanforlines', 0, 'sigtype', 'Channels', 'taperbandwidth', 2, ...
        'tau', 100, 'verb', 1, 'winsize', 4, 'winstep', 3);
    fprintf('  Complete\n\n');

    %% KURTOSIS-BASED CHANNEL REJECTION
    fprintf('STEP 6: Kurtosis-based channel rejection...\n');
    fprintf('  Method: Z-score normalized kurtosis (threshold = 2 SD)\n');
    channels_before_kurt = EEG.nbchan;

    [EEG, removed_kurt_idx, ~] = pop_rejchan(EEG, 'elec', 1:EEG.nbchan, ...
        'threshold', 2, 'norm', 'on', 'measure', 'kurt');

    removed_kurt = channels_before_kurt - EEG.nbchan;
    if removed_kurt > 0
        fprintf('  Removed %d channel(s) based on kurtosis\n', removed_kurt);
        fprintf('  Removed channel indices: %s\n', mat2str(removed_kurt_idx));
    else
        fprintf('  No channels removed\n');
    end
    fprintf('  Remaining channels: %d\n', EEG.nbchan);
    fprintf('  DEBUG: EEG.nbchan=%d, length(EEG.chanlocs)=%d\n\n', EEG.nbchan, length(EEG.chanlocs));

    %% LOAD CHANNEL LOCATIONS (after channel removal to avoid mismatch)
    fprintf('STEP 7: Loading channel locations...\n');
    chanlocs_file = '/home/devon7y/scratch/devon7y/devon_preprocessing/New_Received2025_AdultAverageNet256_v1.sfp';

    if ~exist(chanlocs_file, 'file')
        error('Channel locations file not found: %s', chanlocs_file);
    end

    all_chanlocs = readlocs(chanlocs_file, 'filetype', 'sfp');
    fprintf('  Read %d channel locations from file\n', length(all_chanlocs));

    if length(all_chanlocs) >= EEG.nbchan
        EEG.chanlocs = all_chanlocs(1:EEG.nbchan);
        fprintf('  Assigned first %d channel locations to dataset\n', EEG.nbchan);
    else
        error('Channel location file has %d channels but dataset has %d', length(all_chanlocs), EEG.nbchan);
    end

    num_chans_with_locs = 0;
    for ch = 1:EEG.nbchan
        if ~isempty(EEG.chanlocs(ch).X) && ~isempty(EEG.chanlocs(ch).Y) && ~isempty(EEG.chanlocs(ch).Z)
            num_chans_with_locs = num_chans_with_locs + 1;
        end
    end
    fprintf('  Successfully loaded locations for %d/%d channels\n\n', num_chans_with_locs, EEG.nbchan);

    %% RANSAC-BASED CHANNEL REJECTION
    fprintf('STEP 8: RANSAC-based channel rejection (correlation)...\n');
    fprintf('  Method: Channel correlation with RANSAC reconstruction\n');
    fprintf('  Correlation threshold: 0.8\n');
    fprintf('  Line noise threshold: 4 SD\n');

    num_chans_with_locs = 0;
    for ch = 1:length(EEG.chanlocs)
        if isfield(EEG.chanlocs, 'X') && ~isempty(EEG.chanlocs(ch).X) && ...
           isfield(EEG.chanlocs, 'Y') && ~isempty(EEG.chanlocs(ch).Y) && ...
           isfield(EEG.chanlocs, 'Z') && ~isempty(EEG.chanlocs(ch).Z)
            num_chans_with_locs = num_chans_with_locs + 1;
        end
    end
    fprintf('  Channels with locations: %d/%d\n', num_chans_with_locs, EEG.nbchan);

    if num_chans_with_locs < EEG.nbchan * 0.5
        error('Insufficient channel locations for RANSAC. Most channels need X,Y,Z coordinates.');
    end

    channels_before_ransac = EEG.nbchan;
    [EEG, removed_ransac_idx] = clean_channels(EEG, 0.8, 4, [], 0.5, 50);
    removed_ransac = channels_before_ransac - EEG.nbchan;

    if removed_ransac > 0
        fprintf('  Removed %d channel(s) based on RANSAC\n', removed_ransac);
        if ~isempty(removed_ransac_idx)
            fprintf('  Removed channel indices: %s\n', mat2str(removed_ransac_idx));
        end
    else
        fprintf('  No channels removed\n');
    end
    fprintf('  Remaining channels: %d\n\n', EEG.nbchan);

    % Save checkpoint before ICA
    before_amica_filename = sprintf('before_step9_amica_%d.set', subject_id);
    EEG = pop_saveset(EEG, 'filename', before_amica_filename, 'filepath', interpol_folder);
    fprintf('  Saved checkpoint: %s\n\n', fullfile(interpol_folder, before_amica_filename));

    %% ICA DECOMPOSITION - RUNICA
    fprintf('STEP 9: Running ICA decomposition (Extended Infomax)...\n');
    fprintf('  Method: Extended Infomax ICA (runica)\n');
    fprintf('  Channels: %d\n', EEG.nbchan);
    fprintf('  NOTE: ICA may take 15-30 minutes for high-density EEG\n');
    fprintf('  Start time: %s\n\n', datestr(now));

    tic;
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, ...
        'interrupt', 'off', 'pca', EEG.nbchan);
    amica_time = toc;

    fprintf('  ICA decomposition complete!\n');
    fprintf('  Time elapsed: %.2f minutes\n', amica_time/60);
    fprintf('  Components computed: %d\n\n', size(EEG.icaweights, 1));

    if isempty(EEG.icachansind) || length(EEG.icachansind) ~= EEG.nbchan
        EEG.icachansind = 1:EEG.nbchan;
        fprintf('  Set icachansind to current channels (1:%d)\n\n', EEG.nbchan);
    end

    %% IC CLASSIFICATION
    fprintf('STEP 10: Classifying independent components (ICLabel)...\n');
    fprintf('  Method: ICLabel deep learning classifier\n');

    EEG = pop_iclabel(EEG, 'default');

    ic_classes = {'Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise', 'Channel Noise', 'Other'};
    fprintf('\n  IC Classification Summary:\n');
    for i = 1:length(ic_classes)
        class_prob = EEG.etc.ic_classification.ICLabel.classifications(:, i);
        high_conf = sum(class_prob > 0.5);
        fprintf('    %s: %d ICs (>50%% confidence)\n', ic_classes{i}, high_conf);
    end
    fprintf('\n');

    %% IC REJECTION
    fprintf('STEP 11: Removing artifact ICs...\n');
    fprintf('  Criteria: Eye>90%%, Muscle>90%%, Heart>90%%\n');

    % [Brain, Muscle, Eye, Heart, Line Noise, Channel Noise, Other]
    EEG = pop_icflag(EEG, [NaN NaN;    % Brain - keep all
                            0.9 1;      % Muscle - reject if >90%
                            0.9 1;      % Eye - reject if >90%
                            0.9 1;      % Heart - reject if >90%
                            NaN NaN;    % Line Noise - keep
                            NaN NaN;    % Channel Noise - keep
                            NaN NaN]);  % Other - keep

    num_rejected_ics = sum(EEG.reject.gcompreject);
    rejected_ic_indices = find(EEG.reject.gcompreject);

    fprintf('  Flagged %d artifact components for removal\n', num_rejected_ics);
    if num_rejected_ics > 0
        fprintf('  Rejected IC indices: %s\n', mat2str(rejected_ic_indices));
    end

    EEG = pop_subcomp(EEG, rejected_ic_indices, 0);
    fprintf('  Remaining components: %d\n', size(EEG.icaweights, 1));
    fprintf('  Back-projected to channel space\n\n');

    %% ASR - ARTIFACT SUBSPACE RECONSTRUCTION
    fprintf('STEP 12: ASR - Artifact Subspace Reconstruction...\n');
    fprintf('  Method: Temporal burst artifact correction\n');
    fprintf('  Cutoff threshold: 20 SD (conservative)\n');

    original_data_rms = sqrt(mean(EEG.data(:).^2));
    EEG = clean_asr(EEG, 20, [], [], [], 0.075, [-inf 5.5]);
    corrected_data_rms = sqrt(mean(EEG.data(:).^2));
    percent_changed = 100 * abs(original_data_rms - corrected_data_rms) / original_data_rms;

    fprintf('  ASR correction applied\n');
    fprintf('  Data RMS change: %.2f%%\n', percent_changed);
    fprintf('  Channels preserved: %d\n\n', EEG.nbchan);

    %% REMOVE BAD TIME WINDOWS - SKIPPED FOR NOW
    % NOTE: This step is currently disabled because removing time windows creates
    % discontinuities in the data. Enable for epoched/trial-based analyses.
    fprintf('STEP 13: Skipped (bad window removal disabled - see comments)\n\n');

    %% SAVE DATASET
    fprintf('STEP 14: Saving preprocessed dataset...\n');
    output_filename = sprintf('preprocessed_full_%d.set', subject_id);
    EEG.setname = sprintf('preprocessed_full_%d', subject_id);
    EEG = pop_saveset(EEG, 'filename', output_filename, 'filepath', interpol_folder);
    fprintf('  Saved: %s\n\n', output_filename);

    %% BUILD SUMMARY
    summary_lines = {
        '==============================================';
        '  PREPROCESSING COMPLETE!';
        '==============================================';
        sprintf('Output file: %s', output_filename);
        sprintf('Location:    %s', interpol_folder);
        '';
        'Processing summary:';
        sprintf('  Initial channels:        %d', initial_channels);
        sprintf('  Flatlined channels:      %d removed', removed_flat);
        sprintf('  Kurtosis rejection:      %d removed', removed_kurt);
        sprintf('  RANSAC rejection:        %d removed', removed_ransac);
        sprintf('  Final channels:          %d', EEG.nbchan);
        '';
        '  Filters applied:';
        '    High-pass:             0.1 Hz';
        '    Low-pass:              50 Hz';
        '    Line noise removed:    60 Hz, 120 Hz';
        '';
        '  ICA decomposition:';
        '    Method:                Extended Infomax (runica)';
        sprintf('    Components computed:   %d', size(EEG.icaweights, 1) + num_rejected_ics);
        sprintf('    Artifact ICs removed:  %d', num_rejected_ics);
        sprintf('    Retained ICs:          %d', size(EEG.icaweights, 1));
        '';
        sprintf('  ASR correction:          %.2f%% RMS change', percent_changed);
        '  Time windows removed:    NONE (step disabled)';
        '';
        sprintf('  Final duration:          %.2f seconds', EEG.pnts / EEG.srate);
        sprintf('  Events preserved:        %d', length(EEG.event));
        '==============================================';
        ''
    };

    % Print summary to console
    for line_idx = 1:length(summary_lines)
        fprintf('%s\n', summary_lines{line_idx});
    end

    % Write summary to .txt file
    summary_txt_file = fullfile(interpol_folder, sprintf('preprocessing_summary_%d.txt', subject_id));
    fid = fopen(summary_txt_file, 'w');
    if fid ~= -1
        for line_idx = 1:length(summary_lines)
            fprintf(fid, '%s\n', summary_lines{line_idx});
        end
        fclose(fid);
        fprintf('  Summary written to: %s\n', summary_txt_file);
    else
        warning('Could not write summary file: %s', summary_txt_file);
    end

end

fprintf('\n==============================================\n');
fprintf('  ALL FILES PROCESSED!\n');
fprintf('==============================================\n');
fprintf('Total files processed: %d\n', length(files_to_process));
fprintf('Output folder: %s\n', interpol_folder);
fprintf('==============================================\n\n');

end
