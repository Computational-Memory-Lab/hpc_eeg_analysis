function hpc_set_to_interpol(input_folder, subject_filter)
% HPC_SET_TO_INTERPOL - Apply filtering, channel rejection, and interpolation to EEG data
%
% Usage:
%   hpc_set_to_interpol(input_folder)
%   hpc_set_to_interpol(input_folder, subject_filter)
%
% Inputs:
%   input_folder - Path to folder containing behaviorally aligned .set files to process
%   subject_filter - Optional numeric subject ID. When provided, process
%                    only that subject.
%
% Outputs:
%   - <parent_of_input_folder>/interpol/preprocessed_full_<ID>.set  Filtered and cleaned EEG datasets
%   - <parent_of_input_folder>/interpol/preprocessing_summary_<ID>.txt  Per-subject summary log
%
% Processing Steps:
%   1. Flatline removal: Remove broken electrodes
%   2. High-pass filter: 0.1 Hz cutoff
%   3. Low-pass filter: 50 Hz cutoff
%   4. Line noise removal: 60 Hz and 120 Hz (pop_cleanline)
%   5. Kurtosis-based channel rejection: Z-score threshold = 2 SD
%   6. Load template channel locations for surviving channels
%   7. RANSAC-based channel rejection: Correlation threshold = 0.8
%   8. ICA decomposition: Extended Infomax ICA (runica)
%   9. IC classification: ICLabel deep learning classifier
%   10. IC rejection: Remove artifact components (Eye>90%, Muscle>90%, Heart>90%)
%   11. ASR: Artifact Subspace Reconstruction, cutoff = 20 SD
%   12. Channel interpolation: Reconstruct removed channels (spherical)
%   13. Re-reference: Average reference
%
% NOTE: EEGLAB should be loaded by the SLURM script before calling this function
%   (via: eeglab nogui; in the matlab -r command)

if nargin < 1 || nargin > 2
    error('Usage: hpc_set_to_interpol(input_folder [, subject_filter])');
end

if nargin < 2
    subject_filter = [];
else
    subject_filter = parse_optional_subject_filter(subject_filter);
end
input_folder = normalize_input_folder(input_folder);

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
if ~exist('pop_iclabel', 'file')
    error('ICLabel plugin not found. Install ICLabel for IC classification.');
end
if ~exist('pop_interp', 'file')
    error('pop_interp not found. Ensure EEGLAB interpolation functions are available.');
end
if ~exist('pop_reref', 'file')
    error('pop_reref not found. Ensure EEGLAB rereferencing functions are available.');
end

chanlocs_file = '/home/devon7y/scratch/devon7y/devon_preprocessing/New_Received2025_AdultAverageNet256_v1.sfp';
if ~exist(chanlocs_file, 'file')
    error('Channel locations file not found: %s', chanlocs_file);
end
raw_chanlocs = readlocs(chanlocs_file, 'filetype', 'sfp');
if isempty(raw_chanlocs)
    error('No channel locations loaded from: %s', chanlocs_file);
end
[all_chanlocs, template_index_info] = select_data_chanlocs(raw_chanlocs);
fprintf(['Loaded %d channel-location entries from template and selected %d EEG data channels ' ...
         'for pipeline use.\n'], numel(raw_chanlocs), numel(all_chanlocs));
if numel(raw_chanlocs) ~= numel(all_chanlocs)
    fprintf('  Excluded non-data template entries at rows: %s\n\n', mat2str(template_index_info.excluded_rows));
else
    fprintf('\n');
end

%% CREATE OUTPUT DIRECTORY
pipeline_root = fileparts(input_folder);
if isempty(pipeline_root)
    pipeline_root = pwd;
end
interpol_folder = fullfile(pipeline_root, 'interpol');
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
processed_count = 0;
skipped_existing_count = 0;
filtered_out_count = 0;
matched_subject_count = 0;
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

    if ~isempty(subject_filter) && subject_id ~= subject_filter
        filtered_out_count = filtered_out_count + 1;
        continue;
    end
    matched_subject_count = matched_subject_count + 1;

    output_filename = sprintf('preprocessed_full_%d.set', subject_id);
    output_filepath = fullfile(interpol_folder, output_filename);
    if exist(output_filepath, 'file')
        fprintf('\nSkipping subject %d (already exists): %s\n', subject_id, output_filepath);
        skipped_existing_count = skipped_existing_count + 1;
        continue;
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

    if length(all_chanlocs) ~= EEG.nbchan
        error(['Data-channel template has %d channels but dataset has %d. ' ...
               'Check fiducial/non-data entries in %s.'], length(all_chanlocs), EEG.nbchan, chanlocs_file);
    end

    %% REMOVE FLATLINED CHANNELS
    fprintf('STEP 2: Removing flatlined channels (broken electrodes)...\n');

    initial_channels = EEG.nbchan;
    original_channel_indices = 1:initial_channels;
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
        original_channel_indices = remove_index_positions(original_channel_indices, bad_channels);
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
        original_channel_indices = remove_index_positions(original_channel_indices, removed_kurt_idx);
    else
        fprintf('  No channels removed\n');
    end
    fprintf('  Remaining channels: %d\n', EEG.nbchan);
    fprintf('\n');

    %% LOAD CHANNEL LOCATIONS (aligned to surviving original channel indices)
    fprintf('STEP 7: Loading channel locations...\n');
    if numel(original_channel_indices) ~= EEG.nbchan
        error('Channel index tracking mismatch before chanloc assignment: tracked=%d, EEG.nbchan=%d', ...
            numel(original_channel_indices), EEG.nbchan);
    end
    EEG.chanlocs = all_chanlocs(original_channel_indices);
    fprintf('  Loaded %d template channel locations\n', length(all_chanlocs));
    fprintf('  Assigned location subset for %d surviving channels\n', EEG.nbchan);

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
        original_channel_indices = remove_index_positions(original_channel_indices, removed_ransac_idx);
    else
        fprintf('  No channels removed\n');
    end
    fprintf('  Remaining channels: %d\n\n', EEG.nbchan);

    % Rebuild all channel-location metadata after RANSAC. clean_channels
    % removes channels from EEG.data but can leave chanloc-related metadata
    % inconsistent, which later causes eeg_checkset/pop_saveset to clear
    % locations and makes ICLabel fail inside pop_reref.
    EEG = restore_chanloc_metadata(EEG, all_chanlocs, original_channel_indices, initial_channels);
    EEG = eeg_checkset(EEG);
    assert_chanloc_metadata(EEG, 'post-RANSAC checkpoint');
    fprintf('  Rebuilt channel metadata for %d retained channels and %d removed channels\n\n', ...
        EEG.nbchan, initial_channels - EEG.nbchan);

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
    % Silence runica's step-by-step stdout logging to avoid SLURM/MATLAB
    % batch-stream failures during long ICA runs on Fir.
    EEG = pop_runica(EEG, 'icatype', 'runica', 'options', ...
        {'extended', 1, 'interrupt', 'off', 'pca', EEG.nbchan, 'verbose', 'off'});
    amica_time = toc;

    fprintf('  ICA decomposition complete!\n');
    fprintf('  Time elapsed: %.2f minutes\n', amica_time/60);
    fprintf('  Components computed: %d\n\n', size(EEG.icaweights, 1));

    if isempty(EEG.icachansind)
        EEG.icachansind = 1:EEG.nbchan;
        fprintf('  Set icachansind to current channels (1:%d)\n\n', EEG.nbchan);
    end
    EEG.icachansind = reshape(double(EEG.icachansind), 1, []);
    EEG = eeg_checkset(EEG);
    assert_iclabel_ready(EEG);

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

    % Record ICA stats before clearing the decomposition fields.
    % ICA has been fully applied (back-projected); retaining icaweights/icachansind
    % after channel interpolation causes EEGLAB to silently remove the decomposition
    % in pop_reref because icachansind no longer covers all channels.
    n_ica_computed = size(EEG.icaweights, 1) + num_rejected_ics;
    n_ica_retained = size(EEG.icaweights, 1);
    EEG.icaweights  = [];
    EEG.icasphere   = [];
    EEG.icachansind = [];
    EEG.icawinv     = [];
    if isfield(EEG, 'icaact'), EEG.icaact = []; end
    fprintf('  ICA fields cleared (decomposition already applied; prevents rereferencing mismatch)\n\n');

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

    %% INTERPOLATE REMOVED CHANNELS BACK TO FULL MONTAGE
    fprintf('STEP 13: Interpolating removed channels (spherical)...\n');
    channels_before_interp = EEG.nbchan;
    channels_interpolated = 0;
    if channels_before_interp < initial_channels
        EEG = pop_interp(EEG, all_chanlocs(1:initial_channels), 'spherical');
        channels_interpolated = EEG.nbchan - channels_before_interp;
        fprintf('  Interpolated %d channel(s)\n', channels_interpolated);
        fprintf('  Restored channels: %d -> %d\n', channels_before_interp, EEG.nbchan);
    else
        fprintf('  No interpolation needed (no channels removed)\n');
    end
    fprintf('\n');

    %% RE-REFERENCE TO AVERAGE REFERENCE
    fprintf('STEP 14: Re-referencing to average reference...\n');
    EEG = pop_reref(EEG, []);
    fprintf('  Complete\n\n');

    %% SAVE DATASET
    fprintf('STEP 15: Saving preprocessed dataset...\n');
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
        sprintf('  Interpolated channels:   %d restored', channels_interpolated);
        sprintf('  Final channels:          %d', EEG.nbchan);
        '';
        '  Filters applied:';
        '    High-pass:             0.1 Hz';
        '    Low-pass:              50 Hz';
        '    Line noise removed:    60 Hz, 120 Hz';
        '    Reference:             Average reference';
        '';
        '  ICA decomposition:';
        '    Method:                Extended Infomax (runica)';
        sprintf('    Components computed:   %d', n_ica_computed);
        sprintf('    Artifact ICs removed:  %d', num_rejected_ics);
        sprintf('    Retained ICs:          %d', n_ica_retained);
        '';
        sprintf('  ASR correction:          %.2f%% RMS change', percent_changed);
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

    processed_count = processed_count + 1;
end

if ~isempty(subject_filter) && matched_subject_count == 0
    error('Subject %d not found in input folder: %s', subject_filter, input_folder);
end

fprintf('\n==============================================\n');
fprintf('  ALL FILES PROCESSED!\n');
fprintf('==============================================\n');
fprintf('Processed subjects: %d\n', processed_count);
fprintf('Skipped existing:   %d\n', skipped_existing_count);
if ~isempty(subject_filter)
    fprintf('Filtered out:       %d\n', filtered_out_count);
end
fprintf('Output folder: %s\n', interpol_folder);
fprintf('==============================================\n\n');

end

function out = remove_index_positions(vector_in, positions_to_remove)
out = vector_in;
if isempty(positions_to_remove)
    return;
end

if islogical(positions_to_remove)
    remove_idx = find(positions_to_remove);
else
    remove_idx = unique(round(double(positions_to_remove(:)')));
end

remove_idx = remove_idx(isfinite(remove_idx) & remove_idx >= 1 & remove_idx <= numel(out));
if isempty(remove_idx)
    return;
end
out(remove_idx) = [];
end

function [data_chanlocs, index_info] = select_data_chanlocs(raw_chanlocs)
if ~isstruct(raw_chanlocs) || isempty(raw_chanlocs)
    error('Channel-location template is empty or malformed.');
end
if ~isfield(raw_chanlocs, 'labels')
    error('Channel-location template is missing labels, so data channels cannot be identified.');
end

labels = cell(1, numel(raw_chanlocs));
e_numbers = nan(1, numel(raw_chanlocs));
for idx = 1:numel(raw_chanlocs)
    labels{idx} = strtrim(char(raw_chanlocs(idx).labels));
    token = regexp(labels{idx}, '^E(\d+)$', 'tokens', 'once');
    if ~isempty(token)
        e_numbers(idx) = str2double(token{1});
    end
end

e_rows = find(~isnan(e_numbers));
if ~isempty(e_rows)
    [sorted_numbers, order] = sort(e_numbers(e_rows));
    if numel(unique(sorted_numbers)) ~= numel(sorted_numbers)
        error('Channel-location template contains duplicate E-channel labels.');
    end
    expected_numbers = 1:numel(sorted_numbers);
    if ~isequal(sorted_numbers, expected_numbers)
        error(['Channel-location template E-channel labels must be contiguous from E1 to E%d. ' ...
               'Found: %s'], numel(sorted_numbers), mat2str(sorted_numbers));
    end

    data_rows = e_rows(order);
    data_chanlocs = raw_chanlocs(data_rows);
else
    non_data_mask = false(1, numel(raw_chanlocs));
    for idx = 1:numel(labels)
        upper_label = upper(labels{idx});
        non_data_mask(idx) = startsWith(upper_label, 'FID') || ...
            any(strcmp(upper_label, {'CZ', 'NZ', 'LPA', 'RPA'}));
    end
    data_rows = find(~non_data_mask);
    data_chanlocs = raw_chanlocs(data_rows);
    if isempty(data_chanlocs)
        error('Channel-location template does not contain identifiable EEG data channels.');
    end
end

index_info = struct( ...
    'data_rows', data_rows, ...
    'excluded_rows', setdiff(1:numel(raw_chanlocs), data_rows, 'stable'));
end

function EEG = restore_chanloc_metadata(EEG, all_chanlocs, original_channel_indices, initial_channels)
if numel(original_channel_indices) ~= EEG.nbchan
    error(['Channel metadata rebuild failed: tracked surviving channels (%d) ' ...
           'do not match EEG.nbchan (%d).'], numel(original_channel_indices), EEG.nbchan);
end
if initial_channels > numel(all_chanlocs)
    error(['Channel metadata rebuild failed: template has %d channels but ' ...
           'initial montage expects %d.'], numel(all_chanlocs), initial_channels);
end

EEG.chanlocs = all_chanlocs(original_channel_indices);
EEG.urchanlocs = all_chanlocs(1:initial_channels);

removed_channel_indices = setdiff(1:initial_channels, original_channel_indices, 'stable');
if ~isfield(EEG, 'chaninfo') || ~isstruct(EEG.chaninfo) || isempty(EEG.chaninfo)
    EEG.chaninfo = struct();
end
if isempty(removed_channel_indices)
    EEG.chaninfo.removedchans = [];
else
    EEG.chaninfo.removedchans = all_chanlocs(removed_channel_indices);
end
end

function assert_chanloc_metadata(EEG, context_label)
if nargin < 2 || isempty(context_label)
    context_label = 'channel metadata validation';
end

if isempty(EEG.chanlocs)
    error('%s failed: EEG.chanlocs is empty.', context_label);
end
if numel(EEG.chanlocs) ~= EEG.nbchan
    error(['%s failed: EEG.chanlocs count (%d) does not match EEG.nbchan (%d).'], ...
          context_label, numel(EEG.chanlocs), EEG.nbchan);
end

missing_xyz = false(1, EEG.nbchan);
missing_labels = false(1, EEG.nbchan);
for ch = 1:EEG.nbchan
    missing_labels(ch) = ~isfield(EEG.chanlocs, 'labels') || isempty(EEG.chanlocs(ch).labels);
    missing_xyz(ch) = ~isfield(EEG.chanlocs, 'X') || isempty(EEG.chanlocs(ch).X) || ...
        ~isfield(EEG.chanlocs, 'Y') || isempty(EEG.chanlocs(ch).Y) || ...
        ~isfield(EEG.chanlocs, 'Z') || isempty(EEG.chanlocs(ch).Z);
end

if any(missing_labels)
    error('%s failed: missing channel labels at indices %s.', ...
          context_label, mat2str(find(missing_labels)));
end
if any(missing_xyz)
    error('%s failed: missing X/Y/Z coordinates at channel indices %s.', ...
          context_label, mat2str(find(missing_xyz)));
end
end

function assert_iclabel_ready(EEG)
assert_chanloc_metadata(EEG, 'ICLabel readiness check');

if isempty(EEG.icaweights) || isempty(EEG.icasphere) || isempty(EEG.icawinv)
    error(['ICLabel readiness check failed: ICA decomposition is incomplete ' ...
           '(icaweights/icasphere/icawinv missing).']);
end
if isempty(EEG.icachansind)
    error('ICLabel readiness check failed: EEG.icachansind is empty.');
end

expected_icachansind = 1:EEG.nbchan;
if ~isequal(EEG.icachansind, expected_icachansind)
    error(['ICLabel readiness check failed: EEG.icachansind must equal 1:EEG.nbchan ' ...
           'for this pipeline. Got %s for nbchan=%d.'], mat2str(EEG.icachansind), EEG.nbchan);
end

ica_chanlocs = EEG.chanlocs(EEG.icachansind);
missing_ica_xyz = false(1, numel(ica_chanlocs));
for ch = 1:numel(ica_chanlocs)
    missing_ica_xyz(ch) = isempty(ica_chanlocs(ch).X) || isempty(ica_chanlocs(ch).Y) || isempty(ica_chanlocs(ch).Z);
end
if any(missing_ica_xyz)
    error(['ICLabel readiness check failed: retained ICA channels are missing X/Y/Z ' ...
           'coordinates at indices %s.'], mat2str(find(missing_ica_xyz)));
end
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
