clear; clc;
%% preprocessHaptics.m
% =========================================================================
% Haptic Feedback Preprocessing and CSV/JSON Export Script
% =========================================================================
%
% PURPOSE:
%   Reads raw trial data from the perceptual crossing experiment (PCE)
%   folder structure, extracts the haptic feedback (vibration) time series
%   for each participant in each trial, downsamples to a manageable rate,
%   and exports a consolidated CSV file with a BIDS-inspired JSON metadata
%   sidecar.
%
% EXTRACTION:
%   Each raw trial CSV contains sample-level columns:
%     motor_0_vibrate_software  (participant 1 haptic feedback, 0 or 1)
%     motor_1_vibrate_software  (participant 2 haptic feedback, 0 or 1)
%     timestamp                 (elapsed time in seconds from trial onset)
%
%   The script extracts each participant's haptic signal as a separate
%   row-per-sample time series, identified by DyadID, ParticipantID,
%   and TrialNum.
%
% DOWNSAMPLING:
%   Raw data are typically recorded at ~1000 Hz. Because haptic feedback
%   is a binary (on/off) signal, we downsample by nearest-neighbour
%   decimation: for each output time bin, the raw sample closest to the
%   bin centre is selected, preserving the binary nature of the signal.
%   Default target rate: 100 Hz (preserving binary signal fidelity).
%
% INPUT:
%   Run this script from the root data folder containing pce* subfolders.
%   Each subfolder must contain a trials/ directory with CSV files named:
%     pair_XX_trial_Y.csv
%
% OUTPUT:
%   HapticFeedback.csv   - All trials, long format
%   HapticFeedback.json  - BIDS-style metadata sidecar
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   March 2026
% =========================================================================

%% Parameters
trial_duration   = 60;          % seconds
num_trials_exp   = 18;          % expected trials per dyad
num_participants = 2;           % participants per dyad
fs_target        = 100;         % target sampling rate (Hz) after downsampling

fprintf('==========================================================\n');
fprintf('  Haptic Feedback Preprocessing Script (PCE)\n');
fprintf('==========================================================\n');
fprintf('  Trial duration:           %d s\n', trial_duration);
fprintf('  Expected trials/dyad:     %d\n', num_trials_exp);
fprintf('  Participants per dyad:    %d\n', num_participants);
fprintf('  Target sampling rate:     %d Hz\n', fs_target);
fprintf('==========================================================\n\n');

%% Discover experiment folders
mainDir = pwd;
folders = dir(fullfile(mainDir, 'pce*'));
folders = folders([folders.isdir]);

fprintf('Found %d experiment folders.\n\n', length(folders));

if isempty(folders)
    error(['No pce* folders found in:\n  %s\n' ...
           'Please run this script from the root data directory.'], mainDir);
end

%% Preallocate
sensation_rows = {};
row_count      = 0;

total_trials_found    = 0;
total_trials_skipped  = 0;
dyad_ids_found        = [];
fs_estimates          = [];   % track estimated sampling rates
haptic_means_all      = [];   % track per-segment mean for summary

%% Main loop
for f = 1:length(folders)
    folder_name = folders(f).name;
    trials_path = fullfile(mainDir, folder_name, 'trials');

    % === Extract pair number (XX in pceXXYYMMDD) ===
    dyad_num_str = regexp(folder_name, 'pce(\d+)', 'tokens');
    if isempty(dyad_num_str)
        warning('  Could not parse pair number from folder: %s. Skipping.', folder_name);
        continue;
    end
    dyad_id = str2double(dyad_num_str{1}{1});

    if ~ismember(dyad_id, dyad_ids_found)
        dyad_ids_found = [dyad_ids_found, dyad_id]; %#ok<AGROW>
    end

    fprintf('Processing folder: %s (Dyad %02d)\n', folder_name, dyad_id);

    if ~isfolder(trials_path)
        warning('  No trials/ subfolder found in %s. Skipping.', folder_name);
        continue;
    end

    csv_files = dir(fullfile(trials_path, '*.csv'));

    for j = 1:length(csv_files)
        csv_name = csv_files(j).name;
        csv_path = fullfile(trials_path, csv_name);

        % Parse trial number from filename
        trial_tok = regexp(csv_name, 'trial_?(\d+)', 'tokens');
        if isempty(trial_tok), continue; end
        trial_num = str2double(trial_tok{1}{1});
        total_trials_found = total_trials_found + 1;

        % --- Read the raw trial CSV ---
        try
            data = readtable(csv_path, 'VariableNamingRule', 'preserve');
        catch ME
            warning('  Could not read %s: %s. Skipping.', csv_name, ME.message);
            total_trials_skipped = total_trials_skipped + 1;
            continue;
        end

        % Verify required columns exist
        required_cols = {'timestamp', 'motor_0_vibrate_software', 'motor_1_vibrate_software'};
        if ~all(ismember(required_cols, data.Properties.VariableNames))
            warning('  %s: missing required columns. Skipping.', csv_name);
            total_trials_skipped = total_trials_skipped + 1;
            continue;
        end

        timestamps = data.timestamp;
        motor0     = data.motor_0_vibrate_software;
        motor1     = data.motor_1_vibrate_software;

        n_raw = length(timestamps);
        if n_raw < 2
            warning('  %s: too few samples (%d). Skipping.', csv_name, n_raw);
            total_trials_skipped = total_trials_skipped + 1;
            continue;
        end

        % Estimate raw sampling rate from median inter-sample interval
        dt_raw = median(diff(timestamps));
        fs_raw = 1 / dt_raw;
        fs_estimates = [fs_estimates, fs_raw]; %#ok<AGROW>

        % --- Downsample by nearest-neighbour decimation ---
        % For each output time point (bin centre), select the raw sample
        % closest in time. This preserves the binary nature of the signal.
        bin_width = 1 / fs_target;
        t_max     = timestamps(end);
        bin_edges = 0 : bin_width : t_max;
        n_bins    = length(bin_edges) - 1;

        if n_bins < 1
            warning('  %s: trial too short for downsampling. Skipping.', csv_name);
            total_trials_skipped = total_trials_skipped + 1;
            continue;
        end

        % Bin centres as the output time vector
        time_ds = bin_edges(1:end-1) + bin_width / 2;

        % For each bin centre, find the nearest raw sample
        motor0_ds = zeros(1, n_bins);
        motor1_ds = zeros(1, n_bins);

        for b = 1:n_bins
            [~, nearest_idx] = min(abs(timestamps - time_ds(b)));
            motor0_ds(b) = motor0(nearest_idx);
            motor1_ds(b) = motor1(nearest_idx);
        end

        % --- Store one row per participant ---
        % Participant 1 (motor_0)
        row_count = row_count + 1;
        sensation_rows{row_count} = {dyad_id, 1, trial_num, time_ds, motor0_ds}; %#ok<SAGROW>
        haptic_means_all = [haptic_means_all, mean(motor0_ds)]; %#ok<AGROW>

        % Participant 2 (motor_1)
        row_count = row_count + 1;
        sensation_rows{row_count} = {dyad_id, 2, trial_num, time_ds, motor1_ds}; %#ok<SAGROW>
        haptic_means_all = [haptic_means_all, mean(motor1_ds)]; %#ok<AGROW>
    end
end

fprintf('\n==========================================================\n');
fprintf('  Data collection complete.\n');
fprintf('  Trial files processed:  %d\n', total_trials_found);
fprintf('  Trial files skipped:    %d\n', total_trials_skipped);
fprintf('  Segments extracted:     %d (2 per trial file)\n', row_count);
fprintf('==========================================================\n\n');

%% Sort by DyadID, ParticipantID, TrialNum
sort_keys = zeros(row_count, 3);
for i = 1:row_count
    sort_keys(i, :) = [sensation_rows{i}{1}, sensation_rows{i}{2}, sensation_rows{i}{3}];
end
[~, sort_idx] = sortrows(sort_keys, [1 2 3]);
sensation_rows = sensation_rows(sort_idx);

%% Export CSV
output_csv = 'HapticFeedback.csv';
fprintf('Writing %s ...\n', output_csv);

fid = fopen(output_csv, 'w');
fprintf(fid, 'DyadID,ParticipantID,TrialNum,Time_s,HapticFeedback\n');

total_samples = 0;

for i = 1:row_count
    row       = sensation_rows{i};
    d_id      = row{1};
    p_id      = row{2};
    t_num     = row{3};
    time_vec  = row{4};
    haptic_vec = row{5};

    for s = 1:length(time_vec)
        fprintf(fid, '%d,%d,%d,%.4f,%.6f\n', ...
            d_id, p_id, t_num, time_vec(s), haptic_vec(s));
    end
    total_samples = total_samples + length(time_vec);
end

fclose(fid);

csv_info = dir(output_csv);
fprintf('  -> %s (%.1f MB, %d samples)\n', output_csv, csv_info.bytes / 1e6, total_samples);

%% Summary statistics
if ~isempty(haptic_means_all)
    mean_haptic   = mean(haptic_means_all);
    std_haptic    = std(haptic_means_all);
    median_haptic = median(haptic_means_all);
else
    mean_haptic   = NaN;
    std_haptic    = NaN;
    median_haptic = NaN;
end

if ~isempty(fs_estimates)
    mean_fs_raw = mean(fs_estimates);
    std_fs_raw  = std(fs_estimates);
else
    mean_fs_raw = NaN;
    std_fs_raw  = NaN;
end

%% ========================================================================
%  JSON Metadata Sidecar (BIDS-inspired)
%  ========================================================================

timestamp_str = datestr(now, 'yyyy-mm-ddTHH:MM:SS');

meta = struct();

% Dataset-level
meta.Name = 'Perceptual Crossing Experiment - Haptic Feedback';
meta.Description = ['Haptic feedback (vibration motor activation) time series ' ...
    'extracted from perceptual crossing experiment trials. Each row represents ' ...
    'one downsampled time point for one participant in one trial. The ' ...
    'HapticFeedback column gives the binary state of the vibration motor ' ...
    'at each ' sprintf('%.0f', 1000/fs_target) ' ms time point ' ...
    '(0 = off, 1 = on). Trials last ' num2str(trial_duration) ' seconds.'];

% Sampling
meta.SamplingFrequency     = fs_target;
meta.SamplingFrequencyUnit = 'Hz';
meta.StartTime             = 0;

% Columns
meta.Columns = {'DyadID', 'ParticipantID', 'TrialNum', 'Time_s', 'HapticFeedback'};

meta.DyadID = struct( ...
    'LongName', 'Dyad Identifier', ...
    'Description', ['Numeric identifier for the experimental dyad (pair number 1-32). ' ...
        'Extracted from the pce folder name (e.g., pce01230807 -> pair 01 -> DyadID 1).']);

meta.ParticipantID = struct( ...
    'LongName', 'Participant Identifier', ...
    'Description', ['Participant number within the dyad. ' ...
        '1 = motor_0_vibrate_software, 2 = motor_1_vibrate_software.']);

meta.TrialNum = struct( ...
    'LongName', 'Trial Number', ...
    'Description', ['Trial number within the session (1-' num2str(num_trials_exp) ').']);

meta.Time_s = struct( ...
    'LongName', 'Time', ...
    'Description', 'Elapsed time from trial onset (bin centre after downsampling).', ...
    'Units', 's');

meta.HapticFeedback = struct( ...
    'LongName', 'Haptic Feedback State', ...
    'Description', ['Binary state of the vibration motor at each downsampled time point ' ...
        '(software-triggered). Values are 0 (motor off) or 1 (motor on). ' ...
        'Downsampling uses nearest-neighbour decimation: for each output time bin, ' ...
        'the raw sample closest to the bin centre is selected, preserving the ' ...
        'binary nature of the signal.'], ...
    'Units', 'binary (0 or 1)');

% Source data
meta.SourceData = struct( ...
    'Description', 'Raw CSV trial recordings from perceptual crossing experiment.', ...
    'FileFormat', '.csv', ...
    'TrialDuration', trial_duration, ...
    'TrialDurationUnit', 's', ...
    'TrialsPerDyad', num_trials_exp, ...
    'ParticipantsPerDyad', num_participants, ...
    'FolderStructure', 'pceXXYYMMDD/trials/pair_XX_trial_Y.csv', ...
    'RawColumns', {{'timestamp', 'motor_0_vibrate_software', 'motor_1_vibrate_software'}}, ...
    'EstimatedRawSamplingRate_Hz', round(mean_fs_raw, 1), ...
    'EstimatedRawSamplingRate_Std', round(std_fs_raw, 2));

% Extraction method
meta.ExtractionMethod = struct( ...
    'Description', 'Nearest-neighbour decimation of binary haptic feedback signal.', ...
    'Steps', {{...
        '1. Discover pce* folders and parse DyadID from folder name', ...
        '2. Read each trial CSV (pair_XX_trial_Y.csv)', ...
        '3. Extract timestamp, motor_0_vibrate_software, motor_1_vibrate_software', ...
        '4. Downsample by nearest-neighbour: select raw sample closest to each bin centre', ...
        ['5. Target rate: ' num2str(fs_target) ' Hz (bin width = ' ...
         sprintf('%.0f', 1000/fs_target) ' ms)'], ...
        '6. Store one time series per participant per trial', ...
        '7. Sort by DyadID, ParticipantID, TrialNum and export'}}, ...
    'DownsamplingMethod', 'Nearest-neighbour decimation (preserves binary values)', ...
    'TargetSamplingRate', fs_target, ...
    'TargetSamplingRateUnit', 'Hz');

% Data summary
meta.DataSummary = struct( ...
    'NumDyads', length(dyad_ids_found), ...
    'NumTrialFilesProcessed', total_trials_found, ...
    'NumTrialFilesSkipped', total_trials_skipped, ...
    'NumSegmentsExtracted', row_count, ...
    'TotalSamples', total_samples, ...
    'MeanHapticFeedback', mean_haptic, ...
    'StdHapticFeedback', std_haptic, ...
    'MedianHapticFeedback', median_haptic);

meta.GeneratedBy = struct( ...
    'Name', 'preprocessHaptics.m', ...
    'Description', 'MATLAB preprocessing script for PCE haptic feedback data.', ...
    'GenerationDateTime', timestamp_str);

%% Write JSON
output_json = 'HapticFeedback.json';
json_str = jsonencode(meta);
json_str = prettify_json(json_str);
fid = fopen(output_json, 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

fprintf('Writing %s ... done.\n', output_json);

%% Final Summary
fprintf('\n==========================================================\n');
fprintf('  SUMMARY\n');
fprintf('==========================================================\n');
fprintf('  Dyads:                    %d\n', length(dyad_ids_found));
fprintf('  Trial files processed:    %d\n', total_trials_found);
fprintf('  Trial files skipped:      %d\n', total_trials_skipped);
fprintf('  Segments extracted:       %d\n', row_count);
fprintf('  Total samples (CSV):      %d\n', total_samples);
fprintf('  Target sampling rate:     %d Hz\n', fs_target);
fprintf('  Estimated raw rate:       %.1f +/- %.1f Hz\n', mean_fs_raw, std_fs_raw);
fprintf('  Mean haptic feedback:     %.4f\n', mean_haptic);
fprintf('  Output files created:\n');
fprintf('    %s\n', output_csv);
fprintf('    %s\n', output_json);
fprintf('==========================================================\n');
fprintf('Done.\n');

%% ========================================================================
%  LOCAL FUNCTION: prettify_json
%  ========================================================================

function pretty = prettify_json(json_str)
    indent = 0;
    indent_str = '    ';
    pretty = '';
    in_string = false;
    i = 1;
    n = length(json_str);

    while i <= n
        ch = json_str(i);

        if ch == '"' && (i == 1 || json_str(i-1) ~= '\')
            in_string = ~in_string;
            pretty = [pretty, ch]; %#ok<AGROW>
            i = i + 1;
            continue;
        end

        if in_string
            pretty = [pretty, ch]; %#ok<AGROW>
            i = i + 1;
            continue;
        end

        switch ch
            case {'{', '['}
                pretty = [pretty, ch, newline]; %#ok<AGROW>
                indent = indent + 1;
                pretty = [pretty, repmat(indent_str, 1, indent)]; %#ok<AGROW>

            case {'}', ']'}
                pretty = [pretty, newline]; %#ok<AGROW>
                indent = indent - 1;
                pretty = [pretty, repmat(indent_str, 1, indent), ch]; %#ok<AGROW>

            case ','
                pretty = [pretty, ',', newline]; %#ok<AGROW>
                pretty = [pretty, repmat(indent_str, 1, indent)]; %#ok<AGROW>

            case ':'
                pretty = [pretty, ': ']; %#ok<AGROW>

            otherwise
                if ~isspace(ch)
                    pretty = [pretty, ch]; %#ok<AGROW>
                end
        end
        i = i + 1;
    end
end