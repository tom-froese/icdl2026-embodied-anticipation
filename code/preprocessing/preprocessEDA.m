%% preprocessEDA.m
% =========================================================================
% EDA Preprocessing and CSV Export Script (with JSON Metadata Sidecars)
% =========================================================================
%
% PURPOSE:
%   Reads raw EDA data from the perceptual crossing experiment (PCE) folder
%   structure, applies standard preprocessing following published guidelines
%   (Boucsein et al., 2012; Braithwaite et al., 2013; Society for
%   Psychophysiological Research recommendations), and exports consolidated
%   CSV files for task and rest conditions with BIDS-inspired JSON metadata
%   sidecar files.
%
% PREPROCESSING PIPELINE:
%   1. Low-pass filter: 4th-order Butterworth, 5 Hz cutoff (zero-phase)
%      - Removes high-frequency noise, EMG artifacts, and mains interference
%      - 5 Hz is the standard recommendation for EDA (Boucsein et al., 2012;
%        Bach et al., 2013; Braithwaite et al., 2013)
%   2. Downsample: 1000 Hz -> 25 Hz (factor of 40)
%      - 25 Hz preserves all EDA-relevant information (phasic SCRs have
%        rise times >= 1 s; Nyquist for 5 Hz content requires >= 10 Hz)
%      - Reduces file size by 40x while retaining full signal fidelity
%   3. Artifact flagging: physiological range check (0.01-100 uS) and
%      excessive gradient detection (> 10 uS/s)
%      - Flagged samples marked in a separate column (not removed)
%      - Allows downstream decisions about exclusion
%
% INPUT:
%   Run this script from the root EDA folder containing pce* subfolders.
%   Each subfolder contains .mat files named:
%     pceXX_P{1,2}_Trial{1..18}.mat  (variable: pce_P{1,2}_Trial)
%     pceXX_P{1,2}_Rest{1..4}.mat    (variable: pce_P{1,2}_Rest)
%
% OUTPUT:
%   EDA_Task_Preprocessed.csv   - All task trials (long format)
%   EDA_Task_Preprocessed.json  - BIDS-style metadata sidecar for task CSV
%   EDA_Rest_Preprocessed.csv   - All rest periods (long format)
%   EDA_Rest_Preprocessed.json  - BIDS-style metadata sidecar for rest CSV
%
% REFERENCES:
%   Boucsein, W. et al. (2012). Publication recommendations for
%     electrodermal measurements. Psychophysiology, 49, 1017-1034.
%   Braithwaite, J.J. et al. (2013). A guide for analysing electrodermal
%     activity (EDA) & skin conductance responses (SCRs) for psychological
%     experiments. Psychophysiology, 49, 1017-1034.
%   Bach, D.R. et al. (2013). An improved algorithm for model-based
%     analysis of evoked skin conductance responses. Biological Psychology,
%     94, 490-497.
%   BIDS Specification v1.10 - Physiological recordings:
%     https://bids-specification.readthedocs.io/en/stable/
%
% AUTHOR: Preprocessing script for PCE experiment
% DATE:   February 2026
% =========================================================================

%% Parameters
fs_original   = 1000;       % Original sampling rate (Hz)
fs_target     = 25;         % Target sampling rate after downsampling (Hz)
downsample_factor = fs_original / fs_target;  % = 40

% Low-pass filter parameters
filter_order  = 4;          % Butterworth filter order (4th order = 24 dB/oct)
filter_cutoff = 5;          % Cutoff frequency (Hz) - standard for EDA

% Artifact detection thresholds
eda_min       = 0.01;       % Minimum plausible EDA (uS)
eda_max       = 100;        % Maximum plausible EDA (uS)
gradient_max  = 10;         % Maximum plausible gradient (uS/s)

% Trial/rest parameters
num_trials    = 18;
num_rests     = 4;
task_duration = 60;         % seconds
rest_duration = 180;        % seconds
task_samples  = task_duration * fs_original;    % 60000
rest_samples  = rest_duration * fs_original;    % 180000

%% Design low-pass filter (once, for original sampling rate)
[b_lp, a_lp] = butter(filter_order, filter_cutoff / (fs_original / 2), 'low');

fprintf('==========================================================\n');
fprintf('  EDA Preprocessing Script for PCE Experiment\n');
fprintf('==========================================================\n');
fprintf('  Low-pass filter: %d-order Butterworth, %.0f Hz cutoff\n', filter_order, filter_cutoff);
fprintf('  Downsampling:    %d Hz -> %d Hz (factor %d)\n', fs_original, fs_target, downsample_factor);
fprintf('  Artifact range:  [%.2f, %.0f] uS, gradient < %.0f uS/s\n', eda_min, eda_max, gradient_max);
fprintf('==========================================================\n\n');

%% Discover experiment folders
folders = dir('pce*');
folders = folders([folders.isdir]);  % Keep only directories
fprintf('Found %d experiment folders.\n\n', length(folders));

if isempty(folders)
    error('No pce* folders found. Please run this script from the root EDA data directory.');
end

%% Preallocate output cell arrays
task_rows = {};
rest_rows = {};

task_count = 0;
rest_count = 0;
task_skipped = 0;
rest_skipped = 0;

% Track per-segment summary statistics for metadata
task_eda_means = [];
rest_eda_means = [];
dyad_ids_found = [];

%% Main loop
for f = 1:length(folders)
    folder_name = folders(f).name;
    
    % Extract dyad/experiment number (first 2 digits after 'pce')
    exp_num_str = folder_name(4:5);
    dyad_id = str2double(exp_num_str);
    
    fprintf('Processing folder: %s (Dyad %02d)\n', folder_name, dyad_id);
    
    if ~ismember(dyad_id, dyad_ids_found)
        dyad_ids_found = [dyad_ids_found, dyad_id]; %#ok<AGROW>
    end
    
    for p = 1:2
        participant_id = sprintf('P%d', p);
        
        % ----- TASK TRIALS -----
        for t = 1:num_trials
            filename = fullfile(folder_name, ...
                sprintf('pce%s_%s_Trial%d.mat', exp_num_str, participant_id, t));
            
            if ~isfile(filename)
                task_skipped = task_skipped + 1;
                continue;
            end
            
            try
                data = load(filename);
                var_name = sprintf('pce_%s_Trial', participant_id);
                raw_eda = double(data.(var_name));
                raw_eda = raw_eda(:)';  % Ensure row vector
                
                % Validate expected length
                if length(raw_eda) ~= task_samples
                    warning('  %s: unexpected length %d (expected %d). Skipping.', ...
                        filename, length(raw_eda), task_samples);
                    task_skipped = task_skipped + 1;
                    continue;
                end
                
                % --- Preprocessing ---
                % 1. Zero-phase low-pass filter
                eda_filtered = filtfilt(b_lp, a_lp, raw_eda);
                
                % 2. Downsample (take every Nth sample after filtering)
                eda_ds = eda_filtered(1:downsample_factor:end);
                n_ds = length(eda_ds);
                
                % 3. Construct time vector (downsampled)
                time_ds = (0:n_ds-1) / fs_target;
                
                % 4. Artifact detection
                range_flag = (eda_ds < eda_min) | (eda_ds > eda_max);
                gradient_eda = [0, abs(diff(eda_ds)) * fs_target];
                gradient_flag = gradient_eda > gradient_max;
                artifact_flag = double(range_flag | gradient_flag);
                
                % Store row
                task_count = task_count + 1;
                task_rows{task_count} = {dyad_id, p, t, time_ds, eda_ds, artifact_flag};
                
                % Track statistics
                task_eda_means = [task_eda_means, mean(eda_ds)]; %#ok<AGROW>
                
            catch ME
                warning('  Error loading %s: %s. Skipping.', filename, ME.message);
                task_skipped = task_skipped + 1;
            end
        end
        
        % ----- REST PERIODS -----
        for r = 1:num_rests
            filename = fullfile(folder_name, ...
                sprintf('pce%s_%s_Rest%d.mat', exp_num_str, participant_id, r));
            
            if ~isfile(filename)
                rest_skipped = rest_skipped + 1;
                continue;
            end
            
            try
                data = load(filename);
                var_name = sprintf('pce_%s_Rest', participant_id);
                raw_eda = double(data.(var_name));
                raw_eda = raw_eda(:)';  % Ensure row vector
                
                % Validate expected length
                if length(raw_eda) ~= rest_samples
                    warning('  %s: unexpected length %d (expected %d). Skipping.', ...
                        filename, length(raw_eda), rest_samples);
                    rest_skipped = rest_skipped + 1;
                    continue;
                end
                
                % --- Preprocessing ---
                eda_filtered = filtfilt(b_lp, a_lp, raw_eda);
                eda_ds = eda_filtered(1:downsample_factor:end);
                n_ds = length(eda_ds);
                time_ds = (0:n_ds-1) / fs_target;
                
                range_flag = (eda_ds < eda_min) | (eda_ds > eda_max);
                gradient_eda = [0, abs(diff(eda_ds)) * fs_target];
                gradient_flag = gradient_eda > gradient_max;
                artifact_flag = double(range_flag | gradient_flag);
                
                rest_count = rest_count + 1;
                rest_rows{rest_count} = {dyad_id, p, r, time_ds, eda_ds, artifact_flag};
                
                rest_eda_means = [rest_eda_means, mean(eda_ds)]; %#ok<AGROW>
                
            catch ME
                warning('  Error loading %s: %s. Skipping.', filename, ME.message);
                rest_skipped = rest_skipped + 1;
            end
        end
    end
end

fprintf('\n==========================================================\n');
fprintf('  Data collection complete.\n');
fprintf('  Task trials:  %d loaded, %d skipped\n', task_count, task_skipped);
fprintf('  Rest periods: %d loaded, %d skipped\n', rest_count, rest_skipped);
fprintf('==========================================================\n\n');

%% Export Task CSV
fprintf('Writing EDA_Task_Preprocessed.csv ...\n');

fid = fopen('EDA_Task_Preprocessed.csv', 'w');
fprintf(fid, 'DyadID,ParticipantID,TrialNum,Time_s,EDA_uS,ArtifactFlag\n');

for i = 1:task_count
    row = task_rows{i};
    dyad_id = row{1};
    part_id = row{2};
    trial_num = row{3};
    time_vec = row{4};
    eda_vec  = row{5};
    flag_vec = row{6};
    
    for s = 1:length(time_vec)
        fprintf(fid, '%d,%d,%d,%.4f,%.6f,%d\n', ...
            dyad_id, part_id, trial_num, time_vec(s), eda_vec(s), flag_vec(s));
    end
end

fclose(fid);

task_info = dir('EDA_Task_Preprocessed.csv');
fprintf('  -> %s (%.1f MB)\n', task_info.name, task_info.bytes / 1e6);

%% Export Rest CSV
fprintf('Writing EDA_Rest_Preprocessed.csv ...\n');

fid = fopen('EDA_Rest_Preprocessed.csv', 'w');
fprintf(fid, 'DyadID,ParticipantID,RestNum,Time_s,EDA_uS,ArtifactFlag\n');

for i = 1:rest_count
    row = rest_rows{i};
    dyad_id = row{1};
    part_id = row{2};
    rest_num = row{3};
    time_vec = row{4};
    eda_vec  = row{5};
    flag_vec = row{6};
    
    for s = 1:length(time_vec)
        fprintf(fid, '%d,%d,%d,%.4f,%.6f,%d\n', ...
            dyad_id, part_id, rest_num, time_vec(s), eda_vec(s), flag_vec(s));
    end
end

fclose(fid);

rest_info = dir('EDA_Rest_Preprocessed.csv');
fprintf('  -> %s (%.1f MB)\n', rest_info.name, rest_info.bytes / 1e6);

%% Count artifacts for summary
total_task_samples = 0;
total_task_artifacts = 0;
for i = 1:task_count
    flags = task_rows{i}{6};
    total_task_samples = total_task_samples + length(flags);
    total_task_artifacts = total_task_artifacts + sum(flags);
end

total_rest_samples = 0;
total_rest_artifacts = 0;
for i = 1:rest_count
    flags = rest_rows{i}{6};
    total_rest_samples = total_rest_samples + length(flags);
    total_rest_artifacts = total_rest_artifacts + sum(flags);
end

%% ========================================================================
%  Generate JSON Metadata Sidecars (BIDS-inspired)
%  ========================================================================
%  Following BIDS Specification conventions for physiological recordings:
%    - SamplingFrequency, StartTime, Columns  (REQUIRED per BIDS physio)
%    - Per-column Description and Units        (RECOMMENDED per BIDS)
%    - CamelCase keys                          (RECOMMENDED per BIDS)
%    - GeneratedBy provenance block            (BIDS derivatives convention)
%    - Full preprocessing parameter documentation
%  ========================================================================

timestamp_str = datestr(now, 'yyyy-mm-ddTHH:MM:SS');

% ---- Task JSON ----
task_meta = struct();

% Dataset-level fields
task_meta.Name = 'Perceptual Crossing Experiment - EDA Task Condition (Preprocessed)';
task_meta.Description = 'Preprocessed electrodermal activity (EDA) data from task trials of a perceptual coupling experiment. Each trial is 60 seconds. Data are organized in long format with one row per sample per trial.';

% BIDS Physiological Recording fields (REQUIRED)
task_meta.SamplingFrequency = fs_target;
task_meta.SamplingFrequencyUnit = 'Hz';
task_meta.StartTime = 0;
task_meta.Columns = {'DyadID', 'ParticipantID', 'TrialNum', 'Time_s', 'EDA_uS', 'ArtifactFlag'};

% Per-column documentation (BIDS RECOMMENDED)
task_meta.DyadID = struct( ...
    'LongName', 'Dyad Identifier', ...
    'Description', 'Numeric identifier for the experimental dyad (pair of participants). Extracted from the first two digits of the pce folder name.');

task_meta.ParticipantID = struct( ...
    'LongName', 'Participant Identifier', ...
    'Description', 'Participant number within the dyad (1 or 2).');

task_meta.TrialNum = struct( ...
    'LongName', 'Trial Number', ...
    'Description', 'Task trial number within the session (1-18).');

task_meta.Time_s = struct( ...
    'LongName', 'Time', ...
    'Description', 'Elapsed time from the start of the trial.', ...
    'Units', 's');

task_meta.EDA_uS = struct( ...
    'LongName', 'Electrodermal Activity', ...
    'Description', 'Skin conductance level after low-pass filtering and downsampling. No baseline correction or tonic/phasic decomposition applied.', ...
    'Units', 'uS');

task_meta.ArtifactFlag = struct( ...
    'LongName', 'Artifact Flag', ...
    'Description', 'Binary flag indicating samples that failed physiological plausibility checks. 0 = clean, 1 = flagged. Flagged samples are retained in the data for transparency; downstream exclusion decisions are left to the analyst.', ...
    'Levels', struct('x0', 'Clean sample', 'x1', 'Artifact detected'));

% Source data information
task_meta.SourceData = struct( ...
    'Description', 'Raw EDA recordings from perceptual coupling experiment.', ...
    'OriginalSamplingFrequency', fs_original, ...
    'OriginalSamplingFrequencyUnit', 'Hz', ...
    'FileFormat', '.mat (MATLAB)', ...
    'TrialDuration', task_duration, ...
    'TrialDurationUnit', 's', ...
    'TrialsPerParticipant', num_trials, ...
    'ParticipantsPerDyad', 2, ...
    'RecordingDevice', 'Standard skin conductance electrodes');

% Preprocessing provenance (BIDS GeneratedBy convention)
task_meta.Preprocessing = struct( ...
    'Description', 'Standard EDA preprocessing pipeline following SPR publication recommendations.', ...
    'Steps', {{'1. Zero-phase low-pass Butterworth filter (filtfilt)', ...
               '2. Downsampling by decimation (every Nth sample)', ...
               '3. Artifact detection via range and gradient checks'}}, ...
    'LowPassFilter', struct( ...
        'Type', 'Butterworth', ...
        'Order', filter_order, ...
        'CutoffFrequency', filter_cutoff, ...
        'CutoffFrequencyUnit', 'Hz', ...
        'Implementation', 'Zero-phase (MATLAB filtfilt)', ...
        'Rationale', 'Standard recommendation for EDA: removes high-frequency noise, EMG artifacts, and mains interference while preserving all physiologically relevant EDA dynamics (Boucsein et al., 2012; Bach et al., 2013).'), ...
    'Downsampling', struct( ...
        'OriginalRate', fs_original, ...
        'TargetRate', fs_target, ...
        'Factor', downsample_factor, ...
        'Method', 'Decimation (every Nth sample after anti-alias filtering)', ...
        'Rationale', '25 Hz exceeds Nyquist for 5 Hz low-pass filtered signal. Preserves full signal fidelity with 40x data reduction.'), ...
    'ArtifactDetection', struct( ...
        'Method', 'Combined range and gradient check', ...
        'RangeMin', eda_min, ...
        'RangeMax', eda_max, ...
        'RangeUnit', 'uS', ...
        'GradientThreshold', gradient_max, ...
        'GradientUnit', 'uS/s', ...
        'Action', 'Flag only (no removal or interpolation)'));

% References
task_meta.References = { ...
    'Boucsein, W., Fowles, D.C., Grimnes, S., Ben-Shakhar, G., Roth, W.T., Dawson, M.E., & Filion, D.L. (2012). Publication recommendations for electrodermal measurements. Psychophysiology, 49(8), 1017-1034.', ...
    'Braithwaite, J.J., Watson, D.G., Jones, R., & Rowe, M. (2013). A guide for analysing electrodermal activity (EDA) & skin conductance responses (SCRs) for psychological experiments. Psychophysiology, 49, 1017-1034.', ...
    'Bach, D.R., Flandin, G., Friston, K.J., & Dolan, R.J. (2013). An improved algorithm for model-based analysis of evoked skin conductance responses. Biological Psychology, 94(3), 490-497.'};

% Data summary statistics
task_meta.DataSummary = struct( ...
    'NumDyads', length(dyad_ids_found), ...
    'NumTrialsLoaded', task_count, ...
    'NumTrialsSkipped', task_skipped, ...
    'SamplesPerTrial', task_duration * fs_target, ...
    'TotalSamples', total_task_samples, ...
    'TotalArtifactsFlagged', total_task_artifacts, ...
    'ArtifactRate', total_task_artifacts / max(total_task_samples, 1), ...
    'MeanEDA_uS', mean(task_eda_means), ...
    'StdEDA_uS', std(task_eda_means));

% File provenance
task_meta.GeneratedBy = struct( ...
    'Name', 'preprocessEDA.m', ...
    'Description', 'MATLAB preprocessing script for PCE EDA data.', ...
    'GenerationDateTime', timestamp_str);

% Write Task JSON
task_json_str = jsonencode(task_meta);
task_json_str = prettify_json(task_json_str);
fid = fopen('EDA_Task_Preprocessed.json', 'w');
fprintf(fid, '%s', task_json_str);
fclose(fid);
fprintf('Writing EDA_Task_Preprocessed.json ... done.\n');

% ---- Rest JSON ----
rest_meta = struct();

rest_meta.Name = 'Perceptual Crossing Experiment - EDA Rest Condition (Preprocessed)';
rest_meta.Description = 'Preprocessed electrodermal activity (EDA) data from rest periods of a perceptual coupling experiment. Each rest period is 180 seconds. Data are organized in long format with one row per sample per period.';

rest_meta.SamplingFrequency = fs_target;
rest_meta.SamplingFrequencyUnit = 'Hz';
rest_meta.StartTime = 0;
rest_meta.Columns = {'DyadID', 'ParticipantID', 'RestNum', 'Time_s', 'EDA_uS', 'ArtifactFlag'};

rest_meta.DyadID = struct( ...
    'LongName', 'Dyad Identifier', ...
    'Description', 'Numeric identifier for the experimental dyad (pair of participants). Extracted from the first two digits of the pce folder name.');

rest_meta.ParticipantID = struct( ...
    'LongName', 'Participant Identifier', ...
    'Description', 'Participant number within the dyad (1 or 2).');

rest_meta.RestNum = struct( ...
    'LongName', 'Rest Period Number', ...
    'Description', 'Rest period number within the session (1-4).');

rest_meta.Time_s = struct( ...
    'LongName', 'Time', ...
    'Description', 'Elapsed time from the start of the rest period.', ...
    'Units', 's');

rest_meta.EDA_uS = struct( ...
    'LongName', 'Electrodermal Activity', ...
    'Description', 'Skin conductance level after low-pass filtering and downsampling. No baseline correction or tonic/phasic decomposition applied.', ...
    'Units', 'uS');

rest_meta.ArtifactFlag = struct( ...
    'LongName', 'Artifact Flag', ...
    'Description', 'Binary flag indicating samples that failed physiological plausibility checks. 0 = clean, 1 = flagged. Flagged samples are retained in the data for transparency; downstream exclusion decisions are left to the analyst.', ...
    'Levels', struct('x0', 'Clean sample', 'x1', 'Artifact detected'));

rest_meta.SourceData = struct( ...
    'Description', 'Raw EDA recordings from perceptual coupling experiment.', ...
    'OriginalSamplingFrequency', fs_original, ...
    'OriginalSamplingFrequencyUnit', 'Hz', ...
    'FileFormat', '.mat (MATLAB)', ...
    'RestDuration', rest_duration, ...
    'RestDurationUnit', 's', ...
    'RestPeriodsPerParticipant', num_rests, ...
    'ParticipantsPerDyad', 2, ...
    'RecordingDevice', 'Standard skin conductance electrodes');

rest_meta.Preprocessing = struct( ...
    'Description', 'Standard EDA preprocessing pipeline following SPR publication recommendations.', ...
    'Steps', {{'1. Zero-phase low-pass Butterworth filter (filtfilt)', ...
               '2. Downsampling by decimation (every Nth sample)', ...
               '3. Artifact detection via range and gradient checks'}}, ...
    'LowPassFilter', struct( ...
        'Type', 'Butterworth', ...
        'Order', filter_order, ...
        'CutoffFrequency', filter_cutoff, ...
        'CutoffFrequencyUnit', 'Hz', ...
        'Implementation', 'Zero-phase (MATLAB filtfilt)', ...
        'Rationale', 'Standard recommendation for EDA: removes high-frequency noise, EMG artifacts, and mains interference while preserving all physiologically relevant EDA dynamics (Boucsein et al., 2012; Bach et al., 2013).'), ...
    'Downsampling', struct( ...
        'OriginalRate', fs_original, ...
        'TargetRate', fs_target, ...
        'Factor', downsample_factor, ...
        'Method', 'Decimation (every Nth sample after anti-alias filtering)', ...
        'Rationale', '25 Hz exceeds Nyquist for 5 Hz low-pass filtered signal. Preserves full signal fidelity with 40x data reduction.'), ...
    'ArtifactDetection', struct( ...
        'Method', 'Combined range and gradient check', ...
        'RangeMin', eda_min, ...
        'RangeMax', eda_max, ...
        'RangeUnit', 'uS', ...
        'GradientThreshold', gradient_max, ...
        'GradientUnit', 'uS/s', ...
        'Action', 'Flag only (no removal or interpolation)'));

rest_meta.References = { ...
    'Boucsein, W., Fowles, D.C., Grimnes, S., Ben-Shakhar, G., Roth, W.T., Dawson, M.E., & Filion, D.L. (2012). Publication recommendations for electrodermal measurements. Psychophysiology, 49(8), 1017-1034.', ...
    'Braithwaite, J.J., Watson, D.G., Jones, R., & Rowe, M. (2013). A guide for analysing electrodermal activity (EDA) & skin conductance responses (SCRs) for psychological experiments. Psychophysiology, 49, 1017-1034.', ...
    'Bach, D.R., Flandin, G., Friston, K.J., & Dolan, R.J. (2013). An improved algorithm for model-based analysis of evoked skin conductance responses. Biological Psychology, 94(3), 490-497.'};

rest_meta.DataSummary = struct( ...
    'NumDyads', length(dyad_ids_found), ...
    'NumPeriodsLoaded', rest_count, ...
    'NumPeriodsSkipped', rest_skipped, ...
    'SamplesPerPeriod', rest_duration * fs_target, ...
    'TotalSamples', total_rest_samples, ...
    'TotalArtifactsFlagged', total_rest_artifacts, ...
    'ArtifactRate', total_rest_artifacts / max(total_rest_samples, 1), ...
    'MeanEDA_uS', mean(rest_eda_means), ...
    'StdEDA_uS', std(rest_eda_means));

rest_meta.GeneratedBy = struct( ...
    'Name', 'preprocessEDA.m', ...
    'Description', 'MATLAB preprocessing script for PCE EDA data.', ...
    'GenerationDateTime', timestamp_str);

% Write Rest JSON
rest_json_str = jsonencode(rest_meta);
rest_json_str = prettify_json(rest_json_str);
fid = fopen('EDA_Rest_Preprocessed.json', 'w');
fprintf(fid, '%s', rest_json_str);
fclose(fid);
fprintf('Writing EDA_Rest_Preprocessed.json ... done.\n');

%% Final Summary
fprintf('\n==========================================================\n');
fprintf('  SUMMARY\n');
fprintf('==========================================================\n');
fprintf('  Task:  %d trials, %d samples total, %d flagged (%.2f%%)\n', ...
    task_count, total_task_samples, total_task_artifacts, ...
    100 * total_task_artifacts / max(total_task_samples, 1));
fprintf('  Rest:  %d periods, %d samples total, %d flagged (%.2f%%)\n', ...
    rest_count, total_rest_samples, total_rest_artifacts, ...
    100 * total_rest_artifacts / max(total_rest_samples, 1));
fprintf('\n  Output sampling rate: %d Hz\n', fs_target);
fprintf('  Task samples per trial: %d (%.0f s at %d Hz)\n', ...
    task_duration * fs_target, task_duration, fs_target);
fprintf('  Rest samples per period: %d (%.0f s at %d Hz)\n', ...
    rest_duration * fs_target, rest_duration, fs_target);

raw_task_bytes = task_count * task_samples * 4;
raw_rest_bytes = rest_count * rest_samples * 4;
fprintf('\n  Estimated raw data size:  %.1f MB\n', (raw_task_bytes + raw_rest_bytes) / 1e6);
fprintf('  Output CSV total size:    %.1f MB\n', ...
    (task_info.bytes + rest_info.bytes) / 1e6);
fprintf('\n  Output files:\n');
fprintf('    EDA_Task_Preprocessed.csv   (data)\n');
fprintf('    EDA_Task_Preprocessed.json  (metadata sidecar)\n');
fprintf('    EDA_Rest_Preprocessed.csv   (data)\n');
fprintf('    EDA_Rest_Preprocessed.json  (metadata sidecar)\n');
fprintf('==========================================================\n');
fprintf('Done.\n');


%% ========================================================================
%  LOCAL FUNCTION: prettify_json
%  ========================================================================
%  MATLAB's jsonencode produces a single-line string. This function adds
%  indentation and line breaks for human readability, following the BIDS
%  convention of pretty-printed JSON sidecars.
%  ========================================================================

function pretty = prettify_json(json_str)
    indent = 0;
    indent_str = '    ';  % 4 spaces per level
    pretty = '';
    in_string = false;
    i = 1;
    n = length(json_str);
    
    while i <= n
        ch = json_str(i);
        
        % Track whether we are inside a JSON string value
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