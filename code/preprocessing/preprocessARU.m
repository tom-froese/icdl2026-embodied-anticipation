%% preprocessARU.m
% =========================================================================
% Respiration Preprocessing and CSV/JSON Export
% =========================================================================
%
% PURPOSE:
%   Reads raw respiration data (ARU) from the perceptual crossing experiment
%   (PCE) folder structure, applies standard preprocessing following
%   published guidelines for respiratory inductance plethysmography and
%   strain-gauge belt recordings (Lorig, 2012; Kreibig, 2010; Massaroni
%   et al., 2019), and exports consolidated CSV files for task and rest
%   conditions with BIDS-inspired JSON metadata sidecar files.
%
% PREPROCESSING PIPELINE:
%   1. Low-pass filter: 4th-order Butterworth, 1.0 Hz cutoff (zero-phase)
%      - Removes movement artifacts, cardiac contamination, and
%        high-frequency noise (1.0 Hz = 60 br/min, well above
%        physiological maximum of ~40 br/min)
%      - NO high-pass filter applied: baseline drift and slow trends
%        are intentionally preserved as they carry the P(1) dynamics
%        of theoretical interest (task-dependent respiratory flux)
%   2. Downsample: 1000 Hz -> 25 Hz (factor of 40)
%      - 25 Hz preserves all respiratory information (Nyquist for 1 Hz
%        content requires >= 2 Hz; 25 Hz provides generous oversampling
%        for accurate peak/trough detection)
%   3. Signal inversion: multiply by -1
%      - Brain Products respiration belt records lung extension as
%        negative-going deflections; inversion yields the physiological
%        convention where inspiration = positive deflection (peak),
%        expiration = negative deflection (trough)
%   4. Artifact detection:
%      - Amplitude clipping: abs(ARU) exceeding 3x the recording's IQR
%        above/below the median (adaptive per-recording)
%      - Flatline detection: segments where the gradient is near-zero
%        for >= 2 seconds (belt displacement or loss of contact)
%
% INPUT:
%   Run this script from the root ARU folder containing pce* subfolders.
%   Each subfolder contains .mat files named:
%     pceXX_P{1,2}_Trial{1..18}.mat  (variable: pce_P{1,2}_Trial)
%     pceXX_P{1,2}_Rest{1..4}.mat    (variable: pce_P{1,2}_Rest)
%
% OUTPUT:
%   ARU_Task_Preprocessed.csv   - Preprocessed continuous signal (task)
%   ARU_Task_Preprocessed.json  - BIDS-style metadata sidecar
%   ARU_Rest_Preprocessed.csv   - Preprocessed continuous signal (rest)
%   ARU_Rest_Preprocessed.json  - BIDS-style metadata sidecar
%
% REFERENCES:
%   Lorig, T.S. (2012). Respiratory measurement: Techniques. In J.T.
%     Cacioppo, L.G. Tassinary, & G.G. Berntson (Eds.), Handbook of
%     Psychophysiology (3rd ed., pp. 231-238). Cambridge University Press.
%   Kreibig, S.D. (2010). Autonomic nervous system activity in emotion:
%     A review. Biological Psychology, 84(3), 394-421.
%   Massaroni, C. et al. (2019). Contact-based methods for measuring
%     respiratory rate. Sensors, 19(4), 908.
%   BIDS Specification v1.10 - Physiological recordings:
%     https://bids-specification.readthedocs.io/en/stable/
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   February 2026
% =========================================================================

%% Parameters
fs_original   = 1000;       % Original sampling rate (Hz)
fs_target     = 25;         % Target sampling rate after downsampling (Hz)
downsample_factor = fs_original / fs_target;  % = 40

% Low-pass filter parameters
filter_order  = 4;          % Butterworth filter order (4th order = 24 dB/oct)
lp_cutoff     = 1.0;        % Low-pass cutoff (Hz) - 60 br/min ceiling
                             % NOTE: No high-pass applied. Baseline drift and
                             % slow trends are PRESERVED intentionally, as
                             % they carry the P(1) dynamics of interest.

% Artifact detection parameters (adaptive thresholds)
clip_iqr_factor = 3.0;      % Flag amplitudes > median +/- factor*IQR
flatline_dur    = 2.0;      % Minimum flatline duration to flag (s)
flatline_grad   = 0.5;      % Gradient threshold for flatline (ARU/s)

% Trial/rest parameters
num_trials    = 18;
num_rests     = 4;
task_duration = 60;         % seconds
rest_duration = 180;        % seconds
task_samples  = task_duration * fs_original;    % 60000
rest_samples  = rest_duration * fs_original;    % 180000

%% Design low-pass filter (once, for original sampling rate)
%  Low-pass only: preserves baseline drift and slow trends that carry
%  P(1) dynamics. No high-pass is applied.
[b_lp, a_lp] = butter(filter_order, lp_cutoff / (fs_original / 2), 'low');

fprintf('==========================================================\n');
fprintf('  Respiration Preprocessing Script for PCE Experiment\n');
fprintf('==========================================================\n');
fprintf('  Low-pass filter: %d-order Butterworth, %.1f Hz cutoff\n', ...
    filter_order, lp_cutoff);
fprintf('  High-pass:       none (slow trends preserved for P(1) analysis)\n');
fprintf('  Downsampling:    %d Hz -> %d Hz (factor %d)\n', ...
    fs_original, fs_target, downsample_factor);
fprintf('  Signal inversion: yes (convention: inspiration = positive)\n');
fprintf('==========================================================\n\n');

%% Discover experiment folders
folders = dir('pce*');
folders = folders([folders.isdir]);
fprintf('Found %d experiment folders.\n\n', length(folders));

if isempty(folders)
    error('No pce* folders found. Please run this script from the root ARU data directory.');
end

%% Preallocate output cell arrays
task_signal_rows = {};
rest_signal_rows = {};

task_count = 0;
rest_count = 0;
task_skipped = 0;
rest_skipped = 0;

% Track per-segment summary statistics
task_amp_means = [];
rest_amp_means = [];
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
                raw_resp = double(data.(var_name));
                raw_resp = raw_resp(:)';  % Ensure row vector
                
                % Validate expected length
                if length(raw_resp) ~= task_samples
                    warning('  %s: unexpected length %d (expected %d). Skipping.', ...
                        filename, length(raw_resp), task_samples);
                    task_skipped = task_skipped + 1;
                    continue;
                end
                
                % --- Preprocessing ---
                % 1. Zero-phase low-pass filter (no high-pass; preserves trends)
                resp_filtered = filtfilt(b_lp, a_lp, raw_resp);
                
                % 2. Signal inversion (Brain Products convention -> physio convention)
                resp_filtered = -resp_filtered;
                
                % 3. Downsample (take every Nth sample after filtering)
                resp_ds = resp_filtered(1:downsample_factor:end);
                n_ds = length(resp_ds);
                
                % 4. Construct time vector (downsampled)
                time_ds = (0:n_ds-1) / fs_target;
                
                % 5. Artifact detection (adaptive per-recording)
                resp_median = median(resp_ds);
                resp_iqr = iqr(resp_ds);
                clip_upper = resp_median + clip_iqr_factor * resp_iqr;
                clip_lower = resp_median - clip_iqr_factor * resp_iqr;
                clip_flag = (resp_ds > clip_upper) | (resp_ds < clip_lower);
                
                % Flatline detection
                gradient_resp = [0, abs(diff(resp_ds)) * fs_target];
                flatline_samples = round(flatline_dur * fs_target);
                flatline_flag = false(1, n_ds);
                run_count = 0;
                for s = 1:n_ds
                    if gradient_resp(s) < flatline_grad
                        run_count = run_count + 1;
                    else
                        if run_count >= flatline_samples
                            flatline_flag(s-run_count:s-1) = true;
                        end
                        run_count = 0;
                    end
                end
                if run_count >= flatline_samples
                    flatline_flag(n_ds-run_count+1:n_ds) = true;
                end
                
                artifact_flag = double(clip_flag | flatline_flag);
                
                % Store continuous signal row
                task_count = task_count + 1;
                task_signal_rows{task_count} = {dyad_id, p, t, time_ds, resp_ds, artifact_flag};
                
                % Track statistics
                task_amp_means = [task_amp_means, mean(abs(resp_ds))]; %#ok<AGROW>
                
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
                raw_resp = double(data.(var_name));
                raw_resp = raw_resp(:)';  % Ensure row vector
                
                % Validate expected length
                if length(raw_resp) ~= rest_samples
                    warning('  %s: unexpected length %d (expected %d). Skipping.', ...
                        filename, length(raw_resp), rest_samples);
                    rest_skipped = rest_skipped + 1;
                    continue;
                end
                
                % --- Preprocessing (same pipeline) ---
                resp_filtered = filtfilt(b_lp, a_lp, raw_resp);
                resp_filtered = -resp_filtered;
                resp_ds = resp_filtered(1:downsample_factor:end);
                n_ds = length(resp_ds);
                time_ds = (0:n_ds-1) / fs_target;
                
                % Artifact detection
                resp_median = median(resp_ds);
                resp_iqr = iqr(resp_ds);
                clip_upper = resp_median + clip_iqr_factor * resp_iqr;
                clip_lower = resp_median - clip_iqr_factor * resp_iqr;
                clip_flag = (resp_ds > clip_upper) | (resp_ds < clip_lower);
                
                gradient_resp = [0, abs(diff(resp_ds)) * fs_target];
                flatline_flag = false(1, n_ds);
                run_count = 0;
                for s = 1:n_ds
                    if gradient_resp(s) < flatline_grad
                        run_count = run_count + 1;
                    else
                        if run_count >= flatline_samples
                            flatline_flag(s-run_count:s-1) = true;
                        end
                        run_count = 0;
                    end
                end
                if run_count >= flatline_samples
                    flatline_flag(n_ds-run_count+1:n_ds) = true;
                end
                
                artifact_flag = double(clip_flag | flatline_flag);
                
                rest_count = rest_count + 1;
                rest_signal_rows{rest_count} = {dyad_id, p, r, time_ds, resp_ds, artifact_flag};
                
                rest_amp_means = [rest_amp_means, mean(abs(resp_ds))]; %#ok<AGROW>
                
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

%% ========================================================================
%  Export CSV: Continuous Signal
%  ========================================================================

% --- Task Signal CSV ---
fprintf('Writing ARU_Task_Preprocessed.csv ...\n');
fid = fopen('ARU_Task_Preprocessed.csv', 'w');
fprintf(fid, 'DyadID,ParticipantID,TrialNum,Time_s,RESP_ARU,ArtifactFlag\n');

for i = 1:task_count
    row = task_signal_rows{i};
    d = row{1}; pid = row{2}; tnum = row{3};
    tvec = row{4}; rvec = row{5}; fvec = row{6};
    for s = 1:length(tvec)
        fprintf(fid, '%d,%d,%d,%.4f,%.4f,%d\n', ...
            d, pid, tnum, tvec(s), rvec(s), fvec(s));
    end
end
fclose(fid);
task_signal_info = dir('ARU_Task_Preprocessed.csv');
fprintf('  -> %s (%.1f MB)\n', task_signal_info.name, task_signal_info.bytes / 1e6);

% --- Rest Signal CSV ---
fprintf('Writing ARU_Rest_Preprocessed.csv ...\n');
fid = fopen('ARU_Rest_Preprocessed.csv', 'w');
fprintf(fid, 'DyadID,ParticipantID,RestNum,Time_s,RESP_ARU,ArtifactFlag\n');

for i = 1:rest_count
    row = rest_signal_rows{i};
    d = row{1}; pid = row{2}; rnum = row{3};
    tvec = row{4}; rvec = row{5}; fvec = row{6};
    for s = 1:length(tvec)
        fprintf(fid, '%d,%d,%d,%.4f,%.4f,%d\n', ...
            d, pid, rnum, tvec(s), rvec(s), fvec(s));
    end
end
fclose(fid);
rest_signal_info = dir('ARU_Rest_Preprocessed.csv');
fprintf('  -> %s (%.1f MB)\n', rest_signal_info.name, rest_signal_info.bytes / 1e6);

%% ========================================================================
%  Count artifacts for summary
%  ========================================================================

total_task_samples = 0; total_task_artifacts = 0;
for i = 1:task_count
    flags = task_signal_rows{i}{6};
    total_task_samples = total_task_samples + length(flags);
    total_task_artifacts = total_task_artifacts + sum(flags);
end

total_rest_samples = 0; total_rest_artifacts = 0;
for i = 1:rest_count
    flags = rest_signal_rows{i}{6};
    total_rest_samples = total_rest_samples + length(flags);
    total_rest_artifacts = total_rest_artifacts + sum(flags);
end

%% ========================================================================
%  Generate JSON Metadata Sidecars (BIDS-inspired)
%  ========================================================================

timestamp_str = datestr(now, 'yyyy-mm-ddTHH:MM:SS');

% ---- Shared preprocessing documentation ----
preproc_doc = struct( ...
    'Description', 'Standard respiration preprocessing pipeline following psychophysiology publication recommendations.', ...
    'Steps', {{'1. Zero-phase low-pass Butterworth filter (filtfilt)', ...
               '2. Signal inversion (-ARU; Brain Products convention -> physio convention)', ...
               '3. Downsampling by decimation (every Nth sample)', ...
               '4. Adaptive artifact detection (amplitude clipping + flatline)'}}, ...
    'LowPassFilter', struct( ...
        'Type', 'Butterworth', ...
        'Order', filter_order, ...
        'CutoffFrequency', lp_cutoff, ...
        'CutoffUnit', 'Hz', ...
        'Implementation', 'Zero-phase (MATLAB filtfilt)', ...
        'HighPassApplied', false, ...
        'Rationale', 'Low-pass at 1.0 Hz removes movement artifacts, cardiac contamination, and high-frequency noise while preserving the full respiratory bandwidth (normal adult range: 12-20 br/min = 0.2-0.33 Hz; physiological maximum ~40 br/min = 0.67 Hz). No high-pass filter is applied: baseline drift and slow trends are intentionally preserved as they carry the P(1) dynamics of theoretical interest (task-dependent respiratory flux). Standard low-pass recommendation per Lorig (2012).'), ...
    'SignalInversion', struct( ...
        'Applied', true, ...
        'Method', 'Multiply by -1', ...
        'Rationale', 'Brain Products respiration belt records lung extension as negative-going deflections. Inversion yields the physiological convention: inspiration peak = positive, expiration trough = negative.'), ...
    'Downsampling', struct( ...
        'OriginalRate', fs_original, ...
        'TargetRate', fs_target, ...
        'Factor', downsample_factor, ...
        'Method', 'Decimation (every Nth sample after anti-alias filtering)', ...
        'Rationale', '25 Hz exceeds Nyquist for 1 Hz low-pass filtered signal. Preserves full signal fidelity with 40x data reduction while retaining sufficient temporal resolution for accurate peak/trough detection.'), ...
    'ArtifactDetection', struct( ...
        'Method', 'Combined adaptive amplitude clipping and flatline detection', ...
        'AmplitudeClipping', struct( ...
            'Threshold', 'median +/- 3.0 * IQR (computed per-recording)', ...
            'Rationale', 'Adaptive threshold accounts for between-participant variability in respiratory amplitude. IQR-based approach is robust to non-Gaussian signal distributions.'), ...
        'FlatlineDetection', struct( ...
            'MinDuration_s', flatline_dur, ...
            'GradientThreshold_ARU_per_s', flatline_grad, ...
            'Rationale', 'Detects periods of belt displacement or loss of contact where signal gradient approaches zero for >= 2 seconds.'), ...
        'Action', 'Flag only (no removal or interpolation)'));

source_data_doc = struct( ...
    'Description', 'Raw respiration recordings from perceptual coupling experiment.', ...
    'OriginalSamplingFrequency', fs_original, ...
    'OriginalSamplingFrequencyUnit', 'Hz', ...
    'FileFormat', '.mat (MATLAB)', ...
    'RecordingDevice', 'Brain Products Respiration Belt', ...
    'SignalUnit', 'ARU (Arbitrary Respiration Units)', ...
    'SignalDescription', 'Strain-gauge belt measuring thoracic circumference changes during breathing. Lung extension produces negative-going deflections in the raw signal.');

references = { ...
    'Lorig, T.S. (2012). Respiratory measurement: Techniques. In J.T. Cacioppo, L.G. Tassinary, & G.G. Berntson (Eds.), Handbook of Psychophysiology (3rd ed., pp. 231-238). Cambridge University Press.', ...
    'Kreibig, S.D. (2010). Autonomic nervous system activity in emotion: A review. Biological Psychology, 84(3), 394-421.', ...
    'Massaroni, C., Nicolò, A., Lo Presti, D., Sacchetti, M., Silvestri, S., & Schena, E. (2019). Contact-based methods for measuring respiratory rate. Sensors, 19(4), 908.'};

generated_by = struct( ...
    'Name', 'preprocessARU.m', ...
    'Description', 'MATLAB preprocessing script for PCE respiration data.', ...
    'GenerationDateTime', timestamp_str);

% ---- Task Signal JSON ----
task_sig_meta = struct();
task_sig_meta.Name = 'Perceptual Crossing Experiment - Respiration Task Condition (Preprocessed Signal)';
task_sig_meta.Description = 'Preprocessed respiration signal (ARU) from task trials of a perceptual coupling experiment. Each trial is 60 seconds. Data are organized in long format with one row per sample per trial. Signal is low-pass filtered, inverted to physiological convention (inspiration = positive), and downsampled.';
task_sig_meta.SamplingFrequency = fs_target;
task_sig_meta.SamplingFrequencyUnit = 'Hz';
task_sig_meta.StartTime = 0;
task_sig_meta.Columns = {'DyadID', 'ParticipantID', 'TrialNum', 'Time_s', 'RESP_ARU', 'ArtifactFlag'};

task_sig_meta.DyadID = struct('LongName', 'Dyad Identifier', ...
    'Description', 'Numeric identifier for the experimental dyad (pair of participants). Extracted from the first two digits of the pce folder name.');
task_sig_meta.ParticipantID = struct('LongName', 'Participant Identifier', ...
    'Description', 'Participant number within the dyad (1 or 2).');
task_sig_meta.TrialNum = struct('LongName', 'Trial Number', ...
    'Description', 'Task trial number within the session (1-18).');
task_sig_meta.Time_s = struct('LongName', 'Time', ...
    'Description', 'Elapsed time from the start of the trial.', 'Units', 's');
task_sig_meta.RESP_ARU = struct('LongName', 'Respiration Signal', ...
    'Description', 'Respiratory signal after low-pass filtering, inversion, and downsampling. Positive deflections correspond to inspiration (lung expansion), negative deflections to expiration. No amplitude normalization applied; units are Arbitrary Respiration Units from the Brain Products belt.', ...
    'Units', 'ARU');
task_sig_meta.ArtifactFlag = struct('LongName', 'Artifact Flag', ...
    'Description', 'Binary flag indicating samples that failed artifact checks (amplitude clipping or flatline). 0 = clean, 1 = flagged. Flagged samples are retained for transparency.', ...
    'Levels', struct('x0', 'Clean sample', 'x1', 'Artifact detected'));

task_sig_meta.SourceData = source_data_doc;
task_sig_meta.SourceData.TrialDuration = task_duration;
task_sig_meta.SourceData.TrialDurationUnit = 's';
task_sig_meta.SourceData.TrialsPerParticipant = num_trials;
task_sig_meta.SourceData.ParticipantsPerDyad = 2;
task_sig_meta.Preprocessing = preproc_doc;
task_sig_meta.References = references;
task_sig_meta.DataSummary = struct( ...
    'NumDyads', length(dyad_ids_found), ...
    'NumTrialsLoaded', task_count, ...
    'NumTrialsSkipped', task_skipped, ...
    'SamplesPerTrial', task_duration * fs_target, ...
    'TotalSamples', total_task_samples, ...
    'TotalArtifactsFlagged', total_task_artifacts, ...
    'ArtifactRate', total_task_artifacts / max(total_task_samples, 1), ...
    'MeanAbsAmplitude_ARU', mean(task_amp_means));
task_sig_meta.GeneratedBy = generated_by;

write_json('ARU_Task_Preprocessed.json', task_sig_meta);

% ---- Rest Signal JSON ----
rest_sig_meta = struct();
rest_sig_meta.Name = 'Perceptual Crossing Experiment - Respiration Rest Condition (Preprocessed Signal)';
rest_sig_meta.Description = 'Preprocessed respiration signal (ARU) from rest periods of a perceptual coupling experiment. Each rest period is 180 seconds. Data are organized in long format with one row per sample per period.';
rest_sig_meta.SamplingFrequency = fs_target;
rest_sig_meta.SamplingFrequencyUnit = 'Hz';
rest_sig_meta.StartTime = 0;
rest_sig_meta.Columns = {'DyadID', 'ParticipantID', 'RestNum', 'Time_s', 'RESP_ARU', 'ArtifactFlag'};

rest_sig_meta.DyadID = struct('LongName', 'Dyad Identifier', ...
    'Description', 'Numeric identifier for the experimental dyad (pair of participants). Extracted from the first two digits of the pce folder name.');
rest_sig_meta.ParticipantID = struct('LongName', 'Participant Identifier', ...
    'Description', 'Participant number within the dyad (1 or 2).');
rest_sig_meta.RestNum = struct('LongName', 'Rest Period Number', ...
    'Description', 'Rest period number within the session (1-4).');
rest_sig_meta.Time_s = struct('LongName', 'Time', ...
    'Description', 'Elapsed time from the start of the rest period.', 'Units', 's');
rest_sig_meta.RESP_ARU = struct('LongName', 'Respiration Signal', ...
    'Description', 'Respiratory signal after low-pass filtering, inversion, and downsampling. Positive = inspiration, negative = expiration.', ...
    'Units', 'ARU');
rest_sig_meta.ArtifactFlag = struct('LongName', 'Artifact Flag', ...
    'Description', 'Binary flag indicating samples that failed artifact checks. 0 = clean, 1 = flagged.', ...
    'Levels', struct('x0', 'Clean sample', 'x1', 'Artifact detected'));

rest_sig_meta.SourceData = source_data_doc;
rest_sig_meta.SourceData.RestDuration = rest_duration;
rest_sig_meta.SourceData.RestDurationUnit = 's';
rest_sig_meta.SourceData.RestPeriodsPerParticipant = num_rests;
rest_sig_meta.SourceData.ParticipantsPerDyad = 2;
rest_sig_meta.Preprocessing = preproc_doc;
rest_sig_meta.References = references;
rest_sig_meta.DataSummary = struct( ...
    'NumDyads', length(dyad_ids_found), ...
    'NumPeriodsLoaded', rest_count, ...
    'NumPeriodsSkipped', rest_skipped, ...
    'SamplesPerPeriod', rest_duration * fs_target, ...
    'TotalSamples', total_rest_samples, ...
    'TotalArtifactsFlagged', total_rest_artifacts, ...
    'ArtifactRate', total_rest_artifacts / max(total_rest_samples, 1), ...
    'MeanAbsAmplitude_ARU', mean(rest_amp_means));
rest_sig_meta.GeneratedBy = generated_by;

write_json('ARU_Rest_Preprocessed.json', rest_sig_meta);

%% ========================================================================
%  Final Summary
%  ========================================================================

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
total_csv_bytes = task_signal_info.bytes + rest_signal_info.bytes;
fprintf('  Output CSV total size:    %.1f MB\n', total_csv_bytes / 1e6);

fprintf('\n  Output files:\n');
fprintf('    ARU_Task_Preprocessed.csv   (continuous signal, task)\n');
fprintf('    ARU_Task_Preprocessed.json  (metadata sidecar)\n');
fprintf('    ARU_Rest_Preprocessed.csv   (continuous signal, rest)\n');
fprintf('    ARU_Rest_Preprocessed.json  (metadata sidecar)\n');
fprintf('==========================================================\n');
fprintf('Done.\n');


%% ========================================================================
%  LOCAL FUNCTION: write_json
%  ========================================================================
%  Writes a MATLAB struct as a pretty-printed JSON file.
%  ========================================================================

function write_json(filename, meta_struct)
    json_str = jsonencode(meta_struct);
    json_str = prettify_json(json_str);
    fid = fopen(filename, 'w');
    fprintf(fid, '%s', json_str);
    fclose(fid);
    fprintf('Writing %s ... done.\n', filename);
end


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