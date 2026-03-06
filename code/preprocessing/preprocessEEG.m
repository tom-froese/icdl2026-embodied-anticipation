%% preprocessEEG.m
% =========================================================================
% Beta-Band EEG Coherence Preprocessing and CSV/JSON Export Script
% =========================================================================
%
% PURPOSE:
%   Computes magnitude-squared coherence (MSC) in the beta band (13-30 Hz)
%   for two fixed 4-second epochs per trial, chosen relative to the
%   PAS-crossover confidence interval:
%
%     PRE  epoch : [17 s, 21 s]  -- sampling phase (before CI lower bound)
%     POST epoch : [34 s, 38 s]  -- selection phase (after CI upper bound)
%
%   Crossover context (from plotPASFigures.m):
%     PAS 4 < PAS 3 crossover: 27.73 s  (95% CI: 21.49-33.39 s)
%     Click-time peak:         ~24.4 s   (falls within CI)
%
%   By comparing PRE vs POST we bracket the perceptual transition while
%   avoiding the ambiguous CI window entirely.
%
% WELCH PARAMETERS (per epoch):
%   4000-sample epoch, 1 s segments (1000 pts), 50% overlap -> 7 averages,
%   Hann window, 1 Hz frequency resolution.
%
% INPUT:
%   Run from the directory containing preprocessed .mat files named:
%     pceXX_P{1,2}_Trial{1-18}.mat   (64 x 60000, float32, in Volts)
%   These files are produced by MNE-Python preprocessing (see Methods).
%
% OUTPUT:
%   BetaCoherence.csv      - Full pair-level: DyadID, ParticipantID,
%                            TrialNum, PairIdx, ChA, ChB, Coh_PRE, Coh_POST
%   BetaCoherence_ROI.csv  - ROI-level: 6 ROIs x 2 epochs per trial
%   BetaCoherence.json     - BIDS-style metadata sidecar
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   March 2026
% =========================================================================

clear; clc;

%% ========================================================================
%  1. PARAMETERS
%  ========================================================================

fs         = 1000;      % Sampling rate (Hz)
trial_dur  = 60;        % Trial duration (s)
n_samp     = fs * trial_dur;

% Epoch definitions (seconds)
epoch_pre  = [17 21];   % Sampling phase: 4 s before CI lower bound (21.49 s)
epoch_post = [34 38];   % Selection phase: 4 s after CI upper bound (33.39 s)

pre_idx    = (epoch_pre(1) * fs + 1) : (epoch_pre(2)  * fs);   % 17001-21000
post_idx   = (epoch_post(1)* fs + 1) : (epoch_post(2) * fs);   % 34001-38000
epoch_n    = length(pre_idx);   % 4000 samples (= 4 s)

% Welch spectral estimation settings
welch_n    = 1 * fs;            % 1 s segment length
welch_ov   = floor(welch_n / 2); % 50% overlap

% Beta band
beta_lo    = 13;        % Hz
beta_hi    = 30;        % Hz

fprintf('==========================================================\n');
fprintf('  Beta-Band Coherence Preprocessing Script (PCE)\n');
fprintf('==========================================================\n');
fprintf('  Sampling rate:        %d Hz\n', fs);
fprintf('  Trial duration:       %d s\n', trial_dur);
fprintf('  Beta band:            %d-%d Hz\n', beta_lo, beta_hi);
fprintf('  PRE  epoch:           [%d, %d] s  (sampling phase)\n', ...
    epoch_pre(1), epoch_pre(2));
fprintf('  POST epoch:           [%d, %d] s  (selection phase)\n', ...
    epoch_post(1), epoch_post(2));
fprintf('  Epoch length:         %d s (%d samples)\n', epoch_n/fs, epoch_n);
fprintf('  Welch segment:        %d ms, %d%% overlap\n', ...
    welch_n/fs*1000, round(welch_ov/welch_n*100));
fprintf('==========================================================\n\n');

%% ========================================================================
%  2. CHANNEL MAPPING -- Custom 64Ch actiCAP snap (CACS-64-X22)
%  ========================================================================
%  27 channels selected for coherence analysis (see Putri & Froese, 2025).
%  ch_hw maps from analysis-order index to hardware channel index (1-based)
%  in the 64-channel .mat files.

ch_labels = { ...
    'Fp1','Fp2','F7','F3','Fz','F4','F8', ...
    'FC5','FC1','FC2','FC6', ...
    'T7','C3','Cz','C4','T8', ...
    'CP5','CP1','CP2','CP6', ...
    'P7','P3','Pz','P4','P8', ...
    'O1','O2'};

ch_hw = [ ...
     1, 32,  4,  3,  2, 30, 31, ...   % Fp1 Fp2 F7 F3 Fz F4 F8
     6,  7, 29, 28, ...               % FC5 FC1 FC2 FC6
     9,  8, 24, 25, 26, ...           % T7  C3  Cz  C4  T8
    11, 12, 23, 22, ...               % CP5 CP1 CP2 CP6
    15, 14, 13, 19, 20, ...           % P7  P3  Pz  P4  P8
    16, 18];                           % O1  O2

n_ch    = length(ch_labels);
pairs   = nchoosek(1:n_ch, 2);
n_pairs = size(pairs, 1);  % 351

fprintf('  Channels: %d -> %d unique pairs\n\n', n_ch, n_pairs);

%% ========================================================================
%  3. ROI DEFINITIONS
%  ========================================================================

frontal_idx  = find(ismember(ch_labels, ...
    {'Fp1','Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6'}));
central_idx  = find(ismember(ch_labels, ...
    {'T7','C3','Cz','C4','T8'}));
parietal_idx = find(ismember(ch_labels, ...
    {'CP5','CP1','CP2','CP6','P7','P3','Pz','P4','P8','O1','O2'}));
left_idx     = find(ismember(ch_labels, ...
    {'Fp1','F7','F3','FC5','FC1','T7','C3','CP5','CP1','P7','P3','O1'}));
right_idx    = find(ismember(ch_labels, ...
    {'Fp2','F8','F4','FC6','FC2','T8','C4','CP6','CP2','P8','P4','O2'}));

roi_names = {'Frontal','Central','Parietal','Left','Right','FrontoParietal'};
roi_sets  = {frontal_idx, central_idx, parietal_idx, left_idx, right_idx, []};

% Precompute which channel pairs belong to each ROI
roi_pidx = cell(1, 6);
for r = 1:5
    m = roi_sets{r};
    mask = false(n_pairs, 1);
    for p = 1:n_pairs
        mask(p) = all(ismember(pairs(p,:), m));
    end
    roi_pidx{r} = find(mask);
end

% FrontoParietal: cross-region pairs (one Frontal + one Parietal electrode)
fp = false(n_pairs, 1);
for p = 1:n_pairs
    a = pairs(p,1);  b = pairs(p,2);
    fp(p) = (ismember(a, frontal_idx) && ismember(b, parietal_idx)) || ...
            (ismember(a, parietal_idx) && ismember(b, frontal_idx));
end
roi_pidx{6} = find(fp);

for r = 1:length(roi_names)
    fprintf('  ROI %-16s  %3d pairs\n', roi_names{r}, length(roi_pidx{r}));
end
fprintf('\n');

%% ========================================================================
%  4. DISCOVER FILES
%  ========================================================================

mainDir   = pwd;
mat_files = dir(fullfile(mainDir, 'pce*', 'pce*_P*_Trial*.mat'));

if isempty(mat_files)
    error(['No pce*_P*_Trial*.mat files found in:\n  %s\n' ...
           'Please run this script from the preprocessed EEG directory.'], mainDir);
end

% Parse filenames into structured array
fi = struct('path', {}, 'dyad', {}, 'pid', {}, 'trial', {});
for f = 1:length(mat_files)
    tok = regexp(mat_files(f).name, ...
        'pce(\d+)_P(\d+)_Trial(\d+)\.mat', 'tokens');
    if isempty(tok)
        warning('  Could not parse: %s. Skipping.', mat_files(f).name);
        continue;
    end
    fi(end+1).path  = fullfile(mat_files(f).folder, mat_files(f).name); %#ok<SAGROW>
    fi(end).dyad    = str2double(tok{1}{1});
    fi(end).pid     = str2double(tok{1}{2});
    fi(end).trial   = str2double(tok{1}{3});
end

% Sort numerically by dyad -> participant -> trial
sort_key = [fi.dyad]' * 10000 + [fi.pid]' * 100 + [fi.trial]';
[~, ord] = sort(sort_key);
fi = fi(ord);

n_files    = length(fi);
dyad_ids_found = unique([fi.dyad]);

fprintf('  Found %d files (%d dyads).\n\n', n_files, length(dyad_ids_found));

%% ========================================================================
%  5. PRECOMPUTE FREQUENCY INDICES
%  ========================================================================

[~, fvec] = mscohere(zeros(epoch_n, 1), zeros(epoch_n, 1), ...
    hann(welch_n), welch_ov, welch_n, fs);
beta_idx = find(fvec >= beta_lo & fvec <= beta_hi);

fprintf('  Frequency resolution: %.1f Hz\n', fvec(2) - fvec(1));
fprintf('  Beta bins: %d (%.0f-%.0f Hz)\n\n', ...
    length(beta_idx), fvec(beta_idx(1)), fvec(beta_idx(end)));

%% ========================================================================
%  6. OPEN OUTPUT FILES
%  ========================================================================

output_csv_full = 'BetaCoherence.csv';
output_csv_roi  = 'BetaCoherence_ROI.csv';

fprintf('Writing %s ...\n', output_csv_full);
fid_full = fopen(output_csv_full, 'w');
fprintf(fid_full, 'DyadID,ParticipantID,TrialNum,PairIdx,ChA,ChB,Coh_PRE,Coh_POST\n');

fprintf('Writing %s ...\n', output_csv_roi);
fid_roi = fopen(output_csv_roi, 'w');
fprintf(fid_roi, 'DyadID,ParticipantID,TrialNum,ROI,Coh_PRE,Coh_POST\n');

%% ========================================================================
%  7. MAIN LOOP
%  ========================================================================

n_proc = 0;
n_skip = 0;
expected_bytes = 64 * n_samp * 4;   % float32: 64 channels x 60000 samples x 4 bytes

tic;

for f = 1:n_files
    d = fi(f);

    fprintf('  [%4d/%4d] pce%02d P%d Trial%02d  loading...', ...
        f, n_files, d.dyad, d.pid, d.trial);

    % --- File size pre-check (reject files < 90% of expected size) ---
    finfo = dir(d.path);
    if isempty(finfo) || finfo.bytes < expected_bytes * 0.9
        fprintf(' SKIP (file missing or too small: %d bytes)\n', ...
            finfo.bytes);
        n_skip = n_skip + 1;
        continue;
    end

    % --- Load ---
    try
        S = load(d.path);
        fn = fieldnames(S);
        raw = S.(fn{1});
    catch ME
        fprintf(' SKIP (load error: %s)\n', ME.message);
        n_skip = n_skip + 1;
        continue;
    end
    fprintf(' loaded.');

    % --- Validate dimensions ---
    if size(raw, 1) ~= 64 || size(raw, 2) ~= n_samp
        fprintf(' SKIP (unexpected dims %dx%d)\n', size(raw, 1), size(raw, 2));
        n_skip = n_skip + 1;
        continue;
    end

    % --- Extract 27 channels, slice epochs ---
    eeg      = raw(ch_hw, :);          % 27 x 60000
    seg_pre  = eeg(:, pre_idx);        % 27 x 4000
    seg_post = eeg(:, post_idx);       % 27 x 4000
    clear raw eeg;                     % Free memory early

    fprintf(' coherence...');

    % --- Compute MSC for all 351 pairs ---
    coh = zeros(n_pairs, 2);

    for p = 1:n_pairs
        a = pairs(p, 1);  b = pairs(p, 2);

        Cpre  = mscohere(seg_pre(a,:)',  seg_pre(b,:)',  ...
            hann(welch_n), welch_ov, welch_n, fs);
        Cpost = mscohere(seg_post(a,:)', seg_post(b,:)', ...
            hann(welch_n), welch_ov, welch_n, fs);

        coh(p, 1) = mean(Cpre(beta_idx));
        coh(p, 2) = mean(Cpost(beta_idx));
    end

    % --- Write full pair-level CSV ---
    for p = 1:n_pairs
        fprintf(fid_full, '%d,%d,%d,%d,%s,%s,%.6f,%.6f\n', ...
            d.dyad, d.pid, d.trial, p, ...
            ch_labels{pairs(p, 1)}, ch_labels{pairs(p, 2)}, ...
            coh(p, 1), coh(p, 2));
    end

    % --- Write ROI summary CSV ---
    for r = 1:length(roi_names)
        pidx = roi_pidx{r};
        if isempty(pidx), continue; end
        roi_coh = mean(coh(pidx, :), 1);
        fprintf(fid_roi, '%d,%d,%d,%s,%.6f,%.6f\n', ...
            d.dyad, d.pid, d.trial, roi_names{r}, ...
            roi_coh(1), roi_coh(2));
    end

    n_proc = n_proc + 1;
    fprintf(' done (%.0f s elapsed)\n', toc);
end

fclose(fid_full);
fclose(fid_roi);
elapsed = toc;

csv_info_full = dir(output_csv_full);
csv_info_roi  = dir(output_csv_roi);
fprintf('\n  -> %s (%.1f MB)\n', output_csv_full, csv_info_full.bytes / 1e6);
fprintf('  -> %s (%.1f KB)\n', output_csv_roi, csv_info_roi.bytes / 1e3);

%% ========================================================================
%  8. JSON METADATA SIDECAR (BIDS-inspired)
%  ========================================================================

timestamp_str = datestr(now, 'yyyy-mm-ddTHH:MM:SS');

meta = struct();

% Dataset-level fields
meta.Name = 'Perceptual Crossing Experiment - Beta-Band EEG Coherence';
meta.Description = ['Magnitude-squared coherence (MSC) in the standard beta band ' ...
    '(13-30 Hz) computed for two fixed 4-second epochs that bracket the PAS ' ...
    'crossover confidence interval: PRE [17-21 s] (sampling phase) and POST ' ...
    '[34-38 s] (selection phase). The PAS 4/3 crossover occurs at 27.73 s ' ...
    '(95% CI: 21.49-33.39 s); the click-time peak is ~24.4 s. The entire ' ...
    'CI window is excluded from both epochs.'];

% Column specification
meta.Columns_FullPairs = {'DyadID', 'ParticipantID', 'TrialNum', ...
    'PairIdx', 'ChA', 'ChB', 'Coh_PRE', 'Coh_POST'};
meta.Columns_ROI = {'DyadID', 'ParticipantID', 'TrialNum', ...
    'ROI', 'Coh_PRE', 'Coh_POST'};

meta.DyadID = struct( ...
    'LongName', 'Dyad Identifier', ...
    'Description', ['Numeric identifier for the experimental dyad. ' ...
        'Parsed from the .mat filename (e.g., pce01_P1_Trial01.mat -> Dyad 1).']);

meta.ParticipantID = struct( ...
    'LongName', 'Participant Identifier', ...
    'Description', ['Participant number within the dyad (1 or 2). ' ...
        'Parsed from the .mat filename.']);

meta.TrialNum = struct( ...
    'LongName', 'Trial Number', ...
    'Description', 'Trial number within the session (1-18).');

meta.PairIdx = struct( ...
    'LongName', 'Channel Pair Index', ...
    'Description', ['Index into the 351 unique channel pairs (nchoosek(27,2)). ' ...
        'ChA and ChB identify the two electrodes.']);

meta.Coh_PRE = struct( ...
    'LongName', 'Beta Coherence (PRE epoch)', ...
    'Description', ['Mean magnitude-squared coherence in the beta band ' ...
        '(13-30 Hz) during the PRE epoch [17-21 s].'], ...
    'Units', 'dimensionless (0-1)');

meta.Coh_POST = struct( ...
    'LongName', 'Beta Coherence (POST epoch)', ...
    'Description', ['Mean magnitude-squared coherence in the beta band ' ...
        '(13-30 Hz) during the POST epoch [34-38 s].'], ...
    'Units', 'dimensionless (0-1)');

% Epoch definitions
meta.Epochs = struct( ...
    'PRE',  struct('Start_s', epoch_pre(1),  'End_s', epoch_pre(2), ...
                   'Label', 'Sampling phase', ...
                   'Rationale', '4 s immediately before CI lower bound (21.49 s)'), ...
    'POST', struct('Start_s', epoch_post(1), 'End_s', epoch_post(2), ...
                   'Label', 'Selection phase', ...
                   'Rationale', '4 s immediately after CI upper bound (33.39 s)'));

meta.CrossoverContext = struct( ...
    'CrossoverTime_s',    27.73, ...
    'CI_lower_s',         21.49, ...
    'CI_upper_s',         33.39, ...
    'ClickTimePeak_s',    24.4, ...
    'CIWindow_excluded',  true, ...
    'Source', 'plotPASFigures.m');

% Spectral parameters
meta.FrequencyBand = struct( ...
    'Name', 'beta', ...
    'Range_Hz', [beta_lo beta_hi]);

meta.SamplingRate_Hz  = fs;
meta.TrialDuration_s  = trial_dur;
meta.EpochLength_s    = epoch_n / fs;

meta.WelchParameters = struct( ...
    'SegmentLength_s', welch_n / fs, ...
    'OverlapFraction', welch_ov / welch_n, ...
    'NumAveragingSegments', 7, ...
    'FrequencyResolution_Hz', fs / welch_n, ...
    'WindowFunction', 'Hann', ...
    'Function', 'mscohere (MATLAB)');

% Channel information
meta.Channels = struct( ...
    'NumSelected', n_ch, ...
    'Labels', {ch_labels}, ...
    'HardwareIndices_1based', ch_hw, ...
    'NumPairs', n_pairs, ...
    'Cap', 'Custom 64Ch actiCAP snap (CACS-64-X22)');

% ROI definitions
meta.ROIs = struct( ...
    'Names', {roi_names}, ...
    'Frontal', 'Fp1, Fp2, F7, F3, Fz, F4, F8, FC5, FC1, FC2, FC6', ...
    'Central', 'T7, C3, Cz, C4, T8', ...
    'Parietal', 'CP5, CP1, CP2, CP6, P7, P3, Pz, P4, P8, O1, O2', ...
    'Left', 'Fp1, F7, F3, FC5, FC1, T7, C3, CP5, CP1, P7, P3, O1', ...
    'Right', 'Fp2, F8, F4, FC6, FC2, T8, C4, CP6, CP2, P8, P4, O2', ...
    'FrontoParietal', 'Cross-region: one Frontal + one Parietal electrode');

% Source data
meta.SourceData = struct( ...
    'Description', ['Preprocessed EEG .mat files from MNE-Python pipeline. ' ...
        'Data from pair 31 are not included in the repository due to ' ...
        'incomplete recording caused by a system error.'], ...
    'FileFormat', '.mat (float32)', ...
    'Dimensions', '64 channels x 60000 samples', ...
    'FileNaming', 'pceXX_P{1,2}_Trial{1-18}.mat', ...
    'Preprocessing', ['0.1 Hz high-pass, 60 Hz low-pass (zero-phase FIR), ' ...
        'spherical spline interpolation of bad channels, ' ...
        'ICA artifact removal (30 components, EOG detection via Fp1/Fp2).']);

% Data summary
meta.DataSummary = struct( ...
    'NumDyads', length(dyad_ids_found), ...
    'NumFilesProcessed', n_proc, ...
    'NumFilesSkipped', n_skip, ...
    'TotalPairRows', n_proc * n_pairs, ...
    'TotalROIRows', n_proc * length(roi_names), ...
    'ElapsedTime_s', round(elapsed));

% File provenance
meta.GeneratedBy = struct( ...
    'Name', 'preprocessEEG.m', ...
    'Description', 'MATLAB preprocessing script for PCE beta-band EEG coherence.', ...
    'GenerationDateTime', timestamp_str);

% Output file references
meta.OutputFiles = struct( ...
    'FullPairs', output_csv_full, ...
    'ROISummary', output_csv_roi);

%% Write JSON
output_json = 'BetaCoherence.json';
json_str = jsonencode(meta);
json_str = prettify_json(json_str);
fid = fopen(output_json, 'w');
fprintf(fid, '%s', json_str);
fclose(fid);
fprintf('Writing %s ... done.\n', output_json);

%% ========================================================================
%  9. FINAL SUMMARY
%  ========================================================================

fprintf('\n==========================================================\n');
fprintf('  SUMMARY\n');
fprintf('==========================================================\n');
fprintf('  Beta band:            %d-%d Hz\n', beta_lo, beta_hi);
fprintf('  Epochs:               PRE [%d-%d s] | POST [%d-%d s]\n', ...
    epoch_pre(1), epoch_pre(2), epoch_post(1), epoch_post(2));
fprintf('  Dyads:                %d\n', length(dyad_ids_found));
fprintf('  Files processed:      %d\n', n_proc);
fprintf('  Files skipped:        %d\n', n_skip);
fprintf('  Pairs per file:       %d\n', n_pairs);
fprintf('  Elapsed:              %.0f s (%.1f min)\n', elapsed, elapsed/60);
fprintf('\n');
fprintf('  Output files:\n');
fprintf('    %s\n', output_csv_full);
fprintf('    %s\n', output_csv_roi);
fprintf('    %s\n', output_json);
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