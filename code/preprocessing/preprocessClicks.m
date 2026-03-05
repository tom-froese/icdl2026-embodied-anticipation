%% extractClicks.m
% =========================================================================
% Click Response-Time Preprocessing and CSV/JSON Export Script
% =========================================================================
%
% PURPOSE:
%   Reads raw trial data from the perceptual crossing experiment (PCE)
%   folder structure, extracts the timestamp of the first button press
%   for each participant in each trial, and exports consolidated CSV and
%   BIDS-inspired JSON metadata sidecar files.
%
% EXTRACTION LOGIC:
%   For each trial CSV file, the script identifies the first row where
%   button0 == 1 (participant 1's response) and/or button1 == 1
%   (participant 2's response). The corresponding timestamp is recorded
%   as the participant's click response time for that trial. If no click
%   is detected within the trial, the entry is marked as NaN (missed).
%
% INPUT:
%   Run this script from the root data folder containing pce* subfolders.
%   Each subfolder contains a 'trials/' directory with CSV files named:
%     pair_XX_trial_Y.csv
%   where XX is the pair number (zero-padded) and Y is the trial number.
%
%   CSV columns expected:
%     index, timestamp, static_object_0, static_object_1,
%     motor_0_vibrate_software, motor_1_vibrate_software,
%     pos0, button0, shadow_delta0, pos1, button1, shadow_delta1
%
% OUTPUT:
%   ClickResponseTimes.csv   - One row per participant per trial (long fmt)
%   ClickResponseTimes.json  - BIDS-style metadata sidecar
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   February 2026
% =========================================================================

%% Parameters
trial_duration   = 60;     % Trial duration (s)
num_trials_exp   = 18;     % Expected number of trials per experiment
num_participants = 2;      % Participants per dyad (button0 = P1, button1 = P2)

fprintf('==========================================================\n');
fprintf('  Click Response-Time Preprocessing Script (PCE)\n');
fprintf('==========================================================\n');
fprintf('  Trial duration:           %d s\n', trial_duration);
fprintf('  Expected trials/dyad:     %d\n', num_trials_exp);
fprintf('  Participants per dyad:    %d\n', num_participants);
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

%% Preallocate output storage
% Each row: {DyadID, ParticipantID, TrialNum, ClickTime_s, Clicked}
click_rows = {};
row_count  = 0;

% Counters for summary
total_trials_found    = 0;
total_clicks_detected = 0;
total_clicks_missed   = 0;
dyad_ids_found        = [];
click_times_all       = [];   % For summary statistics

%% Main extraction loop
for f = 1:length(folders)
    folder_name = folders(f).name;
    trials_path = fullfile(mainDir, folder_name, 'trials');

    % Extract dyad number from folder name (digits after 'pce')
    dyad_num_str = regexp(folder_name, 'pce(\d+)', 'tokens');
    if isempty(dyad_num_str)
        warning('  Could not parse dyad ID from folder: %s. Skipping.', folder_name);
        continue;
    end
    dyad_id = str2double(dyad_num_str{1}{1});

    if ~ismember(dyad_id, dyad_ids_found)
        dyad_ids_found = [dyad_ids_found, dyad_id]; %#ok<AGROW>
    end

    fprintf('Processing folder: %s (Dyad %02d)\n', folder_name, dyad_id);

    % Check that trials subfolder exists
    if ~isfolder(trials_path)
        warning('  No trials/ subfolder found in %s. Skipping.', folder_name);
        continue;
    end

    % Get all CSV files in the trials folder
    csv_files = dir(fullfile(trials_path, '*.csv'));

    if isempty(csv_files)
        warning('  No CSV files found in %s/trials/. Skipping.', folder_name);
        continue;
    end

    for j = 1:length(csv_files)
        csv_name = csv_files(j).name;
        csv_path = fullfile(trials_path, csv_name);

        % Parse trial number from filename (e.g., pair_01_trial_3.csv)
        trial_tok = regexp(csv_name, 'trial_?(\d+)', 'tokens');
        if isempty(trial_tok)
            warning('  Could not parse trial number from: %s. Skipping.', csv_name);
            continue;
        end
        trial_num = str2double(trial_tok{1}{1});
        total_trials_found = total_trials_found + 1;

        % Read the CSV file
        try
            data = readtable(csv_path);
        catch ME
            warning('  Error reading %s: %s. Skipping.', csv_path, ME.message);
            continue;
        end

        % Validate required columns exist
        required_cols = {'timestamp', 'button0', 'button1'};
        missing_cols = setdiff(required_cols, data.Properties.VariableNames);
        if ~isempty(missing_cols)
            warning('  %s missing columns: %s. Skipping.', ...
                csv_name, strjoin(missing_cols, ', '));
            continue;
        end

        % --- Extract first click for each participant ---
        % Participant 1 (button0)
        b0_idx = find(data.button0 == 1, 1, 'first');
        if ~isempty(b0_idx)
            click_time_p1 = data.timestamp(b0_idx);
            clicked_p1 = 1;
            total_clicks_detected = total_clicks_detected + 1;
            click_times_all = [click_times_all, click_time_p1]; %#ok<AGROW>
        else
            click_time_p1 = NaN;
            clicked_p1 = 0;
            total_clicks_missed = total_clicks_missed + 1;
        end

        row_count = row_count + 1;
        click_rows{row_count} = {dyad_id, 1, trial_num, click_time_p1, clicked_p1}; %#ok<SAGROW>

        % Participant 2 (button1)
        b1_idx = find(data.button1 == 1, 1, 'first');
        if ~isempty(b1_idx)
            click_time_p2 = data.timestamp(b1_idx);
            clicked_p2 = 1;
            total_clicks_detected = total_clicks_detected + 1;
            click_times_all = [click_times_all, click_time_p2]; %#ok<AGROW>
        else
            click_time_p2 = NaN;
            clicked_p2 = 0;
            total_clicks_missed = total_clicks_missed + 1;
        end

        row_count = row_count + 1;
        click_rows{row_count} = {dyad_id, 2, trial_num, click_time_p2, clicked_p2}; %#ok<SAGROW>
    end
end

fprintf('\n==========================================================\n');
fprintf('  Data extraction complete.\n');
fprintf('  Trial files processed:  %d\n', total_trials_found);
fprintf('  Clicks detected:        %d\n', total_clicks_detected);
fprintf('  Clicks missed (NaN):    %d\n', total_clicks_missed);
fprintf('  Total rows:             %d\n', row_count);
fprintf('==========================================================\n\n');

%% Sort rows by DyadID, ParticipantID, TrialNum
sort_keys = zeros(row_count, 3);
for i = 1:row_count
    sort_keys(i, :) = [click_rows{i}{1}, click_rows{i}{2}, click_rows{i}{3}];
end
[~, sort_idx] = sortrows(sort_keys, [1 2 3]);
click_rows = click_rows(sort_idx);

%% Export CSV
output_csv = 'ClickResponseTimes.csv';
fprintf('Writing %s ...\n', output_csv);

fid = fopen(output_csv, 'w');
fprintf(fid, 'DyadID,ParticipantID,TrialNum,ClickTime_s,Clicked\n');

for i = 1:row_count
    row = click_rows{i};
    dyad_id   = row{1};
    part_id   = row{2};
    trial_num = row{3};
    click_t   = row{4};
    clicked   = row{5};

    if isnan(click_t)
        fprintf(fid, '%d,%d,%d,NaN,%d\n', dyad_id, part_id, trial_num, clicked);
    else
        fprintf(fid, '%d,%d,%d,%.6f,%d\n', dyad_id, part_id, trial_num, click_t, clicked);
    end
end

fclose(fid);

csv_info = dir(output_csv);
fprintf('  -> %s (%.1f KB)\n', csv_info.name, csv_info.bytes / 1e3);

%% Compute summary statistics
if ~isempty(click_times_all)
    mean_click  = mean(click_times_all);
    std_click   = std(click_times_all);
    median_click = median(click_times_all);
    min_click   = min(click_times_all);
    max_click   = max(click_times_all);
    iqr_click   = iqr(click_times_all);
else
    mean_click   = NaN;
    std_click    = NaN;
    median_click = NaN;
    min_click    = NaN;
    max_click    = NaN;
    iqr_click    = NaN;
end

total_possible = total_clicks_detected + total_clicks_missed;
click_rate = total_clicks_detected / max(total_possible, 1);

%% ========================================================================
%  Generate JSON Metadata Sidecar (BIDS-inspired)
%  ========================================================================

timestamp_str = datestr(now, 'yyyy-mm-ddTHH:MM:SS');

meta = struct();

% Dataset-level fields
meta.Name = 'Perceptual Crossing Experiment - Click Response Times';
meta.Description = ['First button-press response times extracted from ' ...
    'perceptual crossing experiment trials. Each row represents one ' ...
    'participant''s response in one trial. Trials last ' ...
    num2str(trial_duration) ' seconds. If no click was recorded within ' ...
    'the trial, ClickTime_s is NaN and Clicked is 0.'];

% Column specification (BIDS-style)
meta.Columns = {'DyadID', 'ParticipantID', 'TrialNum', 'ClickTime_s', 'Clicked'};

meta.DyadID = struct( ...
    'LongName', 'Dyad Identifier', ...
    'Description', ['Numeric identifier for the experimental dyad. ' ...
        'Extracted from the pce folder name (e.g., pce01230807 -> Dyad 1).']);

meta.ParticipantID = struct( ...
    'LongName', 'Participant Identifier', ...
    'Description', ['Participant number within the dyad. ' ...
        '1 = button0 column in raw data; 2 = button1 column.']);

meta.TrialNum = struct( ...
    'LongName', 'Trial Number', ...
    'Description', ['Trial number within the session (1-' ...
        num2str(num_trials_exp) '). Parsed from the CSV filename.']);

meta.ClickTime_s = struct( ...
    'LongName', 'Click Response Time', ...
    'Description', ['Timestamp (in seconds from trial onset) of the ' ...
        'participant''s first button press. NaN if no click was ' ...
        'recorded during the trial.'], ...
    'Units', 's');

meta.Clicked = struct( ...
    'LongName', 'Click Detected', ...
    'Description', ['Binary flag indicating whether a button press ' ...
        'was detected within the trial.'], ...
    'Levels', struct('x0', 'No click detected (missed trial)', ...
                     'x1', 'Click detected'));

% Source data information
meta.SourceData = struct( ...
    'Description', 'Raw CSV trial recordings from perceptual crossing experiment.', ...
    'FileFormat', '.csv', ...
    'TrialDuration', trial_duration, ...
    'TrialDurationUnit', 's', ...
    'TrialsPerDyad', num_trials_exp, ...
    'ParticipantsPerDyad', num_participants, ...
    'FolderStructure', 'pceXXYYMMDD/trials/pair_XX_trial_Y.csv', ...
    'RawColumns', {{'index', 'timestamp', 'static_object_0', ...
        'static_object_1', 'motor_0_vibrate_software', ...
        'motor_1_vibrate_software', 'pos0', 'button0', ...
        'shadow_delta0', 'pos1', 'button1', 'shadow_delta1'}}, ...
    'ButtonMapping', struct( ...
        'button0', 'Participant 1 response', ...
        'button1', 'Participant 2 response'));

% Extraction method
meta.ExtractionMethod = struct( ...
    'Description', ['For each trial CSV, the script finds the first ' ...
        'row where buttonX == 1 and records the corresponding ' ...
        'timestamp value. If no row satisfies this condition, the ' ...
        'click time is set to NaN.'], ...
    'Steps', {{'1. Read trial CSV with readtable()', ...
               '2. Locate first occurrence of button0 == 1 (P1)', ...
               '3. Locate first occurrence of button1 == 1 (P2)', ...
               '4. Record timestamp at first-click index, or NaN if absent'}});

% Data summary statistics
meta.DataSummary = struct( ...
    'NumDyads', length(dyad_ids_found), ...
    'NumTrialFilesProcessed', total_trials_found, ...
    'TotalClicksDetected', total_clicks_detected, ...
    'TotalClicksMissed', total_clicks_missed, ...
    'ClickRate', click_rate, ...
    'MeanClickTime_s', mean_click, ...
    'StdClickTime_s', std_click, ...
    'MedianClickTime_s', median_click, ...
    'MinClickTime_s', min_click, ...
    'MaxClickTime_s', max_click, ...
    'IQRClickTime_s', iqr_click);

% File provenance
meta.GeneratedBy = struct( ...
    'Name', 'preprocessClicks.m', ...
    'Description', 'MATLAB preprocessing script for PCE click response-time extraction.', ...
    'GenerationDateTime', timestamp_str);

% Write JSON
output_json = 'ClickResponseTimes.json';
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
fprintf('  Clicks detected:          %d / %d (%.1f%%)\n', ...
    total_clicks_detected, total_possible, 100 * click_rate);
fprintf('  Clicks missed:            %d\n', total_clicks_missed);
fprintf('\n');
fprintf('  Click time statistics (detected clicks only):\n');
fprintf('    Mean:    %.2f s\n', mean_click);
fprintf('    Std:     %.2f s\n', std_click);
fprintf('    Median:  %.2f s\n', median_click);
fprintf('    Min:     %.2f s\n', min_click);
fprintf('    Max:     %.2f s\n', max_click);
fprintf('    IQR:     %.2f s\n', iqr_click);
fprintf('\n');
fprintf('  Output files:\n');
fprintf('    %s   (data)\n', output_csv);
fprintf('    %s  (metadata sidecar)\n', output_json);
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
