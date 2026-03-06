clear; clc;
%% preprocessPAS.m
% =========================================================================
% PAS Click-Presence Rating Preprocessing and CSV/JSON Export Script
% =========================================================================
%
% PURPOSE:
%   Reads PAS (Perceptual Awareness Scale) questionnaire CSV files from
%   the perceptual crossing experiment (PCE) folder structure, extracts
%   the click_presence rating for each participant in each trial, and
%   exports a consolidated CSV file with a BIDS-inspired JSON metadata
%   sidecar.
%
% EXTRACTION LOGIC:
%   Each PAS CSV file (one per participant per session) contains long-
%   format questionnaire responses with columns: trial_id, question_id,
%   answer.  The script filters for question_id == 'click_presence',
%   pivots to wide format, converts 0-indexed trial_id to 1-indexed
%   TrialNum, and assigns DyadID and ParticipantID parsed from the
%   filename (pair_XX_PY_PAS*.csv).
%
% INPUT:
%   Run this script from the root data folder containing pce* subfolders.
%   Each subfolder contains a questionnaires/ directory with CSV files
%   matching *PAS*.csv, named:
%     pair_XX_PY_PAS*.csv
%   where XX is the pair number and Y is the participant number (1 or 2).
%
% OUTPUT:
%   PASRatings.csv   - One row per participant per trial (wide format)
%   PASRatings.json  - BIDS-style metadata sidecar
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   March 2026
% =========================================================================

%% Parameters
num_trials_exp   = 18;     % Expected number of trials per session
num_participants = 2;      % Participants per dyad

fprintf('==========================================================\n');
fprintf('  PAS Click-Presence Preprocessing Script (PCE)\n');
fprintf('==========================================================\n');
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
pas_rows = {};
row_count = 0;

% Counters for summary
total_files_found    = 0;
total_files_skipped  = 0;
total_ratings        = 0;
total_missing        = 0;
dyad_ids_found       = [];
all_ratings          = [];   % For summary statistics

%% Main extraction loop
for f = 1:length(folders)
    folder_name = folders(f).name;
    quest_path  = fullfile(mainDir, folder_name, 'questionnaires');

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

    % Check that questionnaires subfolder exists
    if ~isfolder(quest_path)
        warning('  No questionnaires/ subfolder found in %s. Skipping.', folder_name);
        continue;
    end

    % Get all PAS CSV files in the questionnaires folder
    pas_files = dir(fullfile(quest_path, '*PAS*.csv'));

    if isempty(pas_files)
        warning('  No *PAS*.csv files found in %s/questionnaires/. Skipping.', folder_name);
        continue;
    end

    for j = 1:length(pas_files)
        csv_name = pas_files(j).name;
        csv_path = fullfile(quest_path, csv_name);
        total_files_found = total_files_found + 1;

        % Parse participant ID from filename (pair_XX_PY_...)
        part_tok = regexp(csv_name, 'pair_\d+_P(\d+)', 'tokens');
        if isempty(part_tok)
            warning('  Could not parse participant ID from: %s. Skipping.', csv_name);
            total_files_skipped = total_files_skipped + 1;
            continue;
        end
        participant_id = str2double(part_tok{1}{1});

        % Read the PAS CSV file
        try
            opts = detectImportOptions(csv_path);
            opts.SelectedVariableNames = {'trial_id', 'question_id', 'answer'};
            opts.VariableNamingRule = 'preserve';
            tbl = readtable(csv_path, opts);
        catch ME
            warning('  Error reading %s: %s. Skipping.', csv_name, ME.message);
            total_files_skipped = total_files_skipped + 1;
            continue;
        end

        if isempty(tbl)
            warning('  %s: empty table. Skipping.', csv_name);
            total_files_skipped = total_files_skipped + 1;
            continue;
        end

        % Keep ONLY click_presence responses
        mask = strcmp(tbl.question_id, 'click_presence');
        tbl = tbl(mask, :);

        if isempty(tbl)
            warning('  %s: no click_presence entries found. Skipping.', csv_name);
            total_files_skipped = total_files_skipped + 1;
            continue;
        end

        % --- Extract one rating per trial ---
        for r = 1:height(tbl)
            trial_id_raw = tbl.trial_id(r);
            trial_num    = trial_id_raw + 1;   % Convert 0-indexed to 1-indexed
            rating       = tbl.answer(r);

            row_count = row_count + 1;
            pas_rows{row_count} = {dyad_id, participant_id, trial_num, rating}; %#ok<SAGROW>

            total_ratings = total_ratings + 1;
            all_ratings = [all_ratings, rating]; %#ok<AGROW>
        end
    end
end

% Count missing trials (expected: 2 participants × num_trials_exp per dyad)
total_expected = length(dyad_ids_found) * num_participants * num_trials_exp;
total_missing  = total_expected - total_ratings;

fprintf('\n==========================================================\n');
fprintf('  Data extraction complete.\n');
fprintf('  PAS files processed:    %d\n', total_files_found);
fprintf('  PAS files skipped:      %d\n', total_files_skipped);
fprintf('  Ratings extracted:      %d\n', total_ratings);
fprintf('  Expected total:         %d\n', total_expected);
fprintf('  Missing ratings:        %d\n', total_missing);
fprintf('  Total rows:             %d\n', row_count);
fprintf('==========================================================\n\n');

%% Sort rows by DyadID, ParticipantID, TrialNum
sort_keys = zeros(row_count, 3);
for i = 1:row_count
    sort_keys(i, :) = [pas_rows{i}{1}, pas_rows{i}{2}, pas_rows{i}{3}];
end
[~, sort_idx] = sortrows(sort_keys, [1 2 3]);
pas_rows = pas_rows(sort_idx);

%% Export CSV
output_csv = 'PASRatings.csv';
fprintf('Writing %s ...\n', output_csv);

fid = fopen(output_csv, 'w');
fprintf(fid, 'DyadID,ParticipantID,TrialNum,ClickPresence\n');

for i = 1:row_count
    row       = pas_rows{i};
    dyad_id   = row{1};
    part_id   = row{2};
    trial_num = row{3};
    rating    = row{4};
    fprintf(fid, '%d,%d,%d,%d\n', dyad_id, part_id, trial_num, rating);
end

fclose(fid);

csv_info = dir(output_csv);
fprintf('  -> %s (%.1f KB)\n', csv_info.name, csv_info.bytes / 1e3);

%% Compute summary statistics
if ~isempty(all_ratings)
    mean_rating   = mean(all_ratings);
    std_rating    = std(all_ratings);
    median_rating = median(all_ratings);
    mode_rating   = mode(all_ratings);
else
    mean_rating   = NaN;
    std_rating    = NaN;
    median_rating = NaN;
    mode_rating   = NaN;
end

% Distribution counts
rating_levels = 1:4;
rating_counts = zeros(1, 4);
for lv = 1:4
    rating_counts(lv) = sum(all_ratings == lv);
end

%% ========================================================================
%  Generate JSON Metadata Sidecar (BIDS-inspired)
%  ========================================================================

timestamp_str = datestr(now, 'yyyy-mm-ddTHH:MM:SS');

meta = struct();

% Dataset-level fields
meta.Name = 'Perceptual Crossing Experiment - PAS Click Presence Ratings';
meta.Description = ['Perceptual Awareness Scale (PAS) click-presence ratings ' ...
    'extracted from perceptual crossing experiment questionnaires. Each row ' ...
    'represents one participant''s self-reported awareness rating for one trial. ' ...
    'Ratings are only available for trials in which the participant pressed the button.'];

% Column specification (BIDS-style)
meta.Columns = {'DyadID', 'ParticipantID', 'TrialNum', 'ClickPresence'};

meta.DyadID = struct( ...
    'LongName', 'Dyad Identifier', ...
    'Description', ['Numeric identifier for the experimental dyad. ' ...
        'Extracted from the pce folder name (e.g., pce01230807 -> Dyad 1).']);

meta.ParticipantID = struct( ...
    'LongName', 'Participant Identifier', ...
    'Description', ['Participant number within the dyad. ' ...
        'Parsed from PAS filename (e.g., pair_01_P1_PAS*.csv -> 1).']);

meta.TrialNum = struct( ...
    'LongName', 'Trial Number', ...
    'Description', ['Trial number within the session (1-' ...
        num2str(num_trials_exp) '). Converted from 0-indexed trial_id in the raw PAS file.']);

meta.ClickPresence = struct( ...
    'LongName', 'PAS Click Presence Rating', ...
    'Description', ['Participant''s Perceptual Awareness Scale rating for the ' ...
        'clarity of experiencing the partner''s presence just before pressing ' ...
        'the button. Ordinal scale from 1 to 4.'], ...
    'Levels', struct( ...
        'x1', 'No experience of the other', ...
        'x2', 'Vague impression of the other', ...
        'x3', 'Almost clear experience of the other', ...
        'x4', 'Clear experience of the other'));

% Source data information
meta.SourceData = struct( ...
    'Description', 'PAS questionnaire CSV files from perceptual crossing experiment.', ...
    'FileFormat', '.csv', ...
    'TrialsPerDyad', num_trials_exp, ...
    'ParticipantsPerDyad', num_participants, ...
    'FolderStructure', 'pceXXYYMMDD/questionnaires/pair_XX_PY_PAS*.csv', ...
    'RawColumns', {{'trial_id', 'question_id', 'answer'}}, ...
    'TrialIdIndexing', '0-indexed in raw file; converted to 1-indexed TrialNum');

% Extraction method
meta.ExtractionMethod = struct( ...
    'Description', ['Extracts click_presence ratings from PAS questionnaire files ' ...
        'and converts from long to wide format.'], ...
    'Steps', {{'1. Discover pce* folders and locate questionnaires/ subfolder', ...
               '2. Find *PAS*.csv files and parse DyadID, ParticipantID from filename', ...
               '3. Read trial_id, question_id, answer columns', ...
               '4. Filter for question_id == ''click_presence''', ...
               '5. Convert trial_id + 1 to create 1-indexed TrialNum', ...
               '6. Sort by DyadID, ParticipantID, TrialNum and export'}});

% Data summary statistics
meta.DataSummary = struct( ...
    'NumDyads', length(dyad_ids_found), ...
    'NumPASFilesProcessed', total_files_found, ...
    'NumPASFilesSkipped', total_files_skipped, ...
    'TotalRatingsExtracted', total_ratings, ...
    'TotalExpected', total_expected, ...
    'TotalMissing', total_missing, ...
    'MeanRating', mean_rating, ...
    'StdRating', std_rating, ...
    'MedianRating', median_rating, ...
    'ModeRating', mode_rating, ...
    'RatingDistribution', struct( ...
        'x1_NoExperience', rating_counts(1), ...
        'x2_VagueImpression', rating_counts(2), ...
        'x3_AlmostClear', rating_counts(3), ...
        'x4_ClearExperience', rating_counts(4)));

% File provenance
meta.GeneratedBy = struct( ...
    'Name', 'preprocessPAS.m', ...
    'Description', 'MATLAB preprocessing script for PCE PAS click-presence extraction.', ...
    'GenerationDateTime', timestamp_str);

%% Write JSON
output_json = 'PASRatings.json';
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
fprintf('  PAS files processed:      %d\n', total_files_found);
fprintf('  PAS files skipped:        %d\n', total_files_skipped);
fprintf('  Ratings extracted:        %d / %d (%.1f%%)\n', ...
    total_ratings, total_expected, 100 * total_ratings / max(total_expected, 1));
fprintf('  Missing ratings:          %d\n', total_missing);
fprintf('\n');
fprintf('  Rating statistics:\n');
fprintf('    Mean:    %.2f\n', mean_rating);
fprintf('    Std:     %.2f\n', std_rating);
fprintf('    Median:  %.1f\n', median_rating);
fprintf('    Mode:    %d\n', mode_rating);
fprintf('\n');
fprintf('  Rating distribution:\n');
fprintf('    1 (No experience):       %d\n', rating_counts(1));
fprintf('    2 (Vague impression):    %d\n', rating_counts(2));
fprintf('    3 (Almost clear):        %d\n', rating_counts(3));
fprintf('    4 (Clear experience):    %d\n', rating_counts(4));
fprintf('\n');
fprintf('  Output files:\n');
fprintf('    %s      (data)\n', output_csv);
fprintf('    %s     (metadata sidecar)\n', output_json);
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
