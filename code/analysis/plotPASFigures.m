%% plotPASFigures.m
% =========================================================================
% Publication Figure: PAS Click-Presence Ratings Over Time
% =========================================================================
%
% Pipeline:
%   1. Load click response times and PAS ratings; merge by trial keys
%   2. Compute relative PAS-level percentages using a moving time window
%      with binomial SEM bands
%   3. Fit a binomial logistic regression to the PAS 3 vs PAS 4 contrast
%      to estimate the crossover time where PAS 4 becomes dominant
%
% Panels:
%   A — Relative percentage of each PAS level (1–4) as a function of
%       click response time, with SEM bands (moving window)
%   B — Logistic regression: P(PAS 4 | PAS 3 or 4) as a function of
%       click time, with 95% CI band and crossover time estimate
%
% INPUT:
%   ../../data/ClickTimes/ClickResponseTimes.csv
%   ../../data/PAS/PASRatings.csv
%
% OUTPUT:
%   ../../results/Fig_PAS.png (300 dpi)
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   March 2026
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  1. LOAD AND MERGE DATA
%  ========================================================================

fprintf('Loading data ...\n');

clicks = readtable('../../data/ClickTimes/ClickResponseTimes.csv', ...
    'VariableNamingRule', 'preserve');
pas = readtable('../../data/PAS/PASRatings.csv', ...
    'VariableNamingRule', 'preserve');

% Merge by trial-level keys
mergedData = innerjoin(clicks, pas, 'Keys', {'DyadID', 'ParticipantID', 'TrialNum'});

% Filter: valid click and valid PAS rating (1–4)
validIdx = mergedData.Clicked == 1 & ...
           mergedData.ClickPresence >= 1 & ...
           mergedData.ClickPresence <= 4;
filteredData = mergedData(validIdx, :);

T        = 60;
n_clicks = height(filteredData);
n_dyads  = length(unique(filteredData.DyadID));

fprintf('  %d valid click-PAS pairs from %d dyads\n', n_clicks, n_dyads);

%% ========================================================================
%  2. MOVING-WINDOW PAS RELATIVE PERCENTAGES
%  ========================================================================

windowWidth = 8;       % seconds
halfWidth   = windowWidth / 2;
stepSize    = 1;       % seconds
timeCenters = halfWidth : stepSize : (T - halfWidth);
numPoints   = length(timeCenters);

fprintf('  Moving window: %d s (step = %d s, %d centres)\n', ...
    windowWidth, stepSize, numPoints);

pasCounts = zeros(numPoints, 4);
for i = 1:numPoints
    t = timeCenters(i);
    inWindow = filteredData.ClickTime_s >= (t - halfWidth) & ...
               filteredData.ClickTime_s <  (t + halfWidth);
    windowData = filteredData(inWindow, :);
    for level = 1:4
        pasCounts(i, level) = sum(windowData.ClickPresence == level);
    end
end

totalPerBin = sum(pasCounts, 2);
totalPerBin(totalPerBin == 0) = 1;   % Avoid division by zero
relPct = (pasCounts ./ totalPerBin) * 100;

% Binomial SEM (in percentage points)
semPct = zeros(numPoints, 4);
for level = 1:4
    p = pasCounts(:, level) ./ totalPerBin;
    semPct(:, level) = sqrt(p .* (1 - p) ./ totalPerBin) * 100;
end

%% ========================================================================
%  3. LOGISTIC REGRESSION: PAS 4 vs PAS 3 CROSSOVER
%  ========================================================================

fprintf('Fitting logistic regression (PAS 4 vs PAS 3) ...\n');

% Filter for PAS 3 and PAS 4 only
targetIdx      = filteredData.ClickPresence == 3 | filteredData.ClickPresence == 4;
logreg_times   = filteredData.ClickTime_s(targetIdx);
logreg_ratings = filteredData.ClickPresence(targetIdx);
n_logreg       = length(logreg_times);

% Binary response: PAS 4 = 1, PAS 3 = 0
binaryPAS = double(logreg_ratings == 4);

% Fit binomial logistic regression
mdl = fitglm(logreg_times, binaryPAS, 'Distribution', 'binomial');

beta0 = mdl.Coefficients.Estimate(1);
beta1 = mdl.Coefficients.Estimate(2);

% Crossover: P(PAS 4) = 0.5  =>  t_cross = -beta0 / beta1
t_cross = -beta0 / beta1;

% 95% CI for crossover time via prediction interval
t_fine = linspace(min(logreg_times), max(logreg_times), 10000)';
[p_pred_fine, p_ci_fine] = predict(mdl, t_fine);

% Lower CI bound drops below 0.5 → earliest plausible crossing
idx_lower = find(p_ci_fine(:, 1) < 0.5, 1);
if ~isempty(idx_lower)
    t_cross_lo = t_fine(idx_lower);
else
    t_cross_lo = NaN;
end

% Upper CI bound drops below 0.5 → latest plausible crossing
idx_upper = find(p_ci_fine(:, 2) < 0.5, 1);
if ~isempty(idx_upper)
    t_cross_hi = t_fine(idx_upper);
else
    t_cross_hi = NaN;
end

% Order CI bounds numerically for display
t_ci_lo = min(t_cross_hi, t_cross_lo);
t_ci_hi = max(t_cross_hi, t_cross_lo);

fprintf('  Data points (PAS 3 & 4): %d\n', n_logreg);
fprintf('  Intercept (beta0):  %8.4f  (p = %.4f)\n', beta0, mdl.Coefficients.pValue(1));
fprintf('  Slope     (beta1):  %8.4f  (p = %.4f)\n', beta1, mdl.Coefficients.pValue(2));
fprintf('  Crossover time:     %.2f s\n', t_cross);
fprintf('  95%% CI:             [%.2f, %.2f] s\n', t_ci_lo, t_ci_hi);

%% ========================================================================
%  4. STYLING PARAMETERS
%  ========================================================================

% PAS level colours (1 = grey, 2 = amber, 3 = orange, 4 = black)
col_pas = [0.65 0.65 0.65; ...
           0.93 0.69 0.13; ...
           0.85 0.33 0.10; ...
           0.00 0.00 0.00];
lw_pas  = [2.0  2.0  2.0  3.0];

col_fit  = [0.80 0.15 0.15];    % Red   — logistic curve
col_data = [0.20 0.40 0.73];    % Blue  — CI band
col_grey = [0.50 0.50 0.50];    % Grey  — annotations

font_sz       = 12;
font_sz_label = 13;
font_sz_title = 14;
font_sz_annot = 11;
font_sz_panel = 18;

%% ========================================================================
%  5. CREATE FIGURE — 2 vertically stacked panels
%  ========================================================================

fprintf('Creating figure ...\n');

fig_w = 7.5;   fig_h = 8.5;
fig = figure('Units', 'inches', 'Position', [0.5 0.5 fig_w fig_h], ...
    'Color', 'w', 'PaperUnits', 'inches', ...
    'PaperSize', [fig_w fig_h], 'PaperPosition', [0 0 fig_w fig_h]);

% Layout geometry (normalised)
ml  = 0.13;  mr = 0.05;  mt = 0.04;  mb = 0.08;  gap = 0.10;
pw  = 1 - ml - mr;

h_a = 0.38;   h_b = 0.38;
y_a = 1 - mt - h_a;
y_b = y_a - gap - h_b;

%% ---- Panel A: PAS relative percentages over time -----------------------

ax_a = axes('Position', [ml y_a pw h_a]);
hold on;

for level = 1:4
    % SEM band
    fill([timeCenters, fliplr(timeCenters)], ...
         [relPct(:, level)' + semPct(:, level)', ...
          fliplr(relPct(:, level)' - semPct(:, level)')], ...
         col_pas(level, :), 'FaceAlpha', 0.13, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');

    % Line
    plot(timeCenters, relPct(:, level), 'Color', col_pas(level, :), ...
         'LineWidth', lw_pas(level), ...
         'DisplayName', sprintf('PAS %d', level));
end

hold off;

xlim([0 T]);
ylim([0 60]);
ylabel('Relative percentage of clicks (%)', 'FontSize', font_sz_label);
title(sprintf('PAS Ratings Over Time  (%d s window, N = %d clicks)', ...
    windowWidth, n_clicks), ...
    'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'on', 'TickDir', 'out');

lgd = legend('Location', 'east', 'Box', 'on', 'FontSize', font_sz_annot);
set(lgd, 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);
title(lgd, 'PAS Rating');

text(-0.10, 1.06, 'A', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ---- Panel B: Logistic regression crossover ----------------------------

ax_b = axes('Position', [ml y_b pw h_b]);
hold on;

% Raw data scatter (jittered for visibility)
jitter = (rand(size(binaryPAS)) - 0.5) * 0.08;
scatter(logreg_times, binaryPAS + jitter, 12, 'k', 'filled', ...
    'MarkerFaceAlpha', 0.10, 'HandleVisibility', 'off');

% Smooth prediction curve
t_plot = linspace(min(logreg_times), max(logreg_times), 500)';
[p_pred, p_ci] = predict(mdl, t_plot);

% 95% CI band
fill([t_plot; flipud(t_plot)], [p_ci(:, 1); flipud(p_ci(:, 2))], ...
    col_data, 'FaceAlpha', 0.20, 'EdgeColor', 'none', ...
    'DisplayName', '95% CI');

% Logistic curve
plot(t_plot, p_pred, '-', 'Color', col_fit, 'LineWidth', 2.5, ...
    'DisplayName', 'P(PAS 4)');

% 50% threshold
yline(0.5, '--k', 'LineWidth', 0.8, 'HandleVisibility', 'off');

% Crossover marker
if ~isnan(t_cross) && t_cross > 0 && t_cross < T
    xline(t_cross, ':', 'Color', col_grey, 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
    plot(t_cross, 0.5, 'o', 'MarkerSize', 9, ...
        'MarkerFaceColor', col_fit, 'MarkerEdgeColor', 'w', ...
        'LineWidth', 1.2, 'HandleVisibility', 'off');
end

% CI bounds for crossover time
if ~isnan(t_cross_lo)
    xline(t_cross_lo, ':', 'Color', col_grey, 'LineWidth', 0.8, ...
        'Alpha', 0.5, 'HandleVisibility', 'off');
end
if ~isnan(t_cross_hi)
    xline(t_cross_hi, ':', 'Color', col_grey, 'LineWidth', 0.8, ...
        'Alpha', 0.5, 'HandleVisibility', 'off');
end

hold off;

xlim([0 T]);
ylim([-0.08 1.08]);
set(gca, 'YTick', [0 0.5 1], 'YTickLabel', {'0', '0.5', '1'});
text(-1.4, 0.95, 'PAS 4', 'FontSize', 9, 'Color', col_grey, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
text(-1.4, -0.05, 'PAS 3', 'FontSize', 9, 'Color', col_grey, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');xlabel('Time (s)', 'FontSize', font_sz_label);
ylabel('P(PAS 4 | PAS 3 or 4)', 'FontSize', font_sz_label);
title('Logistic Regression: PAS 4 vs PAS 3 Crossover', ...
    'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'on', 'TickDir', 'out');

lgd_b = legend('Location', 'east', 'Box', 'on', 'FontSize', font_sz_annot);
set(lgd_b, 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);

% Annotation box
annStr = { ...
    sprintf('Crossover = %.1f s', t_cross), ...
    sprintf('95%% CI [%.1f, %.1f] s', t_ci_lo, t_ci_hi), ...
    sprintf('\\beta_0 = %.3f,  \\beta_1 = %.4f', beta0, beta1), ...
    sprintf('n = %d  (PAS 3 & 4 only)', n_logreg)};
text(0.97, 0.84, annStr, 'Units', 'normalized', 'FontSize', 9, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'Margin', 4);

text(-0.10, 1.06, 'B', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ========================================================================
%  6. SAVE
%  ========================================================================

print(fig, '../../results/Fig_PAS', '-dpng', '-r300');
fprintf('  Saved: Fig_PAS.png (300 dpi)\n');

%% ========================================================================
%  7. PRINT FULL SUMMARY
%  ========================================================================

fprintf('\n==========================================================\n');
fprintf('  FIGURE SUMMARY\n');
fprintf('==========================================================\n');
fprintf('  Panel A: PAS Relative Percentages\n');
fprintf('  --\n');
fprintf('  Valid click-PAS pairs:    %d (from %d dyads)\n', n_clicks, n_dyads);
fprintf('  Moving window width:      %d s (step = %d s)\n', windowWidth, stepSize);
fprintf('  Window centres:           %d points\n', numPoints);
fprintf('\n');
fprintf('  Panel B: Logistic Regression Crossover\n');
fprintf('  --\n');
fprintf('  Data points (PAS 3 & 4):  %d\n', n_logreg);
fprintf('  Intercept (beta0):        %8.4f  (p = %.4f)\n', beta0, mdl.Coefficients.pValue(1));
fprintf('  Slope     (beta1):        %8.4f  (p = %.4f)\n', beta1, mdl.Coefficients.pValue(2));
fprintf('  Crossover time:           %.2f s\n', t_cross);
fprintf('  95%% CI:                   [%.2f, %.2f] s\n', t_ci_lo, t_ci_hi);
fprintf('  --\n');
fprintf('  Output: Fig_PAS.png (300 dpi)\n');
fprintf('==========================================================\n');
fprintf('Done.\n');