%% plotHapticFeedbackFigure.m
% =========================================================================
% Publication Figure: Haptic Feedback Activation with Exponential Rise Fit
% and Gap-Filling Sensitivity Analysis
% =========================================================================
%
% Pipeline:
%   1. Load preprocessed haptic feedback time series
%   2. Compute proportion of participant-trials with active feedback at
%      each time point (binary → proportion curve)
%   3. Smooth with Gaussian kernel and fit exponential rise to saturation:
%          p(t) = A * (1 - exp(-(t - t0) / tau))       for t >= t0
%   4. Sweep gap-filling thresholds to test whether bridging brief motor-
%      off gaps shifts the asymptote A toward 1/e
%
% Panels:
%   A — Smoothed proportion + exponential fit + SEM band
%   B — Residuals (smoothed − fit) showing periodic structure
%   C — Gap-filling sensitivity: asymptote A as a function of the
%       gap-filling threshold, with 95% CI band and 1/e reference
%
% INPUT:  ../../data/Haptics/HapticFeedback.csv
% OUTPUT: ../../results/Fig_HapticFeedback.png (300 dpi)
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   March 2026
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  1. LOAD DATA
%  ========================================================================

fprintf('Loading haptic feedback data ...\n');
data = readtable('../../data/Haptics/HapticFeedback.csv', ...
    'VariableNamingRule', 'preserve');
fprintf('  %d rows loaded.\n', height(data));

%% ========================================================================
%  2. BUILD TIME × SEGMENT MATRIX
%  ========================================================================

segments = unique(data(:, {'DyadID','ParticipantID','TrialNum'}), 'rows');
nSeg = height(segments);

firstIdx = data.DyadID == segments.DyadID(1) & ...
           data.ParticipantID == segments.ParticipantID(1) & ...
           data.TrialNum == segments.TrialNum(1);
timeVec = data.Time_s(firstIdx);
nTime   = length(timeVec);
fs      = round(1 / median(diff(timeVec)));

fprintf('  %d segments, %d time points (%.1f s at %d Hz)\n', ...
    nSeg, nTime, timeVec(end), fs);

% Build raw binary matrix (before any gap-filling)
sigMatrix = NaN(nTime, nSeg);
for s = 1:nSeg
    idx = data.DyadID == segments.DyadID(s) & ...
          data.ParticipantID == segments.ParticipantID(s) & ...
          data.TrialNum == segments.TrialNum(s);
    sig = data.HapticFeedback(idx);
    nSamp = min(length(sig), nTime);
    sigMatrix(1:nSamp, s) = sig(1:nSamp);
end

rawBinMatrix = double(sigMatrix >= 0.5);

%% ========================================================================
%  3. HELPER: APPLY GAP-FILLING TO A BINARY MATRIX
%  ========================================================================

    function filled = fillGaps(binMat, gapSamples)
    % FILLGAPS  Bridge OFF gaps <= gapSamples in each column of binMat.
        filled = binMat;
        if gapSamples < 1
            return;
        end
        [~, nCol] = size(binMat);
        for c = 1:nCol
            col = filled(:, c);
            if all(isnan(col)), continue; end
            d = diff(col);
            offStarts  = find(d == -1) + 1;
            onReturns  = find(d ==  1) + 1;
            if isempty(offStarts) || isempty(onReturns), continue; end
            % Drop trailing gap with no return
            if offStarts(end) > onReturns(end)
                offStarts(end) = [];
            end
            for g = 1:length(offStarts)
                nxt = onReturns(find(onReturns > offStarts(g), 1));
                if ~isempty(nxt) && (nxt - offStarts(g)) <= gapSamples
                    col(offStarts(g) : nxt - 1) = 1;
                end
            end
            filled(:, c) = col;
        end
    end

%% ========================================================================
%  4. HELPER: COMPUTE PROPORTION, SMOOTH, AND FIT EXPONENTIAL
%  ========================================================================

    function res = fitProportion(binMat, timeVec, fs, smoothWin)
    % FITPROPORTION  Proportion -> smooth -> exponential rise fit.
        nValid = sum(~isnan(binMat), 2);
        propON = sum(binMat, 2, 'omitnan') ./ nValid;
        semON  = sqrt(propON .* (1 - propON) ./ nValid);

        smoothSamp = round(smoothWin * fs);
        propSmooth = smoothdata(propON, 'gaussian', smoothSamp);

        expRise = @(p, t) p(1) .* max(1 - exp(-(t - p(3)) ./ p(2)), 0);
        p0 = [max(propSmooth), 10, 1];
        lb = [0.01, 0.5, 0];
        ub = [1.0,  60, 15];
        opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
            'MaxFunctionEvaluations', 5000, 'MaxIterations', 1000);

        [pFit, ~, residuals, ~, ~, ~, J] = ...
            lsqcurvefit(expRise, p0, timeVec, propSmooth, lb, ub, opts);

        SS_res = sum(residuals.^2);
        SS_tot = sum((propSmooth - mean(propSmooth)).^2);

        res.A       = pFit(1);
        res.tau     = pFit(2);
        res.t0      = pFit(3);
        res.R2      = 1 - SS_res / SS_tot;
        res.RMSE    = sqrt(mean(residuals.^2));
        res.propON  = propON;
        res.semON   = semON;
        res.propSmooth = propSmooth;
        res.residuals  = residuals;
        res.pFit    = pFit;

        try
            ci = nlparci(pFit, residuals, 'jacobian', J);
            res.A_ci  = ci(1,:);
            res.tau_ci = ci(2,:);
            res.t0_ci  = ci(3,:);
            res.has_ci = true;
        catch
            res.has_ci = false;
        end
    end

%% ========================================================================
%  5. BASELINE FIT (no gap-filling)
%  ========================================================================

smoothWin = 5;   % seconds
fprintf('  Smoothing: Gaussian window = %.1f s\n', smoothWin);
fprintf('Fitting baseline (no gap-filling) ...\n');

base = fitProportion(rawBinMatrix, timeVec, fs, smoothWin);

fprintf('  A = %.4f, tau = %.2f s, R^2 = %.4f\n', base.A, base.tau, base.R2);

%% ========================================================================
%  6. GAP-FILLING SWEEP
%  ========================================================================
%  Sweep gap-fill thresholds from 0 ms to 2000 ms in 10 ms steps
%  (one sample = 10 ms at 100 Hz). At each threshold, gap-fill the raw
%  binary matrix, recompute proportion, smooth, and refit.

target = 1 / exp(1);   % 1/e ~ 0.3679

gapThresh_ms = 0 : 10 : 2000;   % ms
nThresh      = length(gapThresh_ms);

sweep_A      = zeros(1, nThresh);
sweep_A_lo   = zeros(1, nThresh);
sweep_A_hi   = zeros(1, nThresh);
sweep_R2     = zeros(1, nThresh);
sweep_tau    = zeros(1, nThresh);

fprintf('Sweeping %d gap-fill thresholds (0-2000 ms) ...\n', nThresh);

for k = 1:nThresh
    gapSamp = round(gapThresh_ms(k) / 1000 * fs);
    filled  = fillGaps(rawBinMatrix, gapSamp);
    r       = fitProportion(filled, timeVec, fs, smoothWin);

    sweep_A(k)    = r.A;
    sweep_R2(k)   = r.R2;
    sweep_tau(k)  = r.tau;
    if r.has_ci
        sweep_A_lo(k) = r.A_ci(1);
        sweep_A_hi(k) = r.A_ci(2);
    else
        sweep_A_lo(k) = NaN;
        sweep_A_hi(k) = NaN;
    end
end

% --- Find where 1/e is indistinguishable from A ---
% (a) Point-estimate crossing: where A itself equals 1/e
crossIdx_A = find(sweep_A >= target, 1);
if ~isempty(crossIdx_A) && crossIdx_A > 1
    x0 = gapThresh_ms(crossIdx_A - 1);
    x1 = gapThresh_ms(crossIdx_A);
    y0 = sweep_A(crossIdx_A - 1);
    y1 = sweep_A(crossIdx_A);
    crossThresh_A_ms = x0 + (target - y0) / (y1 - y0) * (x1 - x0);
elseif ~isempty(crossIdx_A)
    crossThresh_A_ms = gapThresh_ms(1);
else
    crossThresh_A_ms = NaN;
end

% (b) CI containment: lower bound crosses target
crossIdx_lo = find(sweep_A_lo >= target, 1);
if ~isempty(crossIdx_lo) && crossIdx_lo > 1
    x0 = gapThresh_ms(crossIdx_lo - 1);
    x1 = gapThresh_ms(crossIdx_lo);
    y0 = sweep_A_lo(crossIdx_lo - 1);
    y1 = sweep_A_lo(crossIdx_lo);
    crossThresh_lo_ms = x0 + (target - y0) / (y1 - y0) * (x1 - x0);
elseif ~isempty(crossIdx_lo)
    crossThresh_lo_ms = gapThresh_ms(1);
else
    crossThresh_lo_ms = NaN;
end

% Upper CI first reaches target
crossIdx_hi = find(sweep_A_hi >= target, 1);
if ~isempty(crossIdx_hi) && crossIdx_hi > 1
    x0 = gapThresh_ms(crossIdx_hi - 1);
    x1 = gapThresh_ms(crossIdx_hi);
    y0 = sweep_A_hi(crossIdx_hi - 1);
    y1 = sweep_A_hi(crossIdx_hi);
    crossThresh_hi_ms = x0 + (target - y0) / (y1 - y0) * (x1 - x0);
elseif ~isempty(crossIdx_hi)
    crossThresh_hi_ms = gapThresh_ms(1);
else
    crossThresh_hi_ms = NaN;
end

fprintf('  Point estimate A = 1/e at:       %.1f ms\n', crossThresh_A_ms);
fprintf('  CI contains 1/e between:         [%.1f, %.1f] ms\n', ...
    crossThresh_hi_ms, crossThresh_lo_ms);

%% ========================================================================
%  7. RESIDUAL ZERO CROSSINGS
%  ========================================================================

zc_times = [];
for k = 1:(length(base.residuals) - 1)
    if base.residuals(k) * base.residuals(k+1) < 0
        frac    = abs(base.residuals(k)) / ...
                  (abs(base.residuals(k)) + abs(base.residuals(k+1)));
        zc_time = timeVec(k) + frac * (timeVec(k+1) - timeVec(k));
        zc_times(end+1) = zc_time; %#ok<SAGROW>
    end
end

fprintf('  Zero crossings:   ');
if ~isempty(zc_times)
    fprintf('%.1f', zc_times(1));
    for k = 2:length(zc_times)
        fprintf(', %.1f', zc_times(k));
    end
end
fprintf(' s  (%d total)\n', length(zc_times));

%% ========================================================================
%  8. STYLING PARAMETERS
%  ========================================================================

col_data  = [0.20 0.40 0.73];    % Blue  — data / smoothed proportion
col_fit   = [0.80 0.15 0.15];    % Red   — model curve
col_pos   = [0.25 0.65 0.25];    % Green — positive residuals
col_neg   = [0.80 0.27 0.27];    % Red   — negative residuals
col_grey  = [0.50 0.50 0.50];    % Grey  — annotations

font_sz       = 12;
font_sz_label = 13;
font_sz_title = 14;
font_sz_annot = 11;
font_sz_panel = 18;

expRise = @(p, t) p(1) .* max(1 - exp(-(t - p(3)) ./ p(2)), 0);

%% ========================================================================
%  9. CREATE FIGURE — 3 vertically stacked panels
%  ========================================================================

fprintf('Creating figure ...\n');

fig_w = 7.5;   fig_h = 11.5;
fig = figure('Units', 'inches', 'Position', [0.5 0.5 fig_w fig_h], ...
    'Color', 'w', 'PaperUnits', 'inches', ...
    'PaperSize', [fig_w fig_h], 'PaperPosition', [0 0 fig_w fig_h]);

% Layout geometry (normalised)
ml  = 0.13;  mr = 0.05;  mt = 0.04;  mb = 0.06;  gap = 0.08;
pw  = 1 - ml - mr;

h_a = 0.30;   h_b = 0.20;   h_c = 0.22;
y_a = 1 - mt - h_a;
y_b = y_a - gap - h_b;
y_c = y_b - gap - h_c;

%% ---- Panel A: Smoothed proportion + exponential fit --------------------

ax_a = axes('Position', [ml y_a pw h_a]);
hold on;

% Raw proportion (faint)
plot(timeVec, base.propON, 'Color', [col_data 0.25], 'LineWidth', 0.5, ...
    'HandleVisibility', 'off');

% SEM band on smoothed
fill([timeVec; flipud(timeVec)], ...
     [base.propSmooth + base.semON; flipud(base.propSmooth - base.semON)], ...
     col_data, 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
     'HandleVisibility', 'off');

% Smoothed line
plot(timeVec, base.propSmooth, 'Color', col_data, 'LineWidth', 2, ...
    'DisplayName', sprintf('Smoothed (%.0f s Gaussian)', smoothWin));

% Exponential fit
tPlot = linspace(0, 60, 1000)';
pPlot = expRise(base.pFit, tPlot);
plot(tPlot, pPlot, '-', 'Color', col_fit, 'LineWidth', 2.5, ...
    'DisplayName', sprintf('Exp. rise  (R^2 = %.3f)', base.R2));

% Asymptote
yline(base.A, '--', 'Color', col_fit, 'LineWidth', 1, 'Alpha', 0.6, ...
    'HandleVisibility', 'off');
text(1, base.A + 0.008, sprintf('A = %.3f', base.A), ...
    'FontSize', font_sz_annot, 'Color', col_fit);

% 1/e reference line
yline(target, ':', 'Color', col_grey, 'LineWidth', 1, 'Alpha', 0.5, ...
    'HandleVisibility', 'off');
text(57, target + 0.008, '1/e', ...
    'FontSize', 10, 'Color', col_grey, ...
    'HorizontalAlignment', 'right');

% Time constant marker
t63 = base.t0 + base.tau;
if t63 < 60
    plot(t63, base.A * (1 - exp(-1)), 'o', ...
        'MarkerSize', 8, 'MarkerFaceColor', col_fit, ...
        'MarkerEdgeColor', 'w', 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
    text(t63 + 1.5, base.A * (1 - exp(-1)) + 0.008, ...
        sprintf('\\tau = %.1f s', base.tau), ...
        'FontSize', font_sz_annot, 'Color', col_fit);
end

hold off;

xlim([0 60]);
ylim([0 max(base.propON + base.semON) * 1.1]);
ylabel('Proportion active', 'FontSize', font_sz_label);
title(sprintf('Haptic Feedback Activation  (N = %d participant-trials)', nSeg), ...
    'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'on', 'TickDir', 'out');

legend('Location', 'southeast', 'Box', 'off', 'FontSize', font_sz_annot);

% Annotation box
if base.has_ci
    annStr = { ...
        sprintf('A = %.3f  [%.3f, %.3f]', base.A, base.A_ci(1), base.A_ci(2)), ...
        sprintf('\\tau = %.2f s  [%.2f, %.2f]', base.tau, base.tau_ci(1), base.tau_ci(2)), ...
        sprintf('t_0 = %.2f s  [%.2f, %.2f]', base.t0, base.t0_ci(1), base.t0_ci(2)), ...
        sprintf('R^2 = %.4f,  RMSE = %.4f', base.R2, base.RMSE)};
else
    annStr = {sprintf('A=%.3f, \\tau=%.2f, R^2=%.4f', base.A, base.tau, base.R2)};
end
text(0.97, 0.5, annStr, 'Units', 'normalized', 'FontSize', 9, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'Margin', 4);

text(-0.10, 1.06, 'A', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ---- Panel B: Residuals (smoothed - fit) -------------------------------

ax_b = axes('Position', [ml y_b pw h_b]);
hold on;

% Filled positive / negative regions
fill_pos = max(base.residuals, 0);
fill_neg = min(base.residuals, 0);
area(timeVec, fill_pos, 'FaceColor', col_pos, 'FaceAlpha', 0.45, ...
    'EdgeColor', 'none');
area(timeVec, fill_neg, 'FaceColor', col_neg, 'FaceAlpha', 0.45, ...
    'EdgeColor', 'none');

% Residual trace
plot(timeVec, base.residuals, '-k', 'LineWidth', 0.8);
yline(0, '-k', 'LineWidth', 0.5);

% Mark zero crossings
for k = 1:length(zc_times)
    xline(zc_times(k), ':', 'Color', [0.5 0.5 0.5], ...
        'LineWidth', 0.6, 'Alpha', 0.5);
end

hold off;

xlim([0 60]);
rlim = max(abs(base.residuals)) * 1.4;
ylim([-rlim rlim]);
xlabel('Time (s)', 'FontSize', font_sz_label);
ylabel('Residual', 'FontSize', font_sz_label);
title(sprintf('Residuals (smoothed - exponential fit),   RMSE = %.4f', base.RMSE), ...
    'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'on', 'TickDir', 'out');

text(0.97, 0.05, ...
    sprintf('%d zero crossings', length(zc_times)), ...
    'Units', 'normalized', 'FontSize', font_sz_annot, ...
    'HorizontalAlignment', 'right', 'Color', col_grey);

text(-0.10, 1.08, 'B', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ---- Panel C: Gap-filling sensitivity — A vs threshold -----------------

ax_c = axes('Position', [ml y_c pw h_c]);
hold on;

% 95% CI band
fill([gapThresh_ms, fliplr(gapThresh_ms)], ...
     [sweep_A_hi, fliplr(sweep_A_lo)], ...
     col_data, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

% A estimate
plot(gapThresh_ms, sweep_A, '-', 'Color', col_data, 'LineWidth', 2);

% 1/e reference
yline(target, '-', '1/e', ...
    'Color', col_fit, 'LineWidth', 1.5, ...
    'LabelHorizontalAlignment', 'left', 'FontSize', font_sz_annot, ...
    'LabelVerticalAlignment', 'bottom');

% Point-estimate crossing
if ~isnan(crossThresh_A_ms)
    xline(crossThresh_A_ms, ':', ...
        'Color', col_grey, 'LineWidth', 1.5);
    plot(crossThresh_A_ms, target, 'o', ...
        'MarkerSize', 9, 'MarkerFaceColor', col_fit, ...
        'MarkerEdgeColor', 'w', 'LineWidth', 1.2);

    text(crossThresh_A_ms + 30, target + 0.02, ...
        sprintf('A = 1/e at %.0f ms', crossThresh_A_ms), ...
        'FontSize', 10, 'Color', col_fit, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'bottom');
end

% CI containment window
if ~isnan(crossThresh_hi_ms) && ~isnan(crossThresh_lo_ms)
    text(0.97, 0.05, ...
        sprintf('CI includes 1/e:  [%.0f, %.0f] ms', ...
            crossThresh_hi_ms, crossThresh_lo_ms), ...
        'Units', 'normalized', 'FontSize', 9, ...
        'HorizontalAlignment', 'right', 'Color', col_grey);
end

hold off;

xlim([0 2000]);
xlabel('Gap-fill threshold (ms)', 'FontSize', font_sz_label);
ylabel('Asymptote A', 'FontSize', font_sz_label);
title('Sensitivity of Asymptote to Perceptual Resolution (Gap-Filling)', ...
    'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'on', 'TickDir', 'out');

text(-0.10, 1.08, 'C', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ========================================================================
%  10. SAVE
%  ========================================================================

print(fig, '../../results/Fig_HapticFeedback', '-dpng', '-r300');
fprintf('  Saved: Fig_HapticFeedback.png (300 dpi)\n');

%% ========================================================================
%  11. PRINT FULL SUMMARY
%  ========================================================================

fprintf('\n==========================================================\n');
fprintf('  FIGURE SUMMARY\n');
fprintf('==========================================================\n');
fprintf('  Model:  p(t) = A * (1 - exp(-(t - t0) / tau))\n');
fprintf('  Smoothing:            Gaussian, %.1f s window\n', smoothWin);
fprintf('  --\n');
fprintf('  Participant-trials:   %d\n', nSeg);
fprintf('  Time points:          %d\n', nTime);
fprintf('  Sampling rate:        %d Hz\n', fs);
fprintf('  --\n');
fprintf('  Asymptote (A):        %.4f', base.A);
if base.has_ci, fprintf('   95%% CI [%.4f, %.4f]', base.A_ci(1), base.A_ci(2)); end
fprintf('\n');
fprintf('  Time constant (tau):  %.2f s', base.tau);
if base.has_ci, fprintf('  95%% CI [%.2f, %.2f]', base.tau_ci(1), base.tau_ci(2)); end
fprintf('\n');
fprintf('  Onset delay (t0):     %.2f s', base.t0);
if base.has_ci, fprintf('  95%% CI [%.2f, %.2f]', base.t0_ci(1), base.t0_ci(2)); end
fprintf('\n');
fprintf('  --\n');
fprintf('  R^2:                  %.4f\n', base.R2);
fprintf('  RMSE:                 %.4f\n', base.RMSE);
fprintf('  63%% rise time:        %.1f s  (t0 + tau)\n', base.t0 + base.tau);
fprintf('  95%% rise time:        %.1f s  (t0 + 3*tau)\n', base.t0 + 3*base.tau);
fprintf('  Zero crossings:       %d\n', length(zc_times));
fprintf('  --\n');
fprintf('  1/e target:                  %.4f\n', target);
fprintf('  A = 1/e at gap threshold:    %.1f ms\n', crossThresh_A_ms);
if ~isnan(crossThresh_hi_ms) && ~isnan(crossThresh_lo_ms)
    fprintf('  CI contains 1/e:             [%.0f, %.0f] ms\n', ...
        crossThresh_hi_ms, crossThresh_lo_ms);
end
fprintf('==========================================================\n');
fprintf('Done.\n');
