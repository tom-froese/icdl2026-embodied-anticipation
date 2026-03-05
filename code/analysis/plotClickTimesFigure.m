%% plotClickTimesFigure.m
% =========================================================================
% Publication Figure: Click Response-Time Distribution with P(1) Fit
% =========================================================================
%
% Pipeline:
%   1. Load click response times
%   2. Estimate smooth density via kernel density estimation (KDE)
%   3. Fit the Poisson first-order term
%          f(t) = A * e * x * exp(-e*x)
%      where  x  = (t - tau) / (T - tau),   tau = onset lag,  T = 60 s.
%      Lambda is fixed at e (Euler's number).  A is estimated by scalar
%      projection at each candidate tau.  B is fixed at 0 (the click
%      density must vanish asymptotically).
%   4. Select the onset lag tau that maximises R^2 against the KDE
%      by grid search (0.1 s steps).
%
% Rationale for P(1):
%   EDA conforms to  P(0) = A0 * exp(-e*x) + B0    (the prepared state).
%   Respiration (ARU) to  P(1) = A1 * e*x*exp(-e*x) (the flux of state).
%   Click-response times are modelled as a third readout of the same
%   flux — behavioural responses are triggered by the rate of change
%   of the prepared state rather than by the state itself.
%
% Panels:
%   A. Click-time histogram with KDE overlay and best-fit P(1) curve
%   B. R^2 as a function of onset lag (with plateau analysis)
%   C. Residuals  (KDE - model)  with zero crossings
%
% INPUT:  ../../data/ClickTimes/ClickResponseTimes.csv
% OUTPUT: ../../results/Fig_ClickTimes.png (300 dpi)
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   February 2026
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  1. LOAD AND PREPARE DATA
%  ========================================================================

fprintf('Loading click data ...\n');
data = readtable('../../data/ClickTimes/ClickResponseTimes.csv');

clickTimes = data.ClickTime_s(data.Clicked == 1);
clickTimes = clickTimes(~isnan(clickTimes));
clickTimes = clickTimes(clickTimes >= 0 & clickTimes < 60);

T        = 60;                             % Trial duration (s)
lambda   = exp(1);                         % Rate parameter fixed at e
n_clicks = length(clickTimes);
n_dyads  = length(unique(data.DyadID));

fprintf('  %d valid clicks from %d dyads\n', n_clicks, n_dyads);

%% ========================================================================
%  2. KERNEL DENSITY ESTIMATE
%  ========================================================================
%  Default ksdensity (Gaussian kernel, Scott's bandwidth rule).

n_kde  = 500;
xi_kde = linspace(0, 60, n_kde);
[f_kde, xi_kde] = ksdensity(clickTimes, xi_kde, 'Support', [-0.001 60.001], 'BoundaryCorrection', 'reflection');

fprintf('  KDE: %d points in [%.1f, %.1f] s\n', ...
    n_kde, xi_kde(1), xi_kde(end));

%% ========================================================================
%  3. HISTOGRAM (for visual reference in Panel A only)
%  ========================================================================

bin_width   = 2;
edges       = 0:bin_width:T;
counts      = histcounts(clickTimes, edges);
bin_centres = edges(1:end-1) + bin_width / 2;

%% ========================================================================
%  4. SWEEP ONSET LAGS — fit P(1) to KDE
%  ========================================================================
%  Single grid search at 0.1 s resolution.

fprintf('Fitting P(1) model ...\n');

offsets = 0:0.1:15;
n_off   = length(offsets);
R2_vals = zeros(1, n_off);

for k = 1:n_off
    res         = fit_P1_kde(offsets(k), xi_kde, f_kde, T, lambda);
    R2_vals(k)  = res.R2;
end

[~, best_idx] = max(R2_vals);
best = fit_P1_kde(offsets(best_idx), xi_kde, f_kde, T, lambda);

%% ========================================================================
%  5. GOODNESS OF FIT METRICS
%  ========================================================================

n_fit    = length(best.f_fit);
n_params = 1;                  % A only  (B fixed at 0; tau selected, not jointly optimised)

% AIC and BIC (Gaussian residual assumption)
SSR = sum((best.f_fit(:) - best.y_fitted(:)).^2);
AIC = n_fit * log(SSR / n_fit) + 2 * n_params;
BIC = n_fit * log(SSR / n_fit) + n_params * log(n_fit);

% R^2 plateau analysis (1% tolerance band around peak R^2)
tol      = 0.01 * best.R2;
plateau  = R2_vals >= (best.R2 - tol);
plateau_lags = offsets(plateau);
plateau_lo   = min(plateau_lags);
plateau_hi   = max(plateau_lags);

fprintf('\n  Best onset lag:   %.1f s\n',  best.offset);
fprintf('  R²:               %.4f\n',      best.R2);
fprintf('  RMSE:             %.2e\n',       best.RMSE);
fprintf('  AIC:              %.1f\n',       AIC);
fprintf('  BIC:              %.1f\n',       BIC);
fprintf('  A (amplitude):    %.6f\n',       best.A);
fprintf('  B (baseline):     0  (fixed)\n');
fprintf('  Model peak:       %.1f s\n',     best.peak_time);
fprintf('  1%% R² plateau:   [%.1f, %.1f] s\n', plateau_lo, plateau_hi);

%% ========================================================================
%  6. RESIDUALS AND ZERO CROSSINGS
%  ========================================================================

residual_t = best.t_fit;
residual_v = best.f_fit - best.y_fitted;

% Locate zero crossings by linear interpolation
zc_times = [];
for k = 1:(length(residual_v) - 1)
    if residual_v(k) * residual_v(k+1) < 0
        frac    = abs(residual_v(k)) / ...
                  (abs(residual_v(k)) + abs(residual_v(k+1)));
        zc_time = residual_t(k) + frac * (residual_t(k+1) - residual_t(k));
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
%  7. STYLING PARAMETERS
%  ========================================================================

col_data  = [0.20 0.40 0.73];    % Blue  — histogram
col_kde   = [0.12 0.30 0.55];    % Dark blue — KDE line
col_fit   = [0.80 0.15 0.15];    % Red   — P(1) model curve
col_pos   = [0.25 0.65 0.25];    % Green — positive residuals
col_neg   = [0.80 0.27 0.27];    % Red   — negative residuals
col_grey  = [0.50 0.50 0.50];    % Grey  — onset-lag marker
col_plat  = [0.85 0.85 0.85];    % Light grey — R^2 plateau band

font_sz       = 12;
font_sz_label = 13;
font_sz_title = 14;
font_sz_annot = 11;
font_sz_panel = 18;
lw_fit        = 2.2;
lw_kde        = 1.8;

%% ========================================================================
%  8. CREATE FIGURE — 3 vertically stacked panels
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

%% ---- Panel A: Histogram + KDE + P(1) fit ----------------------------

ax_a = axes('Position', [ml y_a pw h_a]);
hold on;

% Scale factor: map density to click counts for overlay
kde_scale = n_clicks * bin_width;

% Histogram (visual reference only — the fit targets the KDE, not the bars)
bar(bin_centres, counts, 1, ...
    'FaceColor', col_data, 'FaceAlpha', 0.35, ...
    'EdgeColor', 'w', 'LineWidth', 0.5);

% KDE (the actual target of the fit)
plot(xi_kde, f_kde * kde_scale, '-', ...
    'Color', col_kde, 'LineWidth', lw_kde);

% P(1) model curve (smooth, evaluated on a fine grid)
t_mod = linspace(best.offset, T, 500);
x_mod = (t_mod - best.offset) / best.T_eff;
x_mod = max(x_mod, 1e-12);
y_mod = (best.A * lambda .* x_mod .* exp(-lambda .* x_mod)) ...
        * kde_scale;

plot(t_mod, y_mod, '-', 'Color', col_fit, 'LineWidth', lw_fit);

% Reference lines
xline(best.peak_time, ':', 'Color', col_fit,  'LineWidth', 1, 'Alpha', 0.7);
xline(best.offset,   '--', 'Color', col_grey, 'LineWidth', 1, 'Alpha', 0.7);

xlim([0 T]);
ylim([0 max(counts) * 1.15]);
ylabel('Number of clicks', 'FontSize', font_sz_label);
title('Click Response-Time Distribution', ...
    'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'on', 'TickDir', 'out');

legend({sprintf('Clicks  (n = %d)', n_clicks), 'KDE', ...
        '$P_1(x) = A\,e\,x\,\mathrm{e}^{-ex}$'}, ...
    'FontSize', font_sz_annot, 'Location', 'northwest', ...
    'Box', 'off', 'Interpreter', 'latex');

% Annotation box
text(0.97, 0.95, ...
    {sprintf('R^2 = %.3f', best.R2), ...
     sprintf('Onset lag = %.1f s', best.offset), ...
     sprintf('Peak = %.1f s', best.peak_time), ...
     sprintf('A = %.4f', best.A)}, ...
    'Units', 'normalized', 'FontSize', font_sz_annot, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'Margin', 4);

text(-0.10, 1.06, 'A', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');
hold off;

%% ---- Panel B: R² vs onset lag with plateau band ----------------------

ax_b = axes('Position', [ml y_b pw h_b]);
hold on;

% Plateau band (1% tolerance)
fill([plateau_lo plateau_hi plateau_hi plateau_lo], ...
     [0 0 max(R2_vals)*1.08 max(R2_vals)*1.08], ...
     col_plat, 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% R² curve
plot(offsets, R2_vals, '-', 'Color', col_data, 'LineWidth', 2);

% Best-lag marker
xline(best.offset, '--', 'Color', col_fit, 'LineWidth', 1.5, 'Alpha', 0.8);
plot(best.offset, best.R2, 'o', 'Color', col_fit, ...
    'MarkerSize', 8, 'MarkerFaceColor', col_fit, 'LineWidth', 1.2);

% Annotation — best lag (above curve, offset to the right)
text(best.offset + 1.5, best.R2 + 0.02, ...
    sprintf('Best: %.1f s  (R^2 = %.3f)', best.offset, best.R2), ...
    'FontSize', font_sz_annot, 'Color', col_fit, ...
    'VerticalAlignment', 'bottom');

% Annotation — plateau (below, left-aligned with band)
text(plateau_lo, 0.06, ...
    sprintf('1%% plateau: [%.1f, %.1f] s', plateau_lo, plateau_hi), ...
    'FontSize', 9, 'Color', col_grey, ...
    'VerticalAlignment', 'bottom');

xlim([0 15]);
ylim([0 max(R2_vals) * 1.08]);
xlabel('Onset lag (s)', 'FontSize', font_sz_label);
ylabel('R²', 'FontSize', font_sz_label);
title('Goodness of Fit vs. Onset Lag', ...
    'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'on', 'TickDir', 'out', ...
    'XTick', 0:2:14);

text(-0.10, 1.08, 'B', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');
hold off;

%% ---- Panel C: Residuals (KDE − model) with zero crossings -----------

ax_c = axes('Position', [ml y_c pw h_c]);
hold on;

% Filled positive / negative regions
fill_pos = max(residual_v, 0);
fill_neg = min(residual_v, 0);
area(residual_t, fill_pos, 'FaceColor', col_pos, 'FaceAlpha', 0.45, ...
    'EdgeColor', 'none');
area(residual_t, fill_neg, 'FaceColor', col_neg, 'FaceAlpha', 0.45, ...
    'EdgeColor', 'none');

% Residual trace
plot(residual_t, residual_v, '-k', 'LineWidth', 0.8);

yline(0, '-k', 'LineWidth', 0.5);
xline(best.offset, '--', 'Color', col_grey, 'LineWidth', 1, 'Alpha', 0.7);

% Mark zero crossings
rlim = max(abs(residual_v)) * 1.4;
for k = 1:length(zc_times)
    xline(zc_times(k), ':', 'Color', [0.5 0.5 0.5], ...
        'LineWidth', 0.6, 'Alpha', 0.5);
end

xlim([0 T]);
ylim([-rlim rlim]);
xlabel('Time (s)', 'FontSize', font_sz_label);
ylabel('Residual (density)', 'FontSize', font_sz_label);
title(sprintf('Residuals (KDE - P(1) Model),   RMSE = %.2e', best.RMSE), ...
    'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'on', 'TickDir', 'out');

text(0.97, 0.05, ...
    sprintf('%d zero crossings', length(zc_times)), ...
    'Units', 'normalized', 'FontSize', font_sz_annot, ...
    'HorizontalAlignment', 'right', 'Color', col_grey);

text(-0.10, 1.08, 'C', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');
hold off;

%% ========================================================================
%  9. SAVE
%  ========================================================================

print(fig, '../../results/Fig4_ClickTimes', '-dpng', '-r300');
fprintf('  Saved: Fig4_ClickTimes.png (300 dpi)\n');

%% ========================================================================
%  10. PRINT FULL SUMMARY
%  ========================================================================

fprintf('\n==========================================================\n');
fprintf('  FIGURE SUMMARY\n');
fprintf('==========================================================\n');
fprintf('  Valid clicks:        %d (from %d dyads)\n', n_clicks, n_dyads);
fprintf('  Trial duration:      %d s\n', T);
fprintf('  Lambda:              e = %.6f (fixed)\n', lambda);
fprintf('  KDE points:          %d\n', n_kde);
fprintf('  Histogram bin width: %d s (visual only)\n', bin_width);
fprintf('  --\n');
fprintf('  Best onset lag:      %.1f s\n', best.offset);
fprintf('  Effective window:    %.1f s\n', best.T_eff);
fprintf('  Model peak:          %.1f s  (lag + T_eff/e)\n', best.peak_time);
fprintf('  Amplitude (A):       %.6f\n', best.A);
fprintf('  Baseline (B):        0  (fixed)\n');
fprintf('  --\n');
fprintf('  R²:                  %.4f\n', best.R2);
fprintf('  RMSE:                %.2e\n', best.RMSE);
fprintf('  AIC:                 %.1f\n', AIC);
fprintf('  BIC:                 %.1f\n', BIC);
fprintf('  --\n');
fprintf('  1%% R² plateau:      [%.1f, %.1f] s\n', plateau_lo, plateau_hi);
fprintf('  Zero crossings:      %d\n', length(zc_times));
fprintf('==========================================================\n');
fprintf('Done.\n');

%% ========================================================================
%  LOCAL FUNCTION
%  ========================================================================

function res = fit_P1_kde(offset, xi_kde, f_kde, T, lambda)
% FIT_P1_KDE  Fit P(1) model (B = 0 fixed) to KDE density.
%
%   Model:  f(t) = A * lambda * x * exp(-lambda * x)
%           x    = (t - offset) / (T - offset)
%
%   B is fixed at zero: the click density must vanish asymptotically,
%   matching the theoretical prediction that P(1) decays to zero.
%   A is estimated by scalar projection (closed-form OLS with one
%   regressor).  Lambda is fixed (not estimated here).
%
%   Returns a struct with fields: offset, A, R2, RMSE, T_eff,
%   peak_time, t_fit, f_fit, y_fitted.

    T_eff = T - offset;
    mask  = xi_kde > offset;
    t_fit = xi_kde(mask);
    f_fit = f_kde(mask);

    if sum(mask) < 10
        res = empty_res(offset, T_eff);
        return;
    end

    % Normalised time
    x = (t_fit - offset) / T_eff;
    x = max(x, eps);

    % P(1) shape
    y_shape = lambda .* x .* exp(-lambda .* x);

    % Scalar projection:  A = (y_shape · f) / (y_shape · y_shape)
    A        = dot(y_shape(:), f_fit(:)) / dot(y_shape(:), y_shape(:));
    y_fitted = A .* y_shape;

    % Goodness of fit
    SS_res = sum((f_fit(:) - y_fitted(:)).^2);
    SS_tot = sum((f_fit(:) - mean(f_fit)).^2);
    if SS_tot > 0
        R2 = 1 - SS_res / SS_tot;
    else
        R2 = 0;
    end

    res.offset    = offset;
    res.A         = A;
    res.R2        = R2;
    res.RMSE      = sqrt(mean((f_fit(:) - y_fitted(:)).^2));
    res.T_eff     = T_eff;
    res.peak_time = offset + T_eff / lambda;
    res.t_fit     = t_fit;
    res.f_fit     = f_fit;
    res.y_fitted  = y_fitted;
end

function res = empty_res(offset, T_eff)
% EMPTY_RES  Return a valid struct when fitting is not possible.
    res.offset    = offset;
    res.A         = 0;
    res.R2        = 0;
    res.RMSE      = Inf;
    res.T_eff     = T_eff;
    res.peak_time = NaN;
    res.t_fit     = [];
    res.f_fit     = [];
    res.y_fitted  = [];
end