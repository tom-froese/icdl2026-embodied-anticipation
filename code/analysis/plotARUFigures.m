%% plotARUFigures.m
% =========================================================================
% Grand Average Respiration: 3-Panel Vertical Stack
% =========================================================================
%
% Panels:
%   A. Rest grand average -ARU with smoothed trend
%   B. Task grand average -ARU with P(1) fit (variable lag time)
%   C. Task residuals (smoothed trend minus P(1) fit, full trial)
%
% P(1) model with offset (lambda = e, fixed horizon):
%   f(x) = A0 * e * x * exp(-e * x) + B
%   where x = (t - t0) / T,  T = D - t0  (fixed post-onset window)
%
% This mirrors the EDA P(0) approach: time is normalised to [0,1]
% over the remaining trial window after onset, with lambda = e.
%
% Lambda = e (fixed). Only A0 and B are free parameters.
% The optimal lag time t0 is found by grid search (0.1s resolution)
% over [1, 10] seconds, selecting the t0 that maximises R².
%
% =========================================================================

%% Parameters
fs = 25;
task_duration = 60;
rest_duration = 180;
task_samples = task_duration * fs;
rest_samples = rest_duration * fs;

smooth_window_s = 10;
smooth_samples = smooth_window_s * fs;

% P(1) model: f(x) = A0 * e * x * exp(-e * x) + B,  x = (t-t0)/(D-t0)
% Lambda = e (fixed), T = D - t0 (fixed post-onset window).
% Grid search over t0; for each t0, fit A0 and B.
e_const = exp(1);

% Lag time search range
t0_search = 1.0 : 0.1 : 10.0;

% Figure style
color_task = [0.8, 0.2, 0.2];
color_rest = [0.2, 0.3, 0.7];
color_task_dark = [0.6, 0.0, 0.0];
color_rest_dark = [0.1, 0.1, 0.5];
color_model = [0.1, 0.1, 0.1];
shade_alpha = 0.2;
resid_pos = [0.4, 0.8, 0.4];
resid_neg = [0.9, 0.5, 0.5];
font_size = 11;
title_size = 13;

%% Load data
fprintf('Loading data ...\n');
task_data = readtable('../../data/ARU/ARU_Task_Preprocessed.csv');
rest_data = readtable('../../data/ARU/ARU_Rest_Preprocessed.csv');
fprintf('  Task: %d rows, Rest: %d rows\n', height(task_data), height(rest_data));

%% Per-participant centering
task_data.PartKey = task_data.DyadID * 10 + task_data.ParticipantID;
rest_data.PartKey = rest_data.DyadID * 10 + rest_data.ParticipantID;

unique_parts = unique(task_data.PartKey);
n_participants = length(unique_parts);

% Per-participant centering WITHIN each condition separately.
% This removes between-participant amplitude differences while preserving
% the within-trial temporal structure (including P(1) dynamics).
% Using task+rest combined would bias task centering because the rest
% baseline sits lower, lifting the tail of the task P(1) curve.

task_part_means = zeros(height(task_data), 1);
rest_part_means = zeros(height(rest_data), 1);
for i = 1:n_participants
    pk = unique_parts(i);
    t_mask = task_data.PartKey == pk;
    r_mask = rest_data.PartKey == pk;
    task_part_means(t_mask) = nanmean(task_data.RESP_ARU(t_mask));
    rest_part_means(r_mask) = nanmean(rest_data.RESP_ARU(r_mask));
end

task_centered = task_data.RESP_ARU - task_part_means;
rest_centered = rest_data.RESP_ARU - rest_part_means;

%% Reshape into matrices
task_segments = unique(task_data(:, {'DyadID','ParticipantID','TrialNum'}), 'rows');
n_task_seg = height(task_segments);
task_matrix = NaN(task_samples, n_task_seg);

for i = 1:n_task_seg
    mask = task_data.DyadID == task_segments.DyadID(i) & ...
           task_data.ParticipantID == task_segments.ParticipantID(i) & ...
           task_data.TrialNum == task_segments.TrialNum(i);
    vals = task_centered(mask);
    if length(vals) == task_samples
        task_matrix(:, i) = vals;
    end
end

rest_segments = unique(rest_data(:, {'DyadID','ParticipantID','RestNum'}), 'rows');
n_rest_seg = height(rest_segments);
rest_matrix = NaN(rest_samples, n_rest_seg);

for i = 1:n_rest_seg
    mask = rest_data.DyadID == rest_segments.DyadID(i) & ...
           rest_data.ParticipantID == rest_segments.ParticipantID(i) & ...
           rest_data.RestNum == rest_segments.RestNum(i);
    vals = rest_centered(mask);
    if length(vals) == rest_samples
        rest_matrix(:, i) = vals;
    end
end

%% Grand average, SEM, smoothed
task_mean = nanmean(task_matrix, 2);
task_sem  = nanstd(task_matrix, 0, 2) ./ sqrt(sum(~isnan(task_matrix), 2));
task_time = (0:task_samples-1)' / fs;

rest_mean = nanmean(rest_matrix, 2);
rest_sem  = nanstd(rest_matrix, 0, 2) ./ sqrt(sum(~isnan(rest_matrix), 2));
rest_time = (0:rest_samples-1)' / fs;

task_smooth = movmean(task_mean, smooth_samples);
rest_smooth = movmean(rest_mean, smooth_samples);

%% Fit P(1) with fixed horizon T = D - t0 (normalised time)
% For each candidate t0, normalise time to x = (t-t0)/(D-t0) in [0,1],
% then fit A0 and B in: f(x) = A0 * e * x * exp(-e*x) + B.

best_R2  = -Inf;
best_t0  = NaN;
best_A0  = NaN;
best_B   = NaN;

fprintf('\nFitting P(1) with lambda=e, fixed horizon T=D-t0 ...\n');

% Define P(1) model function for lsqcurvefit:
%   params = [A0, B],  xdata = normalised time x in [0,1]
p1_func = @(params, x) params(1) .* e_const .* x .* exp(-e_const .* x) + params(2);

opts = optimoptions('lsqcurvefit', 'Display', 'off');
lb = [-Inf, -Inf];
ub = [ Inf,  Inf];

for ti = 1:length(t0_search)
    t0_cand = t0_search(ti);
    T_eff = task_duration - t0_cand;           % fixed window length
    if T_eff < 10; continue; end
    
    % Fit region: t >= t0
    fit_mask = task_time >= t0_cand;
    t_fit = task_time(fit_mask);
    y_fit = task_smooth(fit_mask);
    x_norm = (t_fit - t0_cand) / T_eff;       % normalised to [0, 1]
    
    % Initial guesses from data
    y_late = mean(y_fit(end-100:end));         % Asymptotic level -> B
    y_peak = max(y_fit);
    A0_guess = (y_peak - y_late) / (e_const * (1/e_const) * exp(-1));
    
    x0 = [A0_guess, y_late];
    
    try
        params = lsqcurvefit(p1_func, x0, x_norm, y_fit, lb, ub, opts);
    catch
        continue;
    end
    
    y_pred = p1_func(params, x_norm);
    SS_res = sum((y_fit - y_pred).^2);
    SS_tot = sum((y_fit - mean(y_fit)).^2);
    R2 = 1 - SS_res / SS_tot;
    
    if R2 > best_R2
        best_R2  = R2;
        best_t0  = t0_cand;
        best_A0  = params(1);
        best_B   = params(2);
    end
end

best_T = task_duration - best_t0;              % effective window length
best_k = e_const / best_T;                    % effective rate in clock time

fprintf('  Best fit: t0=%.1fs, A0=%.1f, B=%.1f, T=%.1fs, R²=%.3f\n', ...
    best_t0, best_A0, best_B, best_T, best_R2);
fprintf('  lambda = e (fixed), k_eff = e/T = e/%.1f = %.4f\n', best_T, best_k);
fprintf('  P(1) peak at t0 + T/e = %.1f + %.1f = %.1fs\n', ...
    best_t0, best_T/e_const, best_t0 + best_T/e_const);

% Generate P(1) prediction only for t >= t0
fit_mask = task_time >= best_t0;
x_norm_full = (task_time(fit_mask) - best_t0) / best_T;
p1_fit_region = best_A0 .* e_const .* x_norm_full .* exp(-e_const .* x_norm_full) + best_B;

% Residuals only in fit region (t >= t0)
task_residuals_fit = task_smooth(fit_mask) - p1_fit_region;
RMSE_fit = sqrt(mean(task_residuals_fit.^2));

fprintf('  RMSE (fit region): %.2f ARU\n', RMSE_fit);

%% ========================================================================
%  FIGURE: 3-panel vertical stack
%  ========================================================================

fig = figure('Position', [100, 50, 900, 750], 'Color', 'w');

% ---- Panel A: Rest Grand Average ----
ax1 = subplot(3, 1, 1);
hold on;

fill([rest_time; flipud(rest_time)], ...
     [rest_mean + rest_sem; flipud(rest_mean - rest_sem)], ...
     color_rest, 'FaceAlpha', shade_alpha, 'EdgeColor', 'none');

plot(rest_time, rest_mean, 'Color', [color_rest, 0.4], 'LineWidth', 0.5);
h_smooth = plot(rest_time, rest_smooth, 'Color', color_rest_dark, 'LineWidth', 2);

yline(0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
hold off;

ylabel('ARU', 'FontSize', font_size);
title(sprintf('A.  Rest: Grand Average (N=%d, %d segments)', ...
    n_participants, n_rest_seg), 'FontSize', title_size, 'FontWeight', 'bold');
set(gca, 'FontSize', font_size, 'Box', 'on');
xlim([0, rest_duration]);
xlabel('Time (s)', 'FontSize', font_size);

legend(h_smooth, sprintf('Smoothed (%ds)', smooth_window_s), ...
    'Location', 'northeast', 'FontSize', 9, 'Box', 'off');

% ---- Panel B: Task Grand Average with P(1) fit ----
ax2 = subplot(3, 1, 2);
hold on;

% SEM shading
fill([task_time; flipud(task_time)], ...
     [task_mean + task_sem; flipud(task_mean - task_sem)], ...
     color_task, 'FaceAlpha', shade_alpha, 'EdgeColor', 'none');

% Raw grand average
plot(task_time, task_mean, 'Color', [color_task, 0.4], 'LineWidth', 0.5);

% Smoothed trend
h_smooth = plot(task_time, task_smooth, 'Color', color_task_dark, 'LineWidth', 2);

% P(1) fit (dashed black) — only from t0 onward
fit_time = task_time(fit_mask);
h_model = plot(fit_time, p1_fit_region, '--', 'Color', color_model, 'LineWidth', 2);

% Peak marker on P(1) fit
[peak_val, peak_idx] = max(p1_fit_region);
peak_t = fit_time(peak_idx);
plot(peak_t, peak_val, 'v', 'Color', color_model, 'MarkerSize', 8, ...
    'MarkerFaceColor', color_model, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
text(peak_t + 1.5, peak_val + 2, sprintf('%.1fs', peak_t), ...
    'FontSize', 9, 'Color', [0.2 0.2 0.2], 'FontWeight', 'bold');

% Offset line
xline(best_t0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);

yline(0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
hold off;

ylabel('ARU', 'FontSize', font_size);
title(sprintf('B.  Task: Grand Average with P(1) fit (N=%d, %d trials)', ...
    n_participants, n_task_seg), 'FontSize', title_size, 'FontWeight', 'bold');
set(gca, 'FontSize', font_size, 'Box', 'on', 'XTickLabel', []);
xlim([0, task_duration]);

% Match y-axes between panels A and B
yl_rest = get(ax1, 'YLim');
yl_task = get(ax2, 'YLim');
common_yl = [min(yl_rest(1), yl_task(1)), max(yl_rest(2), yl_task(2))];
ylim(ax1, common_yl);
ylim(ax2, common_yl);

legend([h_smooth, h_model], ...
    {sprintf('Smoothed (%ds)', smooth_window_s), ...
     sprintf('P(1): A_0=%.1f, B=%.1f, \\lambda=e, T=%.1fs, t_0=%.1fs, R^2=%.3f', ...
             best_A0, best_B, best_T, best_t0, best_R2)}, ...
    'Location', 'southeast', 'FontSize', 9, 'Box', 'off');

% ---- Panel C: Task Residuals (post-onset only) ----
ax3 = subplot(3, 1, 3);
hold on;

% Filled residual areas (post-onset only)
pos_r = max(task_residuals_fit, 0);
neg_r = min(task_residuals_fit, 0);
area(fit_time, pos_r, 'FaceColor', resid_pos, ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none');
area(fit_time, neg_r, 'FaceColor', resid_neg, ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Residual line
plot(fit_time, task_residuals_fit, 'Color', color_task_dark, 'LineWidth', 1);

% Zero and offset lines
yline(0, '-', 'Color', [0.5, 0, 0], 'LineWidth', 1);
xline(best_t0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);

% Zero crossings with minimum gap filter (matching EDA approach)
min_zc_gap = 5;  % seconds
signs = sign(task_residuals_fit);
zc_idx = find(signs(1:end-1) .* signs(2:end) < 0);
zc_all = zeros(1, length(zc_idx));
for z = 1:length(zc_idx)
    i = zc_idx(z);
    zc_all(z) = fit_time(i) - task_residuals_fit(i) * ...
        (fit_time(i+1) - fit_time(i)) / (task_residuals_fit(i+1) - task_residuals_fit(i));
end
% Filter: keep only crossings separated by >= min_zc_gap
if ~isempty(zc_all)
    zc_filt = zc_all(1);
    for z = 2:length(zc_all)
        if zc_all(z) - zc_filt(end) >= min_zc_gap
            zc_filt(end+1) = zc_all(z);  %#ok<AGROW>
        end
    end
    ylims = [min(task_residuals_fit)*1.1, max(task_residuals_fit)*1.1];
    for z = 1:length(zc_filt)
        xline(zc_filt(z), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
        text(zc_filt(z), ylims(2)*0.85, sprintf('%.1fs', zc_filt(z)), ...
            'FontSize', font_size - 2, 'HorizontalAlignment', 'center', ...
            'Color', [0.3 0.3 0.3]);
    end
end

hold off;

xlabel('Time (s)', 'FontSize', font_size);
ylabel('Residual (ARU)', 'FontSize', font_size);
title(sprintf('C.  Task: Residuals (RMSE = %.2f ARU)', ...
    RMSE_fit), 'FontSize', title_size, 'FontWeight', 'bold');
set(gca, 'FontSize', font_size, 'Box', 'on');
xlim([0, task_duration]);

% Link x-axes for task panels
linkaxes([ax2, ax3], 'x');

%% Save
print(fig, '../../results/Fig7_ARU', '-dpng', '-r300');
fprintf('\nSaved Fig7_ARU.png\n');

fprintf('\n==========================================================\n');
fprintf('  P(1) Model: f(x) = A0 * e * x * exp(-e*x) + B\n');
fprintf('  where x = (t - t0) / T,  T = D - t0  (D = %ds)\n', task_duration);
fprintf('==========================================================\n');
fprintf('  Offset t0:    %.1f s\n', best_t0);
fprintf('  Amplitude:    A0 = %.1f ARU\n', best_A0);
fprintf('  Baseline:     B  = %.1f ARU\n', best_B);
fprintf('  Lambda:       e (fixed = %.4f)\n', e_const);
fprintf('  Time-scale:   T = D - t0 = %.1f s\n', best_T);
fprintf('  Eff. rate:    k = e/T = %.4f s^-1\n', best_k);
fprintf('  Peak at:      t0 + T/e = %.1f + %.1f = %.1fs\n', ...
    best_t0, best_T/e_const, best_t0 + best_T/e_const);
fprintf('  R² (fit region):   %.3f\n', best_R2);
fprintf('  RMSE (fit region): %.2f ARU\n', RMSE_fit);
fprintf('  Zero crossings (>%ds gap):', min_zc_gap);
if exist('zc_filt', 'var')
    for z = 1:length(zc_filt)
        fprintf(' %.1fs', zc_filt(z));
    end
end
fprintf('\n');
fprintf('==========================================================\n');
fprintf('Done.\n');