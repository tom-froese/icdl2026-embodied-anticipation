%% plotEDAFigures.m
% =========================================================================
% Consolidated EDA Figures for IEEE ICDL 2026
% =========================================================================
%
% Generates three figures sized for IEEE double-column format:
%
% FIGURE 1: Rest condition (2 vertically stacked panels)
%   A. Rest grand average + P(0) fit
%   B. Rest residuals with green/red shading and zero crossings
%
% FIGURE 2: Onset lag analysis (2 panels, for reference only — not
%           included in manuscript)
%   A. R² vs. onset lag (both conditions)
%   B. Free rate vs. e/(60-tau) crossing
%
% FIGURE 3: Task condition — optimised P(0) fit (2 vertically stacked panels)
%   A. Task EDA with P(0) fit at optimal onset lag
%   B. Residuals of the optimised fit
%
% NOTE ON CENTERING:
%   Unlike respiration (ARU), EDA is not participant-mean-centered before
%   grand averaging. This is justified because the P(0) model includes a
%   free tonic baseline parameter B, which absorbs between-participant
%   differences in absolute EDA level. The two-level aggregation (within-
%   participant trial averaging, then across-participant grand averaging)
%   gives each participant equal weight regardless of absolute EDA level.
%   EDA in µS is a calibrated physiological measure with a meaningful
%   zero, and centering would destroy the interpretability of A₀ and B
%   in absolute units.
%
% INPUT:  ../../data/EDA/EDA_Task_Preprocessed.csv
%         ../../data/EDA/EDA_Rest_Preprocessed.csv
% OUTPUT: ../../results/Fig1_EDA_Rest.png
%         ../../results/Fig2_EDA_OnsetLag.png
%         ../../results/Fig3_EDA_Task.png (all 600 dpi)
% =========================================================================

%% Parameters
fs = 25;
task_duration = 60;
rest_duration = 180;
task_samples = task_duration * fs;
rest_samples = rest_duration * fs;

%% Load and compute grand averages
fprintf('Loading data ...\n');

% --- Task ---
task_data = readtable('../../data/EDA/EDA_Task_Preprocessed.csv');
task_data.PID = task_data.DyadID * 10 + task_data.ParticipantID;
task_pids = unique(task_data.PID);

task_participant_matrix = NaN(length(task_pids), task_samples);
for i = 1:length(task_pids)
    pid = task_pids(i);
    pdata = task_data(task_data.PID == pid, :);
    trials = unique(pdata.TrialNum);
    trial_sum = zeros(1, task_samples); trial_count = 0;
    for t = 1:length(trials)
        eda = pdata.EDA_uS(pdata.TrialNum == trials(t));
        if length(eda) == task_samples
            trial_sum = trial_sum + eda'; trial_count = trial_count + 1;
        end
    end
    if trial_count > 0
        task_participant_matrix(i, :) = trial_sum / trial_count;
    end
end
task_participant_matrix = task_participant_matrix(~any(isnan(task_participant_matrix), 2), :);
n_task = size(task_participant_matrix, 1);
task_time = (0:task_samples-1) / fs;
task_mean = mean(task_participant_matrix, 1);
task_sem  = std(task_participant_matrix, 0, 1) / sqrt(n_task);

% --- Rest ---
rest_data = readtable('../../data/EDA/EDA_Rest_Preprocessed.csv');
rest_data.PID = rest_data.DyadID * 10 + rest_data.ParticipantID;
rest_pids = unique(rest_data.PID);

rest_participant_matrix = NaN(length(rest_pids), rest_samples);
for i = 1:length(rest_pids)
    pid = rest_pids(i);
    pdata = rest_data(rest_data.PID == pid, :);
    periods = unique(pdata.RestNum);
    period_sum = zeros(1, rest_samples); period_count = 0;
    for r = 1:length(periods)
        eda = pdata.EDA_uS(pdata.RestNum == periods(r));
        if length(eda) == rest_samples
            period_sum = period_sum + eda'; period_count = period_count + 1;
        end
    end
    if period_count > 0
        rest_participant_matrix(i, :) = period_sum / period_count;
    end
end
rest_participant_matrix = rest_participant_matrix(~any(isnan(rest_participant_matrix), 2), :);
n_rest = size(rest_participant_matrix, 1);
rest_time = (0:rest_samples-1) / fs;
rest_mean = mean(rest_participant_matrix, 1);
rest_sem  = std(rest_participant_matrix, 0, 1) / sqrt(n_rest);

fprintf('  Task: N=%d,  Rest: N=%d\n', n_task, n_rest);

%% Fit P(0) models
fprintf('Fitting P(0) models ...\n');

% --- Rest: full window, lambda=e ---
x_rest = rest_time / rest_duration;
[A0_rest, B_rest, R2_rest, RMSE_rest, yfit_rest] = fit_p0_lambda_e(x_rest, rest_mean);
resid_rest = (rest_mean - yfit_rest) * 1000;  % mµS

% --- Task: full window, lambda=e (no offset) ---
x_task_full = task_time / task_duration;
[A0_task_full, B_task_full, R2_task_full, RMSE_task_full, yfit_task_full] = ...
    fit_p0_lambda_e(x_task_full, task_mean);
resid_task_full = (task_mean - yfit_task_full) * 1000;

% --- Task: onset lag sweep ---
tau_sweep = 0:0.1:15;
n_tau = length(tau_sweep);
R2_task_sweep = NaN(1, n_tau);
RMSE_task_sweep = NaN(1, n_tau);
A0_task_sweep = NaN(1, n_tau);
B_task_sweep  = NaN(1, n_tau);
R2_rest_sweep = NaN(1, n_tau);

p0_model_free = @(params, t) params(1) * exp(-params(2) * t) + params(3);
opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
    'MaxIterations', 10000, 'MaxFunctionEvaluations', 30000, ...
    'TolFun', 1e-12, 'TolX', 1e-12);

k_free_at_tau = NaN(1, n_tau);
k_eff_at_tau  = NaN(1, n_tau);

for s = 1:n_tau
    tau = tau_sweep(s);
    tau_samp = round(tau * fs);
    
    % Task
    T_eff_task = task_duration - tau;
    if T_eff_task < 10; continue; end
    idx_t = (tau_samp + 1):task_samples;
    y_t = task_mean(idx_t);
    x_t = (task_time(idx_t) - tau) / T_eff_task;
    [A0, B, R2, RMSE, ~] = fit_p0_lambda_e(x_t, y_t);
    R2_task_sweep(s) = R2;
    RMSE_task_sweep(s) = RMSE;
    A0_task_sweep(s) = A0;
    B_task_sweep(s) = B;
    k_eff_at_tau(s) = exp(1) / T_eff_task;
    
    % Free rate at this tau (for Panel 2B)
    t_clock = task_time(idx_t) - tau;
    A0_init = y_t(1) - y_t(end);
    B_init = y_t(end);
    p0 = [A0_init, 0.05, B_init];
    lb = [0, 0.0001, 0]; ub = [50, 1, 50];
    params = lsqcurvefit(p0_model_free, p0, t_clock, y_t, lb, ub, opts);
    k_free_at_tau(s) = params(2);
    
    % Rest
    T_eff_rest = rest_duration - tau;
    if T_eff_rest < 10; continue; end
    idx_r = (tau_samp + 1):rest_samples;
    y_r = rest_mean(idx_r);
    x_r = (rest_time(idx_r) - tau) / T_eff_rest;
    [~, ~, R2, ~, ~] = fit_p0_lambda_e(x_r, y_r);
    R2_rest_sweep(s) = R2;
end

% Find optimal tau for task
[R2_task_max, idx_task_max] = max(R2_task_sweep);
tau_opt = tau_sweep(idx_task_max);
T_eff_opt = task_duration - tau_opt;

% Compute task fit at optimal tau
tau_samp_opt = round(tau_opt * fs);
idx_opt = (tau_samp_opt + 1):task_samples;
x_opt = (task_time(idx_opt) - tau_opt) / T_eff_opt;
yfit_task_opt = A0_task_sweep(idx_task_max) * exp(-exp(1) * x_opt) + B_task_sweep(idx_task_max);

% Compute residuals for the optimised task fit (over full trial for display)
% Pre-onset region: residual = data - B (since the model is undefined there)
resid_task_opt_post = (task_mean(idx_opt) - yfit_task_opt) * 1000;  % mµS

% Find rate crossing point
k_diff = k_free_at_tau - k_eff_at_tau;
valid = ~isnan(k_diff);
k_diff_v = k_diff(valid);
tau_v = tau_sweep(valid);
cross_idx = find(k_diff_v(1:end-1) .* k_diff_v(2:end) < 0, 1, 'first');
if ~isempty(cross_idx)
    tau_cross = tau_v(cross_idx) - k_diff_v(cross_idx) * ...
        (tau_v(cross_idx+1) - tau_v(cross_idx)) / ...
        (k_diff_v(cross_idx+1) - k_diff_v(cross_idx));
    k_cross = interp1(tau_v, k_eff_at_tau(valid), tau_cross);
else
    tau_cross = tau_opt;
    k_cross = k_eff_at_tau(idx_task_max);
end

fprintf('  Rest:  R2 = %.4f, A0 = %.3f, B = %.3f\n', R2_rest, A0_rest, B_rest);
fprintf('  Task (no lag):  R2 = %.4f, A0 = %.3f, B = %.3f\n', R2_task_full, A0_task_full, B_task_full);
fprintf('  Task (tau=%.1fs): R2 = %.4f, A0 = %.3f, B = %.3f\n', ...
    tau_opt, R2_task_max, A0_task_sweep(idx_task_max), B_task_sweep(idx_task_max));
fprintf('  Rate crossing: tau_cross = %.2fs, k = %.4f\n', tau_cross, k_cross);

%% ========================================================================
%  Shared styling — INCREASED FONT SIZES for IEEE double-column
%  ========================================================================

col_rest = [0.12 0.35 0.65];
col_task = [0.85 0.32 0.1];
col_fit  = [0.1 0.1 0.1];
col_pos  = [0.25 0.65 0.25];  % green for positive residuals
col_neg  = [0.75 0.25 0.25];  % red for negative residuals

% Create at 2x target size (7 in wide) with larger fonts.
% Scale down to 3.5 in in Word — fonts will remain legible.
fig_width_in    = 7.0;
font_size       = 14;   % tick labels, general text
font_size_label = 15;   % axis labels
font_size_title = 16;   % panel titles (Rest / Task)
font_size_annot = 12;   % annotation boxes (R², A₀, B)
font_size_legend = 11;  % legend text
font_size_zc    = 12;   % zero-crossing labels
font_size_panel = 18;   % panel letters (A, B, ...)
line_width_data = 1.5;
line_width_fit  = 2.0;

% Minimum gap (seconds) between reported zero crossings
min_zc_gap = 5;

% Shared y-axis ranges across conditions
eda_ymin = min([min(rest_mean - rest_sem), min(task_mean - task_sem)]) - 0.1;
eda_ymax = max([max(rest_mean + rest_sem), max(task_mean + task_sem)]) + 0.1;

%% ========================================================================
%  FIGURE 1: REST — P(0) Fit and Residuals
%  ========================================================================

fprintf('Creating Figure 1 (Rest) ...\n');

fig1_height_in = 8.0;  % 2 stacked panels
fig1 = figure('Units', 'inches', 'Position', [0.5 0.5 fig_width_in fig1_height_in], ...
    'Color', 'w', 'PaperUnits', 'inches', ...
    'PaperSize', [fig_width_in fig1_height_in], ...
    'PaperPosition', [0 0 fig_width_in fig1_height_in]);

% Layout
gap = 0.08;
margin_left = 0.14;
margin_right = 0.03;
margin_top = 0.04;
margin_bot = 0.10;
panel_w = 1 - margin_left - margin_right;
h_main  = 0.42;
h_resid = 0.32;
start_y = 1 - margin_top - h_main;

resid_absmax_rest = max(abs(resid_rest)) * 1.15;

% ---- Panel A: Rest EDA + P(0) fit ----
ax1a = axes('Position', [margin_left, start_y, panel_w, h_main]);
hold on;
fill([rest_time, fliplr(rest_time)], ...
    [rest_mean + rest_sem, fliplr(rest_mean - rest_sem)], ...
    col_rest, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(rest_time, rest_mean, '-', 'Color', col_rest, 'LineWidth', line_width_data);
plot(rest_time, yfit_rest, '-', 'Color', col_fit, 'LineWidth', line_width_fit);
ylabel('EDA (\muS)', 'FontSize', font_size_label);
set(gca, 'FontSize', font_size, 'Box', 'on', 'TickDir', 'out');
title('Rest', 'FontSize', font_size_title, 'FontWeight', 'bold');
ylim([eda_ymin, eda_ymax]);

% Annotation box
text(0.97, 0.95, {sprintf('R^2 = %.3f', R2_rest), ...
    sprintf('A_0 = %.2f \\muS', A0_rest), ...
    sprintf('B = %.2f \\muS', B_rest)}, ...
    'Units', 'normalized', 'FontSize', font_size_annot, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'Margin', 3);
legend({'\pmSEM', 'Grand average', 'P(0): A_0e^{-e\cdotx}+B'}, ...
    'FontSize', font_size_legend, 'Location', 'east', 'Box', 'off');
text(-0.12, 1.06, 'A', 'Units', 'normalized', 'FontSize', font_size_panel, 'FontWeight', 'bold');
hold off;

% ---- Panel B: Rest residuals ----
pos_1b = [margin_left, margin_bot, panel_w, h_resid];
ax1b = axes('Position', pos_1b);
hold on;
plot_residuals(rest_time, resid_rest, col_pos, col_neg, font_size_zc);
set(gca, 'FontSize', font_size, 'Box', 'on', 'TickDir', 'out');
ylabel('Resid. (m\muS)', 'FontSize', font_size_label);
xlabel('Time (s)', 'FontSize', font_size_label);
ylim([-resid_absmax_rest, resid_absmax_rest]);
text(-0.12, 1.08, 'B', 'Units', 'normalized', 'FontSize', font_size_panel, 'FontWeight', 'bold');
linkaxes([ax1a, ax1b], 'x');
xlim(ax1a, [0 180]);

print(fig1, '../../results/Fig5_EDA_Rest', '-dpng', '-r600');
fprintf('  Saved: Fig5_EDA_Rest.png (600 dpi)\n');

%% ========================================================================
%  FIGURE 2: Onset Lag Analysis (reference only — not for manuscript)
%  ========================================================================

fprintf('Creating Figure 2 (Onset Lag — reference only) ...\n');

fig2_height_in = 8.0;  % 2 stacked panels
fig2 = figure('Units', 'inches', 'Position', [4.5 0.5 fig_width_in fig2_height_in], ...
    'Color', 'w', 'PaperUnits', 'inches', ...
    'PaperSize', [fig_width_in fig2_height_in], ...
    'PaperPosition', [0 0 fig_width_in fig2_height_in]);

h_panel = 0.36;
gap2 = 0.10;
margin_bot2 = 0.10;
margin_top2 = 0.04;
start_y2 = 1 - margin_top2 - h_panel;

% ---- Panel 2A: R² vs. onset lag ----
ax2a = axes('Position', [margin_left, start_y2, panel_w, h_panel]);
hold on;
plot(tau_sweep, R2_task_sweep, '-', 'Color', col_task, 'LineWidth', 1.5);
plot(tau_sweep, R2_rest_sweep, '-', 'Color', col_rest, 'LineWidth', 1.5);
xline(tau_opt, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
plot(tau_opt, R2_task_max, 'v', 'Color', col_task, 'MarkerSize', 8, ...
    'MarkerFaceColor', col_task, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
plot(0, R2_rest_sweep(1), 'v', 'Color', col_rest, 'MarkerSize', 8, ...
    'MarkerFaceColor', col_rest, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
xlabel('Onset lag \tau (s)', 'FontSize', font_size_label);
ylabel('R^2', 'FontSize', font_size_label);
legend({sprintf('Task (peak: \\tau = %.1fs)', tau_opt), ...
        sprintf('Rest (peak: \\tau = 0s)')}, ...
    'FontSize', font_size_legend, 'Location', 'south', 'Box', 'off');
text(tau_opt + 0.4, R2_task_max + 0.003, ...
    sprintf('R^2 = %.3f', R2_task_max), ...
    'FontSize', font_size_annot, 'Color', col_task);
text(0.5, R2_rest_sweep(1) + 0.003, ...
    sprintf('R^2 = %.3f', R2_rest_sweep(1)), ...
    'FontSize', font_size_annot, 'Color', col_rest);
set(gca, 'FontSize', font_size, 'Box', 'on', 'TickDir', 'out');
xlim([0 15]);
text(-0.12, 1.06, 'A', 'Units', 'normalized', 'FontSize', font_size_panel, 'FontWeight', 'bold');
hold off;

% ---- Panel 2B: Free rate vs. e/(60-tau) ----
pos_2b = [margin_left, margin_bot2, panel_w, h_panel];
ax2b = axes('Position', pos_2b);
hold on;
plot(tau_sweep, k_free_at_tau, '-', 'Color', col_task, 'LineWidth', 1.5);
plot(tau_sweep, k_eff_at_tau, '-', 'Color', [0.2 0.5 0.8], 'LineWidth', 1.5);

% Vertical line at crossing
xline(tau_cross, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);

% Diamond marker at crossing point
plot(tau_cross, k_cross, 'd', 'Color', [0.3 0.3 0.3], ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.3 0.3], ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);

% Dashed line for e/60
yline(exp(1)/60, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);

% Annotation
text(tau_cross + 0.4, k_cross + 0.004, ...
    sprintf('k \\approx %.3f', k_cross), ...
    'FontSize', font_size_annot, 'Color', [0.3 0.3 0.3]);

xlabel('Onset lag \tau (s)', 'FontSize', font_size_label);
ylabel('k (s^{-1})', 'FontSize', font_size_label);
legend({'k_{free}(\tau)', 'e/(60-\tau)', ...
    sprintf('\\tau_{cross} = %.1fs', tau_cross), '', 'e/60'}, ...
    'FontSize', font_size_legend, 'Location', 'northeast', 'Box', 'off');
set(gca, 'FontSize', font_size, 'Box', 'on', 'TickDir', 'out');
xlim([0 15]); ylim([0 0.08]);
text(-0.12, 1.06, 'B', 'Units', 'normalized', 'FontSize', font_size_panel, 'FontWeight', 'bold');
hold off;

print(fig2, '../../results/FigSI_EDA_OnsetLag', '-dpng', '-r600');
fprintf('  Saved: FigSI_EDA_OnsetLag.png (600 dpi)\n');

%% ========================================================================
%  FIGURE 3: TASK — Optimised P(0) Fit and Residuals
%  ========================================================================

fprintf('Creating Figure 3 (Task — optimised fit) ...\n');

fig3_height_in = 8.0;  % 2 stacked panels
fig3 = figure('Units', 'inches', 'Position', [9.0 0.5 fig_width_in fig3_height_in], ...
    'Color', 'w', 'PaperUnits', 'inches', ...
    'PaperSize', [fig_width_in fig3_height_in], ...
    'PaperPosition', [0 0 fig_width_in fig3_height_in]);

resid_absmax_task = max(abs(resid_task_opt_post)) * 1.15;

% ---- Panel A: Task EDA + P(0) fit at optimal onset lag ----
ax3a = axes('Position', [margin_left, start_y, panel_w, h_main]);
hold on;

% SEM shading (full trial)
fill([task_time, fliplr(task_time)], ...
    [task_mean + task_sem, fliplr(task_mean - task_sem)], ...
    col_task, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Grand average (full trial)
plot(task_time, task_mean, '-', 'Color', col_task, 'LineWidth', line_width_data);

% P(0) fit (from tau onward)
plot(task_time(idx_opt), yfit_task_opt, '-', 'Color', col_fit, 'LineWidth', line_width_fit);

% Onset lag line
xline(tau_opt, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);

ylabel('EDA (\muS)', 'FontSize', font_size_label);
set(gca, 'FontSize', font_size, 'Box', 'on', 'TickDir', 'out');
title('Task', 'FontSize', font_size_title, 'FontWeight', 'bold');
ylim([eda_ymin, eda_ymax]);

% Annotation box
text(0.97, 0.95, {sprintf('\\tau = %.1f s', tau_opt), ...
    sprintf('R^2 = %.3f', R2_task_max), ...
    sprintf('A_0 = %.2f \\muS', A0_task_sweep(idx_task_max)), ...
    sprintf('B = %.2f \\muS', B_task_sweep(idx_task_max))}, ...
    'Units', 'normalized', 'FontSize', font_size_annot, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'Margin', 3);
legend({'\pmSEM', 'Grand average', 'P(0) fit', ...
    sprintf('\\tau = %.1f s', tau_opt)}, ...
    'FontSize', font_size_legend, 'Location', 'southeast', 'Box', 'off');
text(-0.12, 1.06, 'A', 'Units', 'normalized', 'FontSize', font_size_panel, 'FontWeight', 'bold');
hold off;

% ---- Panel B: Task residuals (optimised fit, post-onset only) ----
pos_3b = [margin_left, margin_bot, panel_w, h_resid];
ax3b = axes('Position', pos_3b);
hold on;
plot_residuals(task_time(idx_opt), resid_task_opt_post, col_pos, col_neg, font_size_zc);
xlabel('Time (s)', 'FontSize', font_size_label);
ylabel('Resid. (m\muS)', 'FontSize', font_size_label);
set(gca, 'FontSize', font_size, 'Box', 'on', 'TickDir', 'out');
ylim([-resid_absmax_task, resid_absmax_task]);
text(-0.12, 1.08, 'B', 'Units', 'normalized', 'FontSize', font_size_panel, 'FontWeight', 'bold');
linkaxes([ax3a, ax3b], 'x');
xlim(ax3a, [0 60]);

print(fig3, '../../results/Fig6_EDA_Task', '-dpng', '-r600');
fprintf('  Saved: Fig6_EDA_Task.png (600 dpi)\n');

%% Print summary for figure captions
fprintf('\n==========================================================\n');
fprintf('  FIGURE CAPTION DATA\n');
fprintf('==========================================================\n');
fprintf('\n  Figure 1 (Rest):\n');
fprintf('    Rest: N=%d, R2=%.3f, A0=%.3f uS, B=%.3f uS, RMSE=%.1f muS\n', ...
    n_rest, R2_rest, A0_rest, B_rest, RMSE_rest);

% Zero crossings for rest residuals
zc_rest = filter_zero_crossings(find_zero_crossings(rest_time, resid_rest), min_zc_gap);
fprintf('    Rest zero crossings:');
for z = 1:length(zc_rest)
    fprintf(' %.1fs', zc_rest(z));
end
fprintf('\n');

fprintf('\n  Figure 2 (Onset Lag — reference only):\n');
fprintf('    Optimal tau (R² peak): %.1f s\n', tau_opt);
fprintf('    Rate crossing tau:     %.2f s\n', tau_cross);
fprintf('    Task at tau=0:    R2=%.4f, RMSE=%.1f muS\n', R2_task_full, RMSE_task_full);
fprintf('    Task at tau=%.1f: R2=%.4f, RMSE=%.1f muS\n', ...
    tau_opt, R2_task_max, RMSE_task_sweep(idx_task_max));
fprintf('    k_free at tau=%.1f: %.4f,  e/(60-%.1f) = %.4f\n', ...
    tau_opt, k_free_at_tau(idx_task_max), tau_opt, k_eff_at_tau(idx_task_max));
fprintf('    k at crossing:     %.4f\n', k_cross);

fprintf('\n  Figure 3 (Task — optimised fit):\n');
fprintf('    Task: N=%d, tau=%.1fs, R2=%.3f, A0=%.3f uS, B=%.3f uS, RMSE=%.1f muS\n', ...
    n_task, tau_opt, R2_task_max, ...
    A0_task_sweep(idx_task_max), B_task_sweep(idx_task_max), ...
    RMSE_task_sweep(idx_task_max));

% Zero crossings for task residuals (optimised fit)
zc_task_opt = filter_zero_crossings( ...
    find_zero_crossings(task_time(idx_opt), resid_task_opt_post), min_zc_gap);
fprintf('    Task zero crossings:');
for z = 1:length(zc_task_opt)
    fprintf(' %.1fs', zc_task_opt(z));
end
fprintf('\n');

fprintf('==========================================================\n');
fprintf('Done.\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [A0, B, R2, RMSE, yfit] = fit_p0_lambda_e(x, y)
    model = @(params, x) params(1) * exp(-exp(1) * x) + params(2);
    A0_init = y(1) - y(end);
    B_init = y(end);
    p0 = [A0_init, B_init];
    lb = [0, 0]; ub = [50, 50];
    opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
        'MaxIterations', 10000, 'MaxFunctionEvaluations', 30000, ...
        'TolFun', 1e-12, 'TolX', 1e-12);
    params = lsqcurvefit(model, p0, x, y, lb, ub, opts);
    A0 = params(1); B = params(2);
    yfit = model(params, x);
    SS_res = sum((y - yfit).^2);
    SS_tot = sum((y - mean(y)).^2);
    R2 = 1 - SS_res / SS_tot;
    RMSE = sqrt(mean((y - yfit).^2)) * 1000;
end

function plot_residuals(t, resid, col_pos, col_neg, font_size_zc)
    % Plot residuals with green/red area shading and zero-crossing labels
    % Filters zero crossings with < 5s separation
    
    pos_resid = resid; pos_resid(pos_resid < 0) = 0;
    neg_resid = resid; neg_resid(neg_resid > 0) = 0;
    
    area(t, pos_resid, 'FaceColor', col_pos, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    area(t, neg_resid, 'FaceColor', col_neg, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(t, resid, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.8);
    yline(0, '-k', 'LineWidth', 0.4);
    
    % Zero crossings with minimum gap filter
    zc_all = find_zero_crossings(t, resid);
    zc = filter_zero_crossings(zc_all, 5);
    ylims = [min(resid)*1.1, max(resid)*1.1];
    for z = 1:length(zc)
        xline(zc(z), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
        text(zc(z), ylims(2)*0.85, sprintf('%.1fs', zc(z)), ...
            'FontSize', font_size_zc, 'HorizontalAlignment', 'center', ...
            'Color', [0.3 0.3 0.3]);
    end
end

function zc = find_zero_crossings(t, resid)
    signs = sign(resid);
    zc_idx = find(signs(1:end-1) .* signs(2:end) < 0);
    zc = zeros(1, length(zc_idx));
    for z = 1:length(zc_idx)
        i = zc_idx(z);
        zc(z) = t(i) - resid(i) * (t(i+1) - t(i)) / (resid(i+1) - resid(i));
    end
end

function zc_filt = filter_zero_crossings(zc, min_gap)
    % Keep only zero crossings separated by at least min_gap seconds
    if isempty(zc); zc_filt = []; return; end
    zc_filt = zc(1);
    for i = 2:length(zc)
        if zc(i) - zc_filt(end) >= min_gap
            zc_filt(end+1) = zc(i);  %#ok<AGROW>
        end
    end
end
