%% plotEEGFigures.m
% =========================================================================
% Publication Figure: Beta-Band Intra-Brain Coherence (Sampling vs Selection)
% =========================================================================
%
% Pipeline:
%   1. Load pair-level and ROI-level beta coherence CSVs
%   2. Aggregate across trials to the individual participant level
%   3. Compute paired t-tests at both ROI level (Bonferroni-corrected)
%      and channel-pair level (Benjamini-Hochberg FDR)
%   4. Generate a 3-panel figure:
%        (A) Scalp topography of mean coherence change (POST - PRE)
%        (B) ROI-level grouped bar chart with paired-difference SEM
%        (C) Distribution of t-statistics across all 351 channel pairs
%
% Panels:
%   A - Scalp map: node colour/size encodes mean delta MSC; FDR-surviving
%       pairs shown as thick red lines
%   B - Grouped bars (PRE blue, POST red) for 6 ROIs with significance
%       stars (Bonferroni-corrected in red)
%   C - Histogram of t-statistics with critical thresholds
%
% INPUT:
%   ../../data/EEG/BetaCoherence.csv
%   ../../data/EEG/BetaCoherence_ROI.csv
%
% OUTPUT:
%   ../../results/Fig_EEG.png (300 dpi)
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   March 2026
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  1. LOAD DATA
%  ========================================================================

fprintf('Loading data ...\n');

full_tbl = readtable('../../data/EEG/BetaCoherence.csv', ...
    'VariableNamingRule', 'preserve');
roi_tbl  = readtable('../../data/EEG/BetaCoherence_ROI.csv', ...
    'VariableNamingRule', 'preserve');

fprintf('  Full-pair table: %d rows\n', height(full_tbl));
fprintf('  ROI table:       %d rows\n', height(roi_tbl));

%% ========================================================================
%  2. AGGREGATE TO INDIVIDUAL LEVEL
%  ========================================================================
%  Average across trials to get one PRE and one POST value per participant
%  per channel pair (or ROI). Paired t-tests are then run at this level.

% --- Channel-pair level ---
[G_pair, G_dyad, G_pid, G_pidx, G_chA, G_chB] = findgroups( ...
    full_tbl.DyadID, full_tbl.ParticipantID, ...
    full_tbl.PairIdx, full_tbl.ChA, full_tbl.ChB);

indiv_pair = table(G_dyad, G_pid, G_pidx, G_chA, G_chB, ...
    splitapply(@mean, full_tbl.Coh_PRE,  G_pair), ...
    splitapply(@mean, full_tbl.Coh_POST, G_pair), ...
    'VariableNames', {'DyadID','ParticipantID','PairIdx','ChA','ChB', ...
                      'Coh_PRE','Coh_POST'});
indiv_pair.Delta = indiv_pair.Coh_POST - indiv_pair.Coh_PRE;

% --- ROI level ---
[G_roi, G_dyad_r, G_pid_r, G_roi_name] = findgroups( ...
    roi_tbl.DyadID, roi_tbl.ParticipantID, roi_tbl.ROI);

indiv_roi = table(G_dyad_r, G_pid_r, G_roi_name, ...
    splitapply(@mean, roi_tbl.Coh_PRE,  G_roi), ...
    splitapply(@mean, roi_tbl.Coh_POST, G_roi), ...
    'VariableNames', {'DyadID','ParticipantID','ROI','Coh_PRE','Coh_POST'});
indiv_roi.Delta = indiv_roi.Coh_POST - indiv_roi.Coh_PRE;

n_indiv = height(unique(indiv_pair(:, {'DyadID','ParticipantID'}), 'rows'));
df = n_indiv - 1;
fprintf('  N = %d participants (df = %d)\n\n', n_indiv, df);

%% ========================================================================
%  3. ROI-LEVEL STATISTICS (Bonferroni-corrected)
%  ========================================================================

roi_names  = unique(indiv_roi.ROI);
n_rois     = numel(roi_names);
bonf_alpha = 0.05 / n_rois;

roi_stats = table('Size', [n_rois 12], ...
    'VariableTypes', {'string','double','double','double','double', ...
                      'double','double','double','double','double','double','logical'}, ...
    'VariableNames', {'ROI','Mean_PRE','Mean_POST','SEM_PRE','SEM_POST', ...
                      'Mean_Delta','SD_Delta','SEM_Delta','t','p','dz','sig_bonf'});

for r = 1:n_rois
    mask = strcmp(indiv_roi.ROI, roi_names{r});
    pre  = indiv_roi.Coh_PRE(mask);
    post = indiv_roi.Coh_POST(mask);
    delt = indiv_roi.Delta(mask);
    [~, pval, ~, s] = ttest(post, pre);

    roi_stats.ROI(r)        = roi_names{r};
    roi_stats.Mean_PRE(r)   = mean(pre);
    roi_stats.Mean_POST(r)  = mean(post);
    roi_stats.SEM_PRE(r)    = std(pre)  / sqrt(n_indiv);
    roi_stats.SEM_POST(r)   = std(post) / sqrt(n_indiv);
    roi_stats.Mean_Delta(r) = mean(delt);
    roi_stats.SD_Delta(r)   = std(delt);
    roi_stats.SEM_Delta(r)  = std(delt) / sqrt(n_indiv);
    roi_stats.t(r)          = s.tstat;
    roi_stats.p(r)          = pval;
    roi_stats.dz(r)         = mean(delt) / std(delt);  % Cohen's d_z
    roi_stats.sig_bonf(r)   = pval < bonf_alpha;
end
roi_stats = sortrows(roi_stats, 'p');

%% ========================================================================
%  4. CHANNEL-PAIR LEVEL STATISTICS (FDR-corrected)
%  ========================================================================

pairs_list = sortrows(unique(indiv_pair(:, {'PairIdx','ChA','ChB'}), 'rows'), 'PairIdx');
n_pairs    = height(pairs_list);

PairIdx_out    = zeros(n_pairs, 1);
ChA_out        = strings(n_pairs, 1);
ChB_out        = strings(n_pairs, 1);
Mean_PRE_out   = zeros(n_pairs, 1);
Mean_POST_out  = zeros(n_pairs, 1);
Mean_Delta_out = zeros(n_pairs, 1);
SD_Delta_out   = zeros(n_pairs, 1);
t_out          = zeros(n_pairs, 1);
p_out          = zeros(n_pairs, 1);
d_out          = zeros(n_pairs, 1);

for pp = 1:n_pairs
    pidx = pairs_list.PairIdx(pp);
    mask = indiv_pair.PairIdx == pidx;
    pre  = indiv_pair.Coh_PRE(mask);
    post = indiv_pair.Coh_POST(mask);
    delt = indiv_pair.Delta(mask);
    [~, pval, ~, s] = ttest(post, pre);

    PairIdx_out(pp)    = pidx;
    ChA_out(pp)        = pairs_list.ChA(pp);
    ChB_out(pp)        = pairs_list.ChB(pp);
    Mean_PRE_out(pp)   = mean(pre);
    Mean_POST_out(pp)  = mean(post);
    Mean_Delta_out(pp) = mean(delt);
    SD_Delta_out(pp)   = std(delt);
    t_out(pp)          = s.tstat;
    p_out(pp)          = pval;
    d_out(pp)          = mean(delt) / std(delt);
end

pair_stats = table(PairIdx_out, ChA_out, ChB_out, ...
    Mean_PRE_out, Mean_POST_out, Mean_Delta_out, SD_Delta_out, ...
    t_out, p_out, d_out, ...
    'VariableNames', {'PairIdx','ChA','ChB','Mean_PRE','Mean_POST', ...
                      'Mean_Delta','SD_Delta','t','p','dz'});

% --- Benjamini-Hochberg FDR ---
[p_sorted, sort_idx] = sort(pair_stats.p);
bh_thresh = (1:n_pairs)' / n_pairs * 0.05;
last_sig  = find(p_sorted <= bh_thresh, 1, 'last');
sig_fdr   = false(n_pairs, 1);
if ~isempty(last_sig)
    sig_fdr(sort_idx(1:last_sig)) = true;
end
pair_stats.sig_fdr   = sig_fdr;
pair_stats.sig_unc05 = pair_stats.p < 0.05;

pair_stats = sortrows(pair_stats, 'p');

n_sig_unc05 = sum(pair_stats.sig_unc05);
n_sig_001   = sum(pair_stats.p < 0.001);
n_sig_fdr   = sum(pair_stats.sig_fdr);

fprintf('  Channel-pair statistics: %d pairs\n', n_pairs);
fprintf('    Uncorrected p < .05:  %d pairs\n', n_sig_unc05);
fprintf('    Uncorrected p < .001: %d pairs\n', n_sig_001);
fprintf('    FDR q < .05:          %d pair(s)\n\n', n_sig_fdr);

%% ========================================================================
%  5. CONSOLE SUMMARY
%  ========================================================================

fprintf('--- ROI-level results (sorted by p) ---\n');
fprintf('%-16s  PRE      POST     Delta    SEM      t(%d)    p          d_z     Bonf\n', ...
    'ROI', df);
fprintf('%s\n', repmat('-', 1, 100));
for r = 1:n_rois
    bonf_str = '';
    if roi_stats.sig_bonf(r), bonf_str = '  *'; end
    fprintf('%-16s  %.4f   %.4f   %+.4f   %.4f   %+6.3f   %-10s  %+.3f%s\n', ...
        roi_stats.ROI(r), ...
        roi_stats.Mean_PRE(r), roi_stats.Mean_POST(r), ...
        roi_stats.Mean_Delta(r), roi_stats.SEM_Delta(r), ...
        roi_stats.t(r), format_p(roi_stats.p(r)), ...
        roi_stats.dz(r), bonf_str);
end
fprintf('\n  Bonferroni threshold: alpha = %.4f (0.05 / %d ROIs)\n\n', bonf_alpha, n_rois);

fprintf('--- Top 10 channel pairs (by p-value) ---\n');
fprintf('%-10s  PRE      POST     Delta    t(%d)    p          d_z     FDR\n', 'Pair', df);
fprintf('%s\n', repmat('-', 1, 85));
n_show = min(10, height(pair_stats));
for pp = 1:n_show
    fdr_str = '';
    if pair_stats.sig_fdr(pp), fdr_str = '  *'; end
    fprintf('%-10s  %.4f   %.4f   %+.4f   %+6.3f   %-10s  %+.3f%s\n', ...
        sprintf('%s-%s', pair_stats.ChA(pp), pair_stats.ChB(pp)), ...
        pair_stats.Mean_PRE(pp), pair_stats.Mean_POST(pp), ...
        pair_stats.Mean_Delta(pp), ...
        pair_stats.t(pp), format_p(pair_stats.p(pp)), ...
        pair_stats.dz(pp), fdr_str);
end
fprintf('\n');

%% ========================================================================
%  6. STYLING PARAMETERS
%  ========================================================================

% 27 EEG channels and their standard 10-20 polar projection coordinates
ch_labels = {'Fp1','Fp2','F7','F3','Fz','F4','F8', ...
             'FC5','FC1','FC2','FC6', ...
             'T7','C3','Cz','C4','T8', ...
             'CP5','CP1','CP2','CP6', ...
             'P7','P3','Pz','P4','P8','O1','O2'};
coords = get_10_20_coords(ch_labels);
n_ch   = numel(ch_labels);

col_pre  = [0.35 0.58 0.85];    % Blue  - Sampling epoch
col_post = [0.88 0.28 0.20];    % Red   - Selection epoch
col_fdr  = [0.80 0.10 0.10];    % Deep red - FDR lines
col_grey = [0.55 0.55 0.55];    % Grey  - histogram bars

font_sz       = 12;
font_sz_label = 13;
font_sz_title = 14;
font_sz_annot = 11;
font_sz_panel = 18;

%% ========================================================================
%  7. CREATE FIGURE - 3 vertically stacked panels
%  ========================================================================

fprintf('Creating figure ...\n');

fig_w = 7.5;   fig_h = 14.0;
fig = figure('Units', 'inches', 'Position', [0.5 0.5 fig_w fig_h], ...
    'Color', 'w', 'PaperUnits', 'inches', ...
    'PaperSize', [fig_w fig_h], 'PaperPosition', [0 0 fig_w fig_h]);

% Layout geometry (normalised)
ml  = 0.13;  mr = 0.05;  mt = 0.03;  mb = 0.05;  gap = 0.06;
pw  = 1 - ml - mr;

h_a = 0.32;   h_b = 0.24;   h_c = 0.22;
y_a = 1 - mt - h_a;
y_b = y_a - gap - h_b;
y_c = y_b - gap - h_c;

%% ---- Panel A: Scalp topography ----------------------------------------

ax_a = axes('Position', [ml y_a pw h_a]);
draw_scalp(ax_a);

% FDR-significant pairs: thick coloured lines
fdr_pairs = pair_stats(pair_stats.sig_fdr, :);
for fp = 1:height(fdr_pairs)
    ia = find(strcmp(ch_labels, fdr_pairs.ChA{fp}));
    ib = find(strcmp(ch_labels, fdr_pairs.ChB{fp}));
    if ~isempty(ia) && ~isempty(ib)
        plot(ax_a, coords([ia ib], 1)', coords([ia ib], 2)', '-', ...
            'Color', col_fdr, 'LineWidth', 4.5);
    end
end

% Node-level mean delta (averaged across all pairs involving each node)
node_delta = zeros(n_ch, 1);
node_count = zeros(n_ch, 1);
for pp = 1:height(pair_stats)
    ia = find(strcmp(ch_labels, pair_stats.ChA{pp}));
    ib = find(strcmp(ch_labels, pair_stats.ChB{pp}));
    if ~isempty(ia)
        node_delta(ia) = node_delta(ia) + pair_stats.Mean_Delta(pp);
        node_count(ia) = node_count(ia) + 1;
    end
    if ~isempty(ib)
        node_delta(ib) = node_delta(ib) + pair_stats.Mean_Delta(pp);
        node_count(ib) = node_count(ib) + 1;
    end
end
node_delta = node_delta ./ max(node_count, 1);

% Diverging colormap (blue-white-red)
nc = 256;
cmap_div = [linspace(0.10, 1, nc/2)', linspace(0.35, 1, nc/2)', linspace(0.78, 1, nc/2)'; ...
            linspace(1, 0.85, nc/2)', linspace(1, 0.12, nc/2)', linspace(1, 0.10, nc/2)'];
clim_v = max(abs(node_delta)) * 1.05;
if clim_v == 0, clim_v = 1; end
norm_d  = (node_delta + clim_v) / (2 * clim_v);
node_sz = 45 + 220 * (abs(node_delta) ./ (max(abs(node_delta)) + eps));

for c = 1:n_ch
    ci = max(1, min(nc, round(norm_d(c) * nc)));
    scatter(ax_a, coords(c, 1), coords(c, 2), node_sz(c), cmap_div(ci, :), ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.6);
end
draw_electrodes(ax_a, coords, ch_labels);

% Compact colorbar
cb_x = [1.08 1.13];
cb_y = linspace(-0.75, 0.75, nc);
for k = 1:nc-1
    patch(ax_a, [cb_x(1) cb_x(2) cb_x(2) cb_x(1)], ...
        [cb_y(k) cb_y(k) cb_y(k+1) cb_y(k+1)], ...
        cmap_div(k, :), 'EdgeColor', 'none');
end
text(ax_a, 1.17, 0.78, '+\DeltaMSC', 'FontSize', font_sz_annot, 'FontWeight', 'bold');
text(ax_a, 1.17, 0.00, '0', 'FontSize', font_sz_annot);
text(ax_a, 1.17, -0.78, sprintf('%.4f', -clim_v), 'FontSize', font_sz_annot);

% FDR annotation
if n_sig_fdr > 0
    fdr_label_parts = cell(n_sig_fdr, 1);
    for fp = 1:n_sig_fdr
        fdr_label_parts{fp} = sprintf('%s-%s', ...
            fdr_pairs.ChA{fp}, fdr_pairs.ChB{fp});
    end
    fdr_label = sprintf('FDR q<.05: %s (thick line)', ...
        strjoin(fdr_label_parts, ', '));
else
    fdr_label = 'FDR q<.05: none';
end
text(ax_a, 0, -1.35, fdr_label, ...
    'FontSize', font_sz_annot, 'Color', col_fdr, 'HorizontalAlignment', 'center');

title(ax_a, 'Beta Coherence Change (POST - PRE)', ...
    'FontSize', font_sz_title, 'FontWeight', 'bold');
text(-0.10, 1.06, 'A', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ---- Panel B: ROI grouped bars ----------------------------------------

ax_b = axes('Position', [ml y_b pw h_b]);

% Display order and labels for ROIs
roi_display_order  = {'Parietal','FrontoParietal','Left','Central','Frontal','Right'};
roi_display_labels = {'Parietal','F-P','Left','Central','Frontal','Right'};
n_r = numel(roi_display_order);

% Reorder roi_stats to match display order
roi_order_idx = zeros(n_r, 1);
for r = 1:n_r
    idx = find(strcmp(roi_stats.ROI, roi_display_order{r}));
    if isempty(idx)
        error('ROI "%s" not found in roi_stats.', roi_display_order{r});
    end
    roi_order_idx(r) = idx;
end
roi_plot = roi_stats(roi_order_idx, :);

b = bar(ax_b, 1:n_r, [roi_plot.Mean_PRE roi_plot.Mean_POST], 'grouped');
b(1).FaceColor = col_pre;  b(1).FaceAlpha = 0.85;
b(2).FaceColor = col_post; b(2).FaceAlpha = 0.85;
hold(ax_b, 'on');

% Error bars: SEM of the paired difference, centred between bars
x_pre  = (1:n_r) + b(1).XOffset;
x_post = (1:n_r) + b(2).XOffset;
x_mid  = (x_pre + x_post) / 2;
y_top  = max(roi_plot.Mean_PRE, roi_plot.Mean_POST);
errorbar(ax_b, x_mid, y_top, zeros(n_r, 1), roi_plot.SEM_Delta, ...
    'Color', 'k', 'LineStyle', 'none', 'LineWidth', 1.4, 'CapSize', 8);

% Significance stars
for r = 1:n_r
    stars = pstars(roi_plot.p(r));
    if isempty(stars), continue; end
    col_s = [0 0 0];
    if roi_plot.sig_bonf(r), col_s = [0.75 0 0]; end
    text(ax_b, x_mid(r), ...
        y_top(r) + roi_plot.SEM_Delta(r) + 0.002, ...
        stars, 'HorizontalAlignment', 'center', ...
        'FontSize', 14, 'Color', col_s, 'FontWeight', 'bold');
end

set(ax_b, 'XTick', 1:n_r, 'XTickLabel', roi_display_labels, ...
    'FontSize', font_sz, 'Box', 'on', 'TickDir', 'out');
ylabel(ax_b, 'Mean MSC (beta 13-30 Hz)', 'FontSize', font_sz_label);
title(ax_b, 'ROI-Level Coherence', ...
    'FontSize', font_sz_title, 'FontWeight', 'bold');
legend(ax_b, {'Sampling [17-21 s]', 'Selection [34-38 s]'}, ...
    'Location', 'northwest', 'Box', 'off', 'FontSize', font_sz_annot);
ylim(ax_b, [0.27 0.42]);

% Footnote
text(ax_b, 0.99, 0.98, ...
    {'* p<.05, ** p<.01  (red = Bonf.)', 'error bars = SEM of paired diff.'}, ...
    'Units', 'normalized', 'FontSize', 9, ...
    'Color', [0.3 0.3 0.3], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

text(-0.10, 1.08, 'B', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ---- Panel C: t-statistic distribution --------------------------------

ax_c = axes('Position', [ml y_c pw h_c]);
histogram(ax_c, pair_stats.t, 60, ...
    'FaceColor', col_grey, 'EdgeColor', 'w', 'FaceAlpha', 0.85);
hold(ax_c, 'on');

t_crit_05  = tinv(0.975,  df);
t_crit_001 = tinv(0.9995, df);
yl = ylim(ax_c);

% Shaded region for p < .05
patch(ax_c, [t_crit_05 7 7 t_crit_05], [0 0 yl(2) yl(2)], ...
    [0.90 0.68 0.52], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
xline(ax_c, t_crit_05,  '--', 'Color', [0.82 0.42 0.08], ...
    'LineWidth', 1.5, 'Label', sprintf('p<.05 (%d)', n_sig_unc05), ...
    'FontSize', font_sz_annot);
xline(ax_c, t_crit_001, ':', 'Color', [0.72 0.08 0.08], ...
    'LineWidth', 1.5, 'Label', sprintf('p<.001 (%d)', n_sig_001), ...
    'FontSize', font_sz_annot);
xline(ax_c, 0, 'k-', 'LineWidth', 0.8);

% FDR threshold
if n_sig_fdr > 0
    fdr_pair_names = cell(n_sig_fdr, 1);
    for fp = 1:n_sig_fdr
        fdr_pair_names{fp} = sprintf('%s-%s', ...
            fdr_pairs.ChA{fp}, fdr_pairs.ChB{fp});
    end
    t_fdr_thresh = min(pair_stats.t(pair_stats.sig_fdr));
    xline(ax_c, t_fdr_thresh, '-', 'Color', [0 0.55 0], ...
        'LineWidth', 2, 'Label', sprintf('FDR q<.05 (%d: %s)', ...
        n_sig_fdr, strjoin(fdr_pair_names, ', ')), ...
        'FontSize', font_sz_annot);
end

xlim(ax_c, [-5 7]);
xlabel(ax_c, 't-statistic (Selection - Sampling)', 'FontSize', font_sz_label);
ylabel(ax_c, 'Number of pairs', 'FontSize', font_sz_label);
title(ax_c, sprintf('All %d Channel Pairs  (N = %d)', n_pairs, n_indiv), ...
    'FontSize', font_sz_title, 'FontWeight', 'bold');
set(ax_c, 'FontSize', font_sz, 'Box', 'on', 'TickDir', 'out');

text(-0.10, 1.08, 'C', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ========================================================================
%  8. SAVE
%  ========================================================================

print(fig, '../../results/Fig_EEG', '-dpng', '-r300');
fprintf('  Saved: Fig_EEG.png (300 dpi)\n');

%% ========================================================================
%  9. PRINT FULL SUMMARY
%  ========================================================================

fprintf('\n==========================================================\n');
fprintf('  FIGURE SUMMARY\n');
fprintf('==========================================================\n');
fprintf('  Participants:         N = %d (df = %d)\n', n_indiv, df);
fprintf('  Channel pairs:        %d\n', n_pairs);
fprintf('  ROIs:                 %d (Bonferroni alpha = %.4f)\n', n_rois, bonf_alpha);
fprintf('  --\n');
fprintf('  Uncorrected p < .05:  %d pairs\n', n_sig_unc05);
fprintf('  Uncorrected p < .001: %d pairs\n', n_sig_001);
fprintf('  FDR q < .05:          %d pair(s)\n', n_sig_fdr);
fprintf('  --\n');
fprintf('  Output: Fig_EEG.png (300 dpi)\n');
fprintf('==========================================================\n');
fprintf('Done.\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function s = pstars(p)
% PSTARS  Return significance stars for a p-value.
    if     p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = '';
    end
end

function s = format_p(p)
% FORMAT_P  Format a p-value string for console display.
    if p < 0.001
        s = sprintf('%.2e', p);
    else
        s = sprintf('%.4f', p);
    end
end

function coords = get_10_20_coords(ch_labels)
% GET_10_20_COORDS  Return [x y] for standard 10-20 electrode positions.
    locs = struct( ...
        'Fp1',[-0.31  0.95], 'Fp2',[ 0.31  0.95], ...
        'F7', [-0.81  0.59], 'F3', [-0.41  0.59], 'Fz', [ 0.00  0.59], ...
        'F4', [ 0.41  0.59], 'F8', [ 0.81  0.59], ...
        'FC5',[-0.61  0.37], 'FC1',[-0.22  0.37], 'FC2',[ 0.22  0.37], ...
        'FC6',[ 0.61  0.37], ...
        'T7', [-1.00  0.00], 'C3', [-0.50  0.00], 'Cz', [ 0.00  0.00], ...
        'C4', [ 0.50  0.00], 'T8', [ 1.00  0.00], ...
        'CP5',[-0.61 -0.37], 'CP1',[-0.22 -0.37], 'CP2',[ 0.22 -0.37], ...
        'CP6',[ 0.61 -0.37], ...
        'P7', [-0.81 -0.59], 'P3', [-0.41 -0.59], 'Pz', [ 0.00 -0.59], ...
        'P4', [ 0.41 -0.59], 'P8', [ 0.81 -0.59], ...
        'O1', [-0.31 -0.95], 'O2', [ 0.31 -0.95]);
    coords = zeros(numel(ch_labels), 2);
    for c = 1:numel(ch_labels)
        coords(c,:) = locs.(ch_labels{c});
    end
end

function draw_scalp(ax)
% DRAW_SCALP  Draw scalp outline with nose and ears.
    hold(ax, 'on');
    theta = linspace(0, 2*pi, 100);
    plot(ax, cos(theta), sin(theta), 'k-', 'LineWidth', 1.5);
    % Nose
    plot(ax, [-0.1 0 0.1], [1 1.15 1], 'k-', 'LineWidth', 1.2);
    % Left ear
    plot(ax, [-1.05 -1.12 -1.12 -1.05], [-0.15 -0.08 0.08 0.15], ...
        'k-', 'LineWidth', 1);
    % Right ear
    plot(ax, [1.05 1.12 1.12 1.05], [-0.15 -0.08 0.08 0.15], ...
        'k-', 'LineWidth', 1);
    axis(ax, 'equal');
    xlim(ax, [-1.45 1.45]); ylim(ax, [-1.45 1.30]);
    set(ax, 'XTick', [], 'YTick', []); box(ax, 'off');
end

function draw_electrodes(ax, coords, ch_labels)
% DRAW_ELECTRODES  Overlay electrode markers and labels.
    scatter(ax, coords(:,1), coords(:,2), 28, [0.3 0.3 0.3], ...
        'filled', 'MarkerFaceAlpha', 0.6);
    for c = 1:numel(ch_labels)
        text(ax, coords(c,1), coords(c,2) + 0.11, ch_labels{c}, ...
            'FontSize', 9, 'HorizontalAlignment', 'center', ...
            'Color', [0.15 0.15 0.15]);
    end
end
