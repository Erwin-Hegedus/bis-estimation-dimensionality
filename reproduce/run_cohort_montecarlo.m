%RUN_COHORT_MONTECARLO  Cohort-wide Monte Carlo equivalence analysis.
%   Loads the canonical analysis run and reports the equivalent-set fraction
%   and normalized maximum equivalence radius across the full cohort,
%   replacing the single-patient (Patient 105) numbers in the manuscript.

addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
P = repro_paths();
review_dir = P.outdir;

% This is the one driver that needs the full published results, not just cfg.
canonical = fullfile(P.repo, 'results', 'bis_analysis_results_v6_0.mat');
if ~exist(canonical, 'file')
    error('repro:missingCanonical', ['%s not found. This driver replays the ' ...
        'per-patient fits of the published run; regenerate them with ' ...
        'run_analysis.m first.'], canonical);
end
fprintf('Loading canonical results: %s\n', canonical);
S = load(canonical, 'results', 'cfg');

t0  = tic;
coh = monte_carlo_cohort(S.results, S.cfg);
fprintf('Cohort Monte Carlo done in %.1f s.\n', toc(t0));

n_valid = numel(coh.pids);
n_skip  = sum(coh.skipped);
n_noeq  = sum(coh.frac == 0);

fprintf('\n=== Cohort Monte Carlo (N=%d valid, %d skipped, %d with no equivalent set) ===\n', ...
        n_valid, n_skip, n_noeq);
fprintf('Box diagonal (max possible L2) = %.2f\n', coh.box_diag);
fprintf('Equivalent-set fraction:  median %.1f%%   IQR [%.1f%%, %.1f%%]   range [%.1f%%, %.1f%%]\n', ...
        100 * coh.frac_stats);
fprintf('Fixed +1 BIS tolerance:   median %.1f%%   IQR [%.1f%%, %.1f%%]   range [%.1f%%, %.1f%%]\n', ...
        100 * coh.frac_fixed_stats);
fprintf('\nRadius, raw units (unit-dependent -- do not report):\n');
fprintf('  max radius / box diag:  median %.0f%%   IQR [%.0f%%, %.0f%%]   range [%.0f%%, %.0f%%]\n', ...
        100 * coh.rad_stats);
fprintf('Radius, component-normalised by each box width (unit-invariant):\n');
fprintf('  max radius:             median %.0f%%   IQR [%.0f%%, %.0f%%]   range [%.0f%%, %.0f%%]\n', ...
        100 * coh.norm_max_stats);
fprintf('  95th pct radius:        median %.0f%%   IQR [%.0f%%, %.0f%%]   range [%.0f%%, %.0f%%]\n', ...
        100 * coh.norm_p95_stats);
fprintf('Share of normalised squared distance:  C50P %.1f%%  C50R %.1f%%  gamma %.1f%%  beta %.1f%%\n', ...
        100 * coh.share_median);

% Patient 105 as the illustrative example.
i105 = find(coh.pids == 105, 1);
if ~isempty(i105)
    fprintf('Patient 105 (example):    fraction %.1f%%   normalized max radius %.0f%%\n', ...
            100 * coh.frac(i105), 100 * coh.radius_norm(i105));
end

prov = provenance_stamp(S.cfg, P.data_file, 42);
save(fullfile(review_dir, 'cohort_montecarlo.mat'), 'coh', 'prov');
fprintf('Saved: %s\n', fullfile(review_dir, 'cohort_montecarlo.mat'));

% ---- Figure: distribution across the cohort, Patient 105 marked ----------
fig = figure('Color', 'w', 'Position', [80 80 1100 460]);

subplot(1,2,1);
histogram(100 * coh.frac, 'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'w');
hold on;
if ~isempty(i105)
    xline(100 * coh.frac(i105), 'r--', 'LineWidth', 2, 'DisplayName', 'Patient 105');
    legend('Location', 'best', 'FontSize', 12);
end
xlabel('Equivalent-set fraction (%)', 'FontSize', 14);
ylabel('Number of patients', 'FontSize', 14);
title('(a) Equivalent-set fraction', 'FontSize', 15, 'FontWeight', 'bold');
set(gca, 'FontSize', 12); grid on;

subplot(1,2,2);
histogram(100 * coh.norm_p95, 'FaceColor', [0.8 0.5 0.3], 'EdgeColor', 'w');
hold on;
if ~isempty(i105)
    xline(100 * coh.norm_p95(i105), 'r--', 'LineWidth', 2, 'DisplayName', 'Patient 105');
    legend('Location', 'best', 'FontSize', 12);
end
xlabel('95th-pct radius, component-normalised (%)', 'FontSize', 14);
ylabel('Number of patients', 'FontSize', 14);
title('(b) Equivalence radius', 'FontSize', 15, 'FontWeight', 'bold');
set(gca, 'FontSize', 12); grid on;

sgtitle('Cohort-wide Monte Carlo equivalence (N=209)', 'FontSize', 16, 'FontWeight', 'bold');

fig_png = fullfile(review_dir, 'figure_cohort_montecarlo.png');
export_figure_png(fig, fig_png);
fprintf('Saved: %s\n', fig_png);
