%RUN_KGAMMA_EXPERIMENT  Evaluate the FIM-motivated shared-potency + gamma model.
%   Runs the full N=209 cohort at the canonical operating point and the same
%   six-point frozen-half holdout sweep used in the revised manuscript. The
%   script saves one provenance-stamped experimental artifact and does not
%   replace any canonical paper result.

addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
P = repro_paths();
cfg = P.cfg;
N = get_patient_count(P.data_file);
cohort = P.cohort(:);
q_grid = [0, 1e-5, 3e-5, 1e-4, 1e-3, 1e-2];

fprintf('=== k-gamma experiment: full-case run at q = %.0e ===\n', cfg.q);
t0 = tic;
whole = run_cohort_allmodels(cfg, cfg.q, P.data_file, N, [], cohort);
fprintf('Full-case run completed in %.1f s.\n', toc(t0));

models = {'pop','k1d','kgamma2d','m2d','loglin3d','van4d','const'};
nq = numel(q_grid);
holdout_runs = cell(nq,1);
for m = models
    holdout_mean.(m{1}) = nan(1,nq);
end

for j = 1:nq
    fprintf('\n=== k-gamma holdout %d/%d: q = %.0e ===\n', j, nq, q_grid(j));
    t0 = tic;
    holdout_runs{j} = run_cohort_holdout(cfg, q_grid(j), P.data_file, N, cohort);
    for m = models
        holdout_mean.(m{1})(j) = mean(holdout_runs{j}.mae_out.(m{1})(cohort), 'omitnan');
    end
    fprintf('OUT: k-gamma %.3f | 1D %.3f | old 2D %.3f | 3D %.3f | 4D %.3f | const %.3f (%.1f s)\n', ...
        holdout_mean.kgamma2d(j), holdout_mean.k1d(j), holdout_mean.m2d(j), ...
        holdout_mean.loglin3d(j), holdout_mean.van4d(j), holdout_mean.const(j), toc(t0));
end

whole_mean = struct();
for m = {'pop','k1d','kgamma2d','m2d','loglin3d','van4d'}
    whole_mean.(m{1}) = mean(whole.mae.(m{1})(cohort), 'omitnan');
end

kg = whole.mae.kgamma2d(cohort);
comparison = struct();
for m = {'k1d','m2d','loglin3d','van4d'}
    other = whole.mae.(m{1})(cohort);
    d = kg - other;
    comparison.(m{1}) = struct('mean_difference', mean(d, 'omitnan'), ...
        'median_difference', median(d, 'omitnan'), ...
        'kgamma_better', sum(d < 0), 'ties', sum(d == 0), ...
        'kgamma_worse', sum(d > 0), 'equivalent_within_1', sum(abs(d) <= 1));
end

diag2d = whole.kgamma_diag;
diag2d.bound_cases = sum(any(diag2d.at_bound(:,cohort), 1));
diag2d.final_param_median = median(diag2d.final_params(:,cohort), 2, 'omitnan');
diag2d.final_param_iqr = [prctile(diag2d.final_params(:,cohort), 25, 2), ...
                          prctile(diag2d.final_params(:,cohort), 75, 2)];
diag2d.relative_sd_median = median(sqrt(diag2d.final_P_diag(:,cohort)), 2, 'omitnan');
diag2d.lambda2_share_median = median(diag2d.fim_eigenvalues(2,cohort) ./ ...
    max(sum(diag2d.fim_eigenvalues(:,cohort), 1), eps), 'omitnan');
diag2d.lambda2_cum_share_median = median(diag2d.fim_cum_eigenvalues(2,cohort) ./ ...
    max(sum(diag2d.fim_cum_eigenvalues(:,cohort), 1), eps), 'omitnan');
diag2d.rate_limit_fraction = sum(diag2d.rate_limit_count(cohort)) / ...
    max(sum(diag2d.update_count(cohort)), 1);

experiment = struct('q_grid', q_grid, 'cohort', cohort, 'whole', whole, ...
    'whole_mean', whole_mean, 'comparison', comparison, ...
    'holdout_mean', holdout_mean, 'holdout_runs', {holdout_runs}, ...
    'kgamma_diag', diag2d);
prov = provenance_stamp(cfg, P.data_file);
prov.experiment_note = ['k-gamma model as published; the canonical run in ' ...
    'results/bis_analysis_results_v6_0.mat carries this estimator.'];
save(fullfile(P.outdir, 'kgamma_experiment.mat'), 'experiment', 'prov', '-v7.3');

fprintf('\n=== k-gamma summary ===\n');
fprintf('Whole-case mean MAE: %.3f\n', whole_mean.kgamma2d);
fprintf('Final k median [IQR]: %.3f [%.3f, %.3f]\n', diag2d.final_param_median(1), ...
    diag2d.final_param_iqr(1,1), diag2d.final_param_iqr(1,2));
fprintf('Final gamma median [IQR]: %.3f [%.3f, %.3f]\n', diag2d.final_param_median(2), ...
    diag2d.final_param_iqr(2,1), diag2d.final_param_iqr(2,2));
fprintf('Cases at a bound: %d/%d\n', diag2d.bound_cases, numel(cohort));
fprintf('Saved %s\n', fullfile(P.outdir, 'kgamma_experiment.mat'));
