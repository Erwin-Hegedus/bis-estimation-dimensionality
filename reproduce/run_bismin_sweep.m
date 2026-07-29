%RUN_BISMIN_SWEEP  Held-out accuracy against the fixed lower asymptote BISmin.
%   Every predictor returns BISmin + (E0 - BISmin)*(1 - H), so no model can
%   predict below BISmin and the residual is one-signed by construction wherever
%   the observed BIS falls under it. This sweep measures what that costs, and
%   whether the ranking of the models depends on where the floor is placed.
%
%   Writes reproduce/output/bismin_sweep.mat.

addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
P = repro_paths();
cfg        = P.cfg;
DATA_FILE  = P.data_file;
valid_pids = P.cohort;
review_dir = P.outdir;

q      = 1e-3;
bvals  = [0, 5, 10, 20, 30, 35, 40, 45];
N      = get_patient_count(DATA_FILE);
models = {'pop','k1d','m2d','loglin3d','van4d'};

R = struct();
for b = bvals
    fprintf('\n===== BISmin = %d =====\n', b);
    t0 = tic;
    cfg_b = cfg;
    cfg_b.BISmin_fixed = b;
    out = run_cohort_holdout(cfg_b, q, DATA_FILE, N, valid_pids);

    r = struct('BISmin', b, 'mae_in', out.mae_in, 'mae_out', out.mae_out);
    for m = models
        r.IN.(m{1})  = mean(out.mae_in.(m{1})(valid_pids),  'omitnan');
        r.OUT.(m{1}) = mean(out.mae_out.(m{1})(valid_pids), 'omitnan');
    end
    R.(sprintf('b%d', b)) = r;

    fprintf('BISmin=%2d  held-out 1D %.3f  2D %.3f  3D %.3f  4D %.3f  pop %.3f  in %.0f s\n', ...
        b, r.OUT.k1d, r.OUT.m2d, r.OUT.loglin3d, r.OUT.van4d, r.OUT.pop, toc(t0));

    valid = valid_pids;
    prov = provenance_stamp(cfg, DATA_FILE);
    save(fullfile(review_dir, 'bismin_sweep.mat'), 'R', 'bvals', 'q', 'valid', 'prov');
end

% What the floor costs at the value actually used: no predictor can return less
% than BISmin, so every observed sample below it carries a one-signed residual.
canonical = fullfile(P.repo, 'results', 'bis_analysis_results_v6_0.mat');
if exist(canonical, 'file')
    C = load(canonical, 'results', 'cfg');
    b0 = C.cfg.BISmin_fixed;
    preds = {'pred_kscale','pred_2d','pred_loglin','pred_van','pred_pop'};
    names = {'1D','2D','3D','4D','pop'};
    n_tot = 0; n_below = 0; e_all = zeros(1,5); e_below = zeros(1,5);
    for c = 1:numel(valid_pids)
        r = C.results.raw(valid_pids(c));
        if isempty(r.bis), continue; end
        bis = r.bis(:); below = bis < b0;
        n_tot = n_tot + numel(bis); n_below = n_below + sum(below);
        for k = 1:5
            p = r.(preds{k}); if isempty(p), continue; end
            e = abs(p(:) - bis); ok = ~isnan(e);
            e_all(k) = e_all(k) + sum(e(ok));
            e_below(k) = e_below(k) + sum(e(ok & below));
        end
    end
    fprintf('\n=== cost of the BISmin = %g floor (N=%d) ===\n', b0, numel(valid_pids));
    fprintf('samples below the floor: %d / %d (%.2f%%)\n', n_below, n_tot, 100*n_below/n_tot);
    for k = 1:5
        fprintf('  %-4s %.1f%% of total absolute error lies below the floor\n', ...
            names{k}, 100*e_below(k)/e_all(k));
    end
end

fprintf('\n=== held-out cohort mean MAE vs BISmin (N=%d, q=%.0e) ===\n', numel(valid_pids), q);
fprintf('%-8s %8s %8s %8s %8s %8s\n', 'BISmin', '1D', '2D', '3D', '4D', 'pop');
for b = bvals
    o = R.(sprintf('b%d', b)).OUT;
    fprintf('%-8d %8.4f %8.4f %8.4f %8.4f %8.4f\n', b, o.k1d, o.m2d, o.loglin3d, o.van4d, o.pop);
end

fprintf('\nSaved: %s\n', fullfile(review_dir, 'bismin_sweep.mat'));
