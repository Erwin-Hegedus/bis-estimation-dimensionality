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
    out = run_cohort_holdout(cfg_b, q, DATA_FILE, N);

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

fprintf('\n=== held-out cohort mean MAE vs BISmin (N=%d, q=%.0e) ===\n', numel(valid_pids), q);
fprintf('%-8s %8s %8s %8s %8s %8s\n', 'BISmin', '1D', '2D', '3D', '4D', 'pop');
for b = bvals
    o = R.(sprintf('b%d', b)).OUT;
    fprintf('%-8d %8.4f %8.4f %8.4f %8.4f %8.4f\n', b, o.k1d, o.m2d, o.loglin3d, o.van4d, o.pop);
end

fprintf('\nSaved: %s\n', fullfile(review_dir, 'bismin_sweep.mat'));
