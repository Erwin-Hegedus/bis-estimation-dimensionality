%RUN_C50_VARIANTS  Held-out test of three candidate fixes for the C50R/C50P drift.
%   Arm 1 is the shipped estimator and must reproduce 7.716 on the 4D at q=3e-5.
%   1D/2D/3D are untouched by every arm and act as a contamination control.

here = fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(fileparts(here))));
% Prepended after the repo, so these win over ekf/init_ekf_4d and
% ekf/update_ekf_4d here. Every other driver adds the repo alone, where
% genpath orders ekf/ ahead of this folder and the shipped pair is used.
addpath(fullfile(here, 'variant_ekf'));

P = repro_paths();
cfg = P.cfg; DATA_FILE = P.data_file; valid_pids = P.cohort;
N = get_patient_count(DATA_FILE);

arms = {'shipped', 'fim_cum', 'q_norank', 'freeze_ratio'};
res = struct();

for a = 1:numel(arms)
    cfg_a = cfg;
    cfg_a.c50_variant = arms{a};
    t0 = tic;
    out = run_cohort_holdout(cfg_a, cfg.q, DATA_FILE, N, valid_pids);
    res(a).name = arms{a};
    res(a).mae_out = out.mae_out;
    res(a).mae_in  = out.mae_in;
    fprintf('%-14s held-out  1D %.4f  2D %.4f  3D %.4f  4D %.4f   (%.0f s)\n', ...
        arms{a}, mean(out.mae_out.k1d(valid_pids)), mean(out.mae_out.m2d(valid_pids)), ...
        mean(out.mae_out.loglin3d(valid_pids)), mean(out.mae_out.van4d(valid_pids)), toc(t0));
end

fprintf('\n=== held-out 4D MAE, q = %.0e, N = %d ===\n', cfg.q, numel(valid_pids));
base = res(1).mae_out.van4d(valid_pids);
fprintf('%-14s %8s %9s %10s %8s\n', 'arm', 'mean', 'median', 'd vs ship', 'p');
for a = 1:numel(arms)
    v = res(a).mae_out.van4d(valid_pids);
    k = isfinite(v) & isfinite(base);
    if a == 1
        fprintf('%-14s %8.4f %9.4f %10s %8s\n', arms{a}, mean(v(k)), median(v(k)), '-', '-');
    else
        [~, p] = ttest(v(k), base(k));
        fprintf('%-14s %8.4f %9.4f %+10.4f %8.3g\n', arms{a}, mean(v(k)), median(v(k)), ...
            mean(v(k) - base(k)), p);
    end
end

fprintf('\ncontrol - 1D/2D/3D must be identical across arms:\n');
for m = {'k1d','m2d','loglin3d'}
    d = arrayfun(@(a) mean(res(a).mae_out.(m{1})(valid_pids)), 1:numel(arms));
    fprintf('  %-9s max spread %.3e\n', m{1}, max(d) - min(d));
end

%% ratio drift actually achieved, 25 cases per arm
fprintf('\n=== C50R/C50P drift over the case (25 cases) ===\n');
sub = valid_pids(round(linspace(1, numel(valid_pids), 25)));
for a = 1:numel(arms)
    cfg_a = cfg; cfg_a.c50_variant = arms{a};
    dr = nan(numel(sub),1); dp = nan(numel(sub),1);
    for s = 1:numel(sub)
        pid = sub(s);
        [bis, remi_rate, prop_rate, time] = load_patient_record(pid, DATA_FILE);
        n = min([numel(bis), numel(remi_rate), numel(prop_rate), numel(time)]);
        bis = bis(1:n); prop_rate = prop_rate(1:n); remi_rate = remi_rate(1:n); time = time(1:n);
        demo = get_patient_demo_from_mat(pid, DATA_FILE);
        [pk_prop, pk_remi] = pk_from_demographics_schnider_minto(demo, cfg_a);
        [prop_rate, remi_rate, bis, time, ~, n, ~] = align_bis_to_infusion(bis, prop_rate, remi_rate, time, 0.3/60);
        proc = init_processor(cfg_a, time(1), pk_prop, pk_remi, demo);
        X = nan(4, n);
        for k = 1:n
            [~, ~, ~, ~, ~, proc] = process_sample(proc, time(k), bis(k), prop_rate(k), remi_rate(k));
            X(:,k) = proc.ekf_van.current_params(:);
        end
        X = X(:, all(isfinite(X),1));
        if size(X,2) < 100, continue; end
        rat = X(2,:)./X(1,:); pot = sqrt(X(1,:).*X(2,:));
        dr(s) = max(rat)/min(rat) - 1;
        dp(s) = max(pot)/min(pot) - 1;
    end
    fprintf('%-14s ratio drift median %7.4f   potency drift median %7.4f\n', ...
        arms{a}, median(dr,'omitnan'), median(dp,'omitnan'));
end

save(fullfile(P.outdir, 'c50_variants.mat'), 'res', 'arms', 'valid_pids');
fprintf('\nsaved c50_variants.mat\n');
