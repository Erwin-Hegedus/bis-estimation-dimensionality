function out = run_cohort_allmodels(cfg, q, DATA_FILE, N, save_pids, cohort)
%RUN_COHORT_ALLMODELS  In-sample cohort MAE and parameter roughness vs q.
%   Runs all five models over the cohort with every process noise set to
%   Q = q * S (see PARAM_SCALES) and returns per-patient MAE together with a
%   roughness measure: the mean absolute sample-to-sample parameter change,
%   normalized by each parameter's initial standard deviation. Roughness
%   separates a model that fits by tracking the signal from one that fits by
%   tracking noise.
%
%   Trajectories are retained for the patient indices listed in save_pids.

    S = param_scales();
    cfg.q = q;

    if nargin < 6 || isempty(cohort), cohort = 1:min(300, N); end
    cohort = cohort(:).';

    names = {'pop','k1d','kgamma2d','m2d','loglin3d','van4d'};
    mae   = struct(); rough = struct();
    for f = names, mae.(f{1}) = nan(N,1); rough.(f{1}) = nan(N,1); end
    traj = struct();

    kgamma_diag = struct('final_params', nan(2,N), 'final_P', nan(2,2,N), ...
        'final_P_diag', nan(2,N), 'min_P_eigenvalue', nan(N,1), ...
        'fim_eigenvalues', nan(2,N), 'fim_cum_eigenvalues', nan(2,N), ...
        'at_bound', false(2,N), 'update_count', nan(N,1), ...
        'projection_count', nan(N,1), 'rate_limit_count', nan(N,1));

    for i = cohort
        pid = i;
        try
            [bis, remi_rate, prop_rate, time] = load_patient_record(pid, DATA_FILE);
            n = min([numel(bis), numel(remi_rate), numel(prop_rate), numel(time)]);
            if n < 900, continue; end
            bis = bis(1:n); prop_rate_raw = prop_rate(1:n); remi_rate_raw = remi_rate(1:n); time = time(1:n);

            demo = get_patient_demo_from_mat(pid, DATA_FILE);
            [pk_prop, pk_remi] = pk_from_demographics_schnider_minto(demo, cfg);
            [prop_rate, remi_rate, bis, time, ~, n, ~] = align_bis_to_infusion( ...
                bis, prop_rate_raw, remi_rate_raw, time, 0.3/60);

            processor = init_processor(cfg, time(1), pk_prop, pk_remi, demo);

            p_pop = nan(n,1); p_1d = nan(n,1); p_kg = nan(n,1); p_2d = nan(n,1); p_3d = nan(n,1); p_4d = nan(n,1);
            th1 = nan(1,n); thkg = nan(2,n); th2 = nan(2,n); th3 = nan(3,n); th4 = nan(4,n);
            for k = 1:n
                [p_4d(k), ~, p_1d(k), p_3d(k), p_2d(k), processor, p_kg(k)] = process_sample( ...
                    processor, time(k), bis(k), prop_rate(k), remi_rate(k));
                if processor.initialized
                    CeP = processor.effect_site_P.Ce_delayed_output;
                    CeR = processor.effect_site_R.Ce_delayed_output;
                    p_pop(k) = predict_bis_population(cfg.population_params_van, CeP, CeR, cfg.E0_fixed, 'vanluchene', cfg.BISmin_fixed);
                    th1(k)   = processor.ekf_k.k;
                    thkg(:,k) = processor.ekf_kgamma.current_params(:);
                    th2(:,k) = processor.ekf_2d_fim.current_params(:);
                    th3(:,k) = processor.ekf_loglin.current_params(:);
                    th4(:,k) = processor.ekf_van.current_params(:);
                end
            end

            vi = max(1, processor.online.initialization_samples):n;
            if numel(vi) <= 10, continue; end
            mae.pop(i)      = mean(abs(p_pop(vi) - bis(vi)), 'omitnan');
            mae.k1d(i)      = mean(abs(p_1d(vi)  - bis(vi)), 'omitnan');
            mae.kgamma2d(i) = mean(abs(p_kg(vi)  - bis(vi)), 'omitnan');
            mae.m2d(i)      = mean(abs(p_2d(vi)  - bis(vi)), 'omitnan');
            mae.loglin3d(i) = mean(abs(p_3d(vi)  - bis(vi)), 'omitnan');
            mae.van4d(i)    = mean(abs(p_4d(vi)  - bis(vi)), 'omitnan');

            rough.k1d(i)      = tv(th1(:,vi), S.sd_1d);
            rough.kgamma2d(i) = tv(thkg(:,vi), S.sd_kgamma);
            rough.m2d(i)      = tv(th2(:,vi), S.sd_2d);
            rough.loglin3d(i) = tv(th3(:,vi), S.sd_3d);
            rough.van4d(i)    = tv(th4(:,vi), S.sd_4d);

            if any(save_pids == pid)
                tag = sprintf('p%d', pid);
                traj.(tag) = struct('bis', bis(vi), 'th1', th1(:,vi), ...
                    'thkg', thkg(:,vi), 'th2', th2(:,vi), 'th3', th3(:,vi), 'th4', th4(:,vi));
            end

            kgamma_diag.final_params(:,i) = processor.ekf_kgamma.current_params(:);
            kgamma_diag.final_P(:,:,i) = processor.ekf_kgamma.P;
            kgamma_diag.final_P_diag(:,i) = diag(processor.ekf_kgamma.P);
            kgamma_diag.min_P_eigenvalue(i) = min(eig(processor.ekf_kgamma.P));
            kgamma_diag.fim_eigenvalues(:,i) = processor.ekf_kgamma.FIM_eigenvalues(:);
            kgamma_diag.fim_cum_eigenvalues(:,i) = sort(real(eig( ...
                processor.ekf_kgamma.FIM_cum)), 'descend');
            kgamma_diag.at_bound(:,i) = ...
                abs(processor.ekf_kgamma.current_params(:) - processor.ekf_kgamma.lb(:)) < 1e-9 | ...
                abs(processor.ekf_kgamma.current_params(:) - processor.ekf_kgamma.ub(:)) < 1e-9;
            kgamma_diag.update_count(i) = processor.ekf_kgamma.update_count;
            kgamma_diag.projection_count(i) = processor.ekf_kgamma.projection_count;
            kgamma_diag.rate_limit_count(i) = processor.ekf_kgamma.rate_limit_count;
        catch ME
            fprintf('  patient %d skipped: %s\n', pid, ME.message);
        end
    end

    out.q = q; out.mae = mae; out.rough = rough; out.traj = traj;
    out.cohort = cohort; out.kgamma_diag = kgamma_diag;
end

function r = tv(Th, sc)
    r = mean(mean(abs(diff(Th, 1, 2)), 2, 'omitnan') ./ sc, 'omitnan');
end
