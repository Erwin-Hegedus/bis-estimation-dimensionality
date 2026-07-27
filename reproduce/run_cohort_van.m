function [mae_van, drift_van] = run_cohort_van(cfg, q, DATA_FILE, N)
%RUN_COHORT_VAN  4D Bouillon MAE and parameter drift over the cohort.
%   Runs only the 4D model with process noise Q = q * P0 and returns per-patient
%   MAE and drift, the parameter standard deviation over the case normalized by
%   each parameter's initial standard deviation (see PARAM_SCALES).

    S = param_scales();
    cfg.q = q;
    mae_van   = nan(N, 1);
    drift_van = nan(N, 1);
    scale = S.sd_4d;

    for i = 1:min(300, N)
        pid = i;
        try
            [bis, remi_rate, prop_rate, time] = load_patient_record(pid, DATA_FILE);
            n = min([numel(bis), numel(remi_rate), numel(prop_rate), numel(time)]);
            if n < 900
                continue;
            end
            bis = bis(1:n);
            remi_rate_raw = remi_rate(1:n);
            prop_rate_raw = prop_rate(1:n);
            time = time(1:n);

            demo = get_patient_demo_from_mat(pid, DATA_FILE);
            [pk_prop, pk_remi] = pk_from_demographics_schnider_minto(demo, cfg);
            [prop_rate, remi_rate, bis, time, ~, n, ~] = ...
                align_bis_to_infusion( ...
                bis, prop_rate_raw, remi_rate_raw, time, 0.3 / 60);

            processor = init_processor(cfg, time(1), pk_prop, pk_remi, demo);

            pred_van  = nan(n, 1);
            Xhist_van = nan(4, n);
            for k = 1:n
                [pred_van(k), ~, ~, ~, ~, processor] = process_sample( ...
                    processor, time(k), bis(k), prop_rate(k), remi_rate(k));
                if ~isempty(processor.ekf_van.current_params)
                    Xhist_van(:, k) = processor.ekf_van.current_params(:);
                end
            end

            valid_idx = max(1, processor.online.initialization_samples):n;
            if numel(valid_idx) > 10
                mae_van(i)   = mean(abs(pred_van(valid_idx) - bis(valid_idx)), 'omitnan');
                drift_van(i) = mean(std(Xhist_van(:, valid_idx), 0, 2, 'omitnan') ./ scale, 'omitnan');
            end
        catch ME
            fprintf('  patient %d skipped: %s\n', pid, ME.message);
        end
    end
end
