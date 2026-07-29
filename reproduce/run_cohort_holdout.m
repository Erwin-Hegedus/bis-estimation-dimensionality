function out = run_cohort_holdout(cfg, q, DATA_FILE, N, cohort)
%RUN_COHORT_HOLDOUT  Held-out evaluation of all models over the cohort.
%   Each case is split at the midpoint of its valid window. Every model adapts
%   normally on the first half. Past the split the loop calls PREDICT_SAMPLE,
%   which takes no BIS, so nothing adapts on the held-out half.
%
%   q scales every model's process noise as Q = cfg.q * P0, which makes one
%   scalar comparable across models of different dimension.
%
%   'const' is a drug-free control: the mean BIS of the adaptation half, held
%   constant over the held-out half. It uses no infusion data and no model, so a
%   personalised predictor that does not beat it is not carrying pharmacodynamic
%   information across the split.
%
%   Only COHORT is evaluated, and every case must produce a finite result: a
%   missing one aborts rather than silently dropping a patient from one q row.
%
%   Returns per-patient in-sample and held-out MAE, and the split fraction.

    cfg.q = q;
    cohort = cohort(:).';

    names = {'pop','k1d','m2d','loglin3d','van4d','const'};
    for f = names
        mae_in.(f{1})  = nan(N,1);
        mae_out.(f{1}) = nan(N,1);
    end
    split_frac = nan(N,1);

    for i = cohort
        pid = i;
        if pid > N
            error('holdout:badCohort', 'cohort case %d exceeds the %d records in %s.', pid, N, DATA_FILE);
        end
        [bis, remi_rate, prop_rate, time] = load_patient_record(pid, DATA_FILE);
        n = min([numel(bis), numel(remi_rate), numel(prop_rate), numel(time)]);
        if n < 900
            error('holdout:shortCase', 'cohort case %d has %d samples, under the 900 required.', pid, n);
        end
        bis = bis(1:n); prop_rate_raw = prop_rate(1:n); remi_rate_raw = remi_rate(1:n); time = time(1:n);

        demo = get_patient_demo_from_mat(pid, DATA_FILE);
        [pk_prop, pk_remi] = pk_from_demographics_schnider_minto(demo, cfg);
        [prop_rate, remi_rate, bis, time, ~, n, ~] = align_bis_to_infusion( ...
            bis, prop_rate_raw, remi_rate_raw, time, 0.3/60);

        processor = init_processor(cfg, time(1), pk_prop, pk_remi, demo);

        k0 = max(1, processor.online.initialization_samples);
        T_split = k0 + floor((n - k0) / 2);
        if (n - T_split) < 50
            error('holdout:shortSplit', 'cohort case %d leaves only %d held-out samples.', pid, n - T_split);
        end

        p_pop = nan(n,1); p_1d = nan(n,1); p_2d = nan(n,1); p_3d = nan(n,1); p_4d = nan(n,1);
        for k = 1:T_split
            [p_4d(k), ~, p_1d(k), p_3d(k), p_2d(k), processor] = process_sample( ...
                processor, time(k), bis(k), prop_rate(k), remi_rate(k));
            if processor.initialized
                CeP = processor.effect_site_P.Ce_delayed_output;
                CeR = processor.effect_site_R.Ce_delayed_output;
                p_pop(k) = predict_bis_population(cfg.population_params_van, CeP, CeR, cfg.E0_fixed, 'vanluchene', cfg.BISmin_fixed);
            end
        end
        for k = (T_split+1):n
            [p_4d(k), ~, p_1d(k), p_3d(k), p_2d(k), p_pop(k), processor] = ...
                predict_sample(processor, time(k), prop_rate(k), remi_rate(k));
        end

        vi_in  = k0:T_split;
        vi_out = (T_split+1):n;
        b_const = mean(bis(vi_in), 'omitnan');
        mae_in.const(i)     = mean(abs(b_const      - bis(vi_in)),  'omitnan');
        mae_out.const(i)    = mean(abs(b_const      - bis(vi_out)), 'omitnan');
        mae_in.pop(i)       = mean(abs(p_pop(vi_in)  - bis(vi_in)),  'omitnan');
        mae_in.k1d(i)       = mean(abs(p_1d(vi_in)   - bis(vi_in)),  'omitnan');
        mae_in.m2d(i)       = mean(abs(p_2d(vi_in)   - bis(vi_in)),  'omitnan');
        mae_in.loglin3d(i)  = mean(abs(p_3d(vi_in)   - bis(vi_in)),  'omitnan');
        mae_in.van4d(i)     = mean(abs(p_4d(vi_in)   - bis(vi_in)),  'omitnan');
        mae_out.pop(i)      = mean(abs(p_pop(vi_out) - bis(vi_out)), 'omitnan');
        mae_out.k1d(i)      = mean(abs(p_1d(vi_out)  - bis(vi_out)), 'omitnan');
        mae_out.m2d(i)      = mean(abs(p_2d(vi_out)  - bis(vi_out)), 'omitnan');
        mae_out.loglin3d(i) = mean(abs(p_3d(vi_out)  - bis(vi_out)), 'omitnan');
        mae_out.van4d(i)    = mean(abs(p_4d(vi_out)  - bis(vi_out)), 'omitnan');
        split_frac(i)       = T_split / n;
    end

    for f = names
        bad = cohort(~isfinite(mae_out.(f{1})(cohort)));
        if ~isempty(bad)
            error('holdout:missingResult', '%s has no held-out MAE for case(s) %s.', ...
                f{1}, mat2str(bad));
        end
    end

    out.q = q; out.mae_in = mae_in; out.mae_out = mae_out; out.split_frac = split_frac;
end
