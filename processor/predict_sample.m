function [pred_van, pred_gre, pred_kscale, pred_loglin, pred_2d_fim, pred_pop, processor, pred_kgamma] = ...
    predict_sample(processor, time_k, prop_rate_k, remi_rate_k)
%PREDICT_SAMPLE  Advance the PK and effect-site states from the infusions and
%   predict with the parameters each model currently holds. Takes no BIS, so
%   nothing adapts.

    processor.sample_count = processor.sample_count + 1;
    processor.current_time = time_k;
    cfg = processor.cfg;

    processor.pk_state_P = update_pk_state_online(processor.pk_state_P, time_k, prop_rate_k);
    processor.pk_state_R = update_pk_state_online(processor.pk_state_R, time_k, remi_rate_k);

    processor.effect_site_P = update_effect_site(processor.effect_site_P, time_k, processor.pk_state_P.Cp);
    processor.effect_site_R = update_effect_site(processor.effect_site_R, time_k, processor.pk_state_R.Cp);

    CeP = processor.effect_site_P.Ce_delayed_output;
    CeR = processor.effect_site_R.Ce_delayed_output;

    E0     = cfg.E0_fixed;
    BISmin = cfg.BISmin_fixed;

    pred_van    = predict_bis_4d(processor.ekf_van.current_params, CeP, CeR, E0, BISmin, 'vanluchene');
    pred_gre    = predict_bis_4d(processor.ekf_gre.current_params, CeP, CeR, E0, BISmin, 'greco');
    pred_kscale = predict_bis_1d_internal(processor.ekf_k.k, CeP, CeR, E0, BISmin, processor.ekf_k);
    pred_kgamma = predict_bis_kgamma_model(processor.ekf_kgamma.current_params(1), ...
        processor.ekf_kgamma.current_params(2), CeP, CeR, E0, BISmin, cfg);
    pred_loglin = predict_bis_loglin3d(processor.ekf_loglin.current_params, CeP, CeR, E0, BISmin, processor.ekf_loglin);
    pred_2d_fim = predict_bis_2d_model(processor.ekf_2d_fim.current_params(1), ...
        processor.ekf_2d_fim.current_params(2), CeP, CeR, E0, BISmin, cfg);
    pred_pop    = predict_bis_population(cfg.population_params_van, CeP, CeR, ...
        cfg.E0_fixed, 'vanluchene', cfg.BISmin_fixed);
end
