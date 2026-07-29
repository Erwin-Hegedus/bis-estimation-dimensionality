function [pred_van, pred_gre, pred_kscale, pred_loglin, pred_2d_fim, processor] = process_sample(processor, time_k, bis_k, prop_rate_k, remi_rate_k)
    processor.sample_count = processor.sample_count + 1;
    processor.current_time = time_k;
    cfg = processor.cfg;
    
    % 1. Artifact Gating
    [y_eff, Rmult, processor.art, data_quality] = artifact_gate(processor.art, time_k, bis_k, cfg);
    
    % 2. Initialization
    if ~processor.initialized
        processor = update_initialization(processor, y_eff, time_k);
        if ~processor.initialized
            pred_van = y_eff; 
            pred_gre = y_eff;
            pred_kscale = y_eff;
            pred_loglin = y_eff;
                pred_2d_fim = y_eff;
            return;
        end
    end
    
    % 4. PK Update
    processor.pk_state_P = update_pk_state_online(processor.pk_state_P, time_k, prop_rate_k);
    processor.pk_state_R = update_pk_state_online(processor.pk_state_R, time_k, remi_rate_k);
    Cp_P = processor.pk_state_P.Cp;
    Cp_R = processor.pk_state_R.Cp;
    
    % 5. Effect Site Update
    processor.effect_site_P = update_effect_site(processor.effect_site_P, time_k, Cp_P);
    processor.effect_site_R = update_effect_site(processor.effect_site_R, time_k, Cp_R);
    
    CeP = processor.effect_site_P.Ce_delayed_output;
    CeR = processor.effect_site_R.Ce_delayed_output;
    
    % 6. ENDPOINTS: population values, identical for every model
    E0 = cfg.E0_fixed;
    BISmin = cfg.BISmin_fixed;

    % 7. ESTIMATE PD PARAMETERS
    drug_effect = CeP / cfg.population_params_van(1) ...
                + CeR * 1000 / cfg.population_params_van(2);
    learning_enabled = (data_quality <= 2) && (drug_effect >= cfg.min_drug_effect);


    [pred_van, processor.ekf_van] = update_ekf_4d(...
        processor.ekf_van, y_eff, CeP, CeR, E0, BISmin, ...
        'vanluchene', cfg, learning_enabled, Rmult);

    [pred_gre, processor.ekf_gre] = update_ekf_4d(...
        processor.ekf_gre, y_eff, CeP, CeR, E0, BISmin, ...
        'greco', cfg, learning_enabled, Rmult);

    % 8. UPDATE 1-PARAMETER k_scale MODEL
    [pred_kscale, processor.ekf_k] = update_ekf_kscale( ...
        processor.ekf_k, y_eff, CeP, CeR, E0, BISmin, cfg, ...
        learning_enabled, Rmult);

    % 9. UPDATE 3-PARAMETER LOG-LINEAR MODEL
    [pred_loglin, processor.ekf_loglin] = update_ekf_loglin3d( ...
        processor.ekf_loglin, y_eff, CeP, CeR, E0, BISmin, cfg, ...
        learning_enabled, Rmult);

    % 10. UPDATE 2-PARAMETER (kP, kR) MODEL
    [pred_2d_fim, processor.ekf_2d_fim] = update_ekf_2d(...
        processor.ekf_2d_fim, y_eff, CeP, CeR, E0, BISmin, ...
        cfg, learning_enabled, Rmult);

    % Personalization Tracking
    processor = track_personalization_realtime(processor, time_k, cfg);
end
