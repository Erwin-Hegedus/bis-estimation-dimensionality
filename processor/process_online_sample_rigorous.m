function [pred_van, pred_gre, pred_kscale, pred_loglin, pred_2d, pred_2d_fim, processor] = process_online_sample_rigorous(processor, time_k, bis_k, prop_rate_k, remi_rate_k)
    processor.sample_count = processor.sample_count + 1;
    processor.current_time = time_k;
    cfg = processor.cfg;
    
    % 1. Artifact Gating
    [y_eff, ~, processor.art, data_quality] = artifact_gate_simple(processor.art, time_k, bis_k, cfg);
    
    % 2. Initialization
    if ~processor.initialized
        processor = update_initialization_clean(processor, y_eff, time_k);
        if ~processor.initialized
            pred_van = y_eff; 
            pred_gre = y_eff;
            pred_kscale = y_eff;
            pred_loglin = y_eff;
            pred_2d = y_eff;
            pred_2d_fim = y_eff;
            return;
        end
    end
    
    % 3. Induction Detection
    if data_quality <= 2
        processor = update_induction_detection_simple(processor, y_eff, prop_rate_k, remi_rate_k);
    end
    
    % 4. PK Update
    processor.pk_state_P = update_pk_state_online(processor.pk_state_P, time_k, prop_rate_k);
    processor.pk_state_R = update_pk_state_online(processor.pk_state_R, time_k, remi_rate_k);
    Cp_P = processor.pk_state_P.Cp;
    Cp_R = processor.pk_state_R.Cp;
    
    % 5. Effect Site Update
    processor.effect_site_P = update_effect_site_with_delay_estimation(processor.effect_site_P, time_k, Cp_P, y_eff);
    processor.effect_site_R = update_effect_site_with_delay_estimation(processor.effect_site_R, time_k, Cp_R, y_eff);
    
    CeP = processor.effect_site_P.Ce_delayed_output;
    CeR = processor.effect_site_R.Ce_delayed_output;
    
    % 6. Ensure State Initialization
    if isnan(processor.E0.x)
        processor.E0.x = 93;
        processor.E0.P = 25;
    end
    if isnan(processor.BISmin.x)
        processor.BISmin.x = cfg.bismin.init;
        processor.BISmin.P = cfg.bismin.P0;
    end
    
    % 7. ESTIMATE ENDPOINTS (E0 and BISmin)
    params_van = processor.ekf_van.current_params;
    
    if data_quality <= 2
        processor = update_E0_and_BISmin(processor, params_van, CeP, CeR, y_eff, data_quality, Cp_P, Cp_R);
    end
    
    E0 = processor.E0.x;
    BISmin = processor.BISmin.x;
    
    % 8. ESTIMATE PD PARAMETERS
    learning_enabled = (data_quality <= 2);
    
    CpP = processor.pk_state_P.Cp;
    CpR = processor.pk_state_R.Cp;
    ke0P = processor.effect_site_P.ke0;
    ke0R = processor.effect_site_R.ke0;
    
    [pred_van, processor.ekf_van] = update_level_rate_rls(...
        processor.ekf_van, y_eff, CeP, CeR, E0, BISmin, ...
        'vanluchene', cfg, data_quality, learning_enabled, ...
        CpP, CpR, ke0P, ke0R);
    
    [pred_gre, processor.ekf_gre] = update_level_rate_rls(...
        processor.ekf_gre, y_eff, CeP, CeR, E0, BISmin, ...
        'greco', cfg, data_quality, learning_enabled, ...
        CpP, CpR, ke0P, ke0R);

    % 9. UPDATE 1-PARAMETER k_scale MODEL
    [pred_kscale, processor.ekf_k] = update_ekf_kscale( ...
        processor.ekf_k, y_eff, CeP, CeR, E0, cfg.BISmin_fixed, cfg, CpP, CpR);

    % 10. UPDATE 3-PARAMETER LOG-LINEAR MODEL
    [pred_loglin, processor.ekf_loglin] = update_ekf_loglin3d( ...
        processor.ekf_loglin, y_eff, CeP, CeR, E0, cfg.BISmin_fixed, cfg, CpP, CpR);

    % 11. UPDATE 2-PARAMETER (kP, kR) MODEL
    [pred_2d, processor.ekf_2d] = update_ekf_2d( ...
        processor.ekf_2d, y_eff, CeP, CeR, E0, cfg.BISmin_fixed, cfg, CpP, CpR);

    [pred_2d_fim, processor.ekf_2d_fim] = update_ekf_2d_fim(...
        processor.ekf_2d_fim, y_eff, CeP, CeR, E0, BISmin, ...
        cfg, data_quality, learning_enabled, CpP, CpR, ke0P, ke0R);

    % Personalization Tracking
    processor = track_personalization_realtime(processor, time_k, cfg);
end
