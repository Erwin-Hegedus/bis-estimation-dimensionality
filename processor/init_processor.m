function processor = init_processor(cfg, start_time, pk_prop, pk_remi, ~)
    processor = struct();
    processor.cfg = cfg;
    processor.sample_count = 0;
    processor.start_time = start_time;
    processor.current_time = start_time;
    processor.initialized = false;
    processor.art = init_artifact_gate(cfg);
    processor.fast_init = struct('bis_buffer',[],'samples_needed',10,'ready',false);
    processor.ekf_k = init_ekf_kscale(cfg);
    processor.ekf_loglin = init_ekf_loglin3d(cfg);
    processor.ekf_2d_fim = init_ekf_2d(cfg);
    
    processor.personalization.van = struct('detected',false(4,1), 'time',nan(4,1));
    processor.personalization.gre = struct('detected',false(4,1), 'time',nan(4,1));
    processor.personalization.baseline_van = [];
    processor.personalization.baseline_gre = [];
    
    if nargin < 3 || isempty(pk_prop), pk_prop = cfg.pk.prop; end
    if nargin < 4 || isempty(pk_remi), pk_remi = cfg.pk.remi; end

    processor.effect_site_P = init_effect_site(pk_prop.ke0, cfg.bis_delay);
    processor.effect_site_R = init_effect_site(pk_remi.ke0, cfg.bis_delay);

    processor.pk_state_P = init_pk_state_online(pk_prop);
    processor.pk_state_R = init_pk_state_online(pk_remi);
    
    processor.ekf_van = init_ekf_4d(cfg, 'vanluchene');
    processor.ekf_gre = init_ekf_4d(cfg, 'greco');
    
    processor.online = struct('initialization_samples',10);
    
    processor.E0 = struct('x', cfg.E0_fixed);
    processor.BISmin = struct('x', cfg.BISmin_fixed);
    
    processor.induction_delay = struct('time',[],'bis',[],'Cp_P',[],'estimated',false,'ke0_from_delay',NaN);
end
