function processor = init_online_processor_rigorous(cfg, start_time, pk_prop, pk_remi, demographics)
    processor = struct();
    processor.cfg = cfg;
    processor.sample_count = 0;
    processor.start_time = start_time;
    processor.current_time = start_time;
    processor.initialized = false;
    processor.art = init_artifact_gate_simple(cfg);
    processor.fast_init = struct('bis_buffer',[],'samples_needed',10,'ready',false);
    processor.ekf_k = init_ekf_kscale(cfg);
    processor.ekf_loglin = init_ekf_loglin3d(cfg);
    processor.ekf_2d = init_ekf_2d(cfg); 
    processor.ekf_2d_fim = init_ekf_2d_fim(cfg);
    
    processor.personalization.van = struct('detected',false(4,1), 'time',nan(4,1));
    processor.personalization.gre = struct('detected',false(4,1), 'time',nan(4,1));
    processor.personalization.baseline_van = [];
    processor.personalization.baseline_gre = [];
    
    [kP_demo, kR_demo] = compute_demographic_ke0_literature(demographics, cfg);
    
    processor.effect_site_P = init_effect_site_rigorous(kP_demo);
    processor.effect_site_R = init_effect_site_rigorous(kR_demo);
    
    processor.induction_status = struct('phase','pre_induction', ...
        'induction_start_sample',NaN,'induction_end_sample',NaN,'bis_history',[]);
    
    if nargin < 3 || isempty(pk_prop), pk_prop = cfg.pk.prop; end
    if nargin < 4 || isempty(pk_remi), pk_remi = cfg.pk.remi; end
    processor.pk_state_P = init_pk_state_online(pk_prop);
    processor.pk_state_R = init_pk_state_online(pk_remi);
    
    processor.ekf_van = init_rls_state(cfg, 'vanluchene');
    processor.ekf_gre = init_rls_state(cfg, 'greco');
    
    processor.online = struct('bis_buffer',[],'time_buffer',[], ...
        'initialization_samples',10,'E0_estimate',NaN,'R0_estimate',NaN);
    
    processor.E0 = struct('x', 93, 'buffer', [], 'locked', false);
    processor.BISmin = struct('x', cfg.bismin.init, 'P', cfg.bismin.P0);
    
    processor.induction_delay = struct('time',[],'bis',[],'Cp_P',[],'estimated',false,'ke0_from_delay',NaN);
end
