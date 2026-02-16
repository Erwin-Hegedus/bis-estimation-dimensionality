function processor = update_initialization_clean(processor, bis_k, time_k)
    if processor.initialized, return; end
    
    processor.fast_init.bis_buffer(end+1) = bis_k;
    processor.online.bis_buffer(end+1) = bis_k;
    
    if length(processor.fast_init.bis_buffer) >= processor.fast_init.samples_needed
        E0 = get_initial_E0_robust(processor.fast_init.bis_buffer);
        processor.online.E0_estimate = E0;
        processor.E0.x = E0;
        processor.E0.P = processor.cfg.e0.P0;
        
        processor.ekf_van = init_rls_state(processor.cfg, 'vanluchene');
        processor.ekf_gre = init_rls_state(processor.cfg, 'greco');
        
        processor.initialized = true;
    end
end
