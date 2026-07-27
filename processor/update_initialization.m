function processor = update_initialization(processor, bis_k, time_k)
    if processor.initialized, return; end
    
    processor.fast_init.bis_buffer(end+1) = bis_k;

    if length(processor.fast_init.bis_buffer) >= processor.fast_init.samples_needed
        processor.ekf_van = init_ekf_4d(processor.cfg, 'vanluchene');
        processor.ekf_gre = init_ekf_4d(processor.cfg, 'greco');
        
        processor.initialized = true;
    end
end
