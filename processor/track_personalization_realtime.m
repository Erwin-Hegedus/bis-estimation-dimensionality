function processor = track_personalization_realtime(processor, time_k, cfg)
    if isempty(processor.personalization.baseline_van)
        processor.personalization.baseline_van = cfg.population_params_van(:);
        processor.personalization.baseline_gre = cfg.population_params_gre(:);
    end
    if processor.sample_count < 180, return; end
    
    current_van = processor.ekf_van.current_params;
    for ii = 1:4
        if ~processor.personalization.van.detected(ii)
            delta = abs(current_van(ii) - processor.personalization.baseline_van(ii));
            pct = delta / abs(processor.personalization.baseline_van(ii));
            if pct > 0.05
                processor.personalization.van.detected(ii) = true;
                processor.personalization.van.time(ii) = (time_k - processor.start_time) / 60;
            end
        end
    end
    
    current_gre = processor.ekf_gre.current_params;
    for ii = 1:4
        if ~processor.personalization.gre.detected(ii)
            delta = abs(current_gre(ii) - processor.personalization.baseline_gre(ii));
            pct = delta / abs(processor.personalization.baseline_gre(ii));
            if pct > 0.05
                processor.personalization.gre.detected(ii) = true;
                processor.personalization.gre.time(ii) = (time_k - processor.start_time) / 60;
            end
        end
    end
end
