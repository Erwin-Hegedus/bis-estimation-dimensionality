function es = update_effect_site_with_delay_estimation(es, time_k, Cp_k, y_eff)
    if isnan(es.last_time)
        es.last_time = time_k;
        es.Ce_current = 0;
        es.Ce_delayed_output = 0;
        return;
    end
    
    dt = time_k - es.last_time;
    if dt <= 0, return; end
    
    alpha = exp(-es.ke0 * dt);
    es.Ce_current = alpha * es.Ce_current + (1 - alpha) * Cp_k;
    
    es.buffer_idx = mod(es.buffer_idx, es.buffer_size) + 1;
    es.Ce_buffer(es.buffer_idx) = es.Ce_current;
    es.BIS_buffer(es.buffer_idx) = y_eff;
    es.buffer_fill = min(es.buffer_fill + 1, es.buffer_size);
    
    es.Ce_history_idx = mod(es.Ce_history_idx, 60) + 1;
    es.Ce_history(es.Ce_history_idx) = es.Ce_current;
    
    if ~es.tau_d_locked && es.buffer_fill >= 120 && mod(es.buffer_fill, 60) == 0
        [tau_new, conf] = estimate_delay_xcorr(es);
        if ~isnan(tau_new)
            es.tau_d = 0.7 * es.tau_d + 0.3 * tau_new;
            es.tau_d = max(0, min(50, es.tau_d));
            es.tau_d_confidence = conf;
            if conf > 0.85
                es.tau_d_locked = true;
            end
        end
    end
    
    delay_samples = round(es.tau_d);
    delay_samples = max(0, min(59, delay_samples));
    lookup_idx = mod(es.Ce_history_idx - delay_samples - 1, 60) + 1;
    es.Ce_delayed_output = es.Ce_history(lookup_idx);
    es.last_time = time_k;
end
