function es = init_effect_site_rigorous(ke0_per_min)
    es = struct();
    es.ke0 = ke0_per_min / 60;
    es.Ce_current = 0;
    es.tau_d = 20;
    es.tau_d_prior = 20;
    es.tau_d_locked = false;
    es.tau_d_confidence = 0;
    es.buffer_size = 300;
    es.Ce_buffer = zeros(es.buffer_size, 1);
    es.BIS_buffer = zeros(es.buffer_size, 1);
    es.buffer_idx = 0;
    es.buffer_fill = 0;
    es.Ce_history = zeros(60, 1);
    es.Ce_history_idx = 0;
    es.last_time = NaN;
    es.Ce_delayed_output = 0;
end
