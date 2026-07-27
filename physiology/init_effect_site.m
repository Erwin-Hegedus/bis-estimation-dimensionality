function es = init_effect_site(ke0_per_min, delay_samples)
%INIT_EFFECT_SITE_RIGOROUS  First-order effect site with a fixed output delay.
%   ke0_per_min is the equilibration rate constant; delay_samples is the lag
%   applied to the effect-site output to account for BIS monitor smoothing.
    es = struct();
    es.ke0 = ke0_per_min / 60;
    es.Ce_current = 0;
    es.tau_d = delay_samples;
    es.Ce_history = zeros(es.tau_d + 1, 1);
    es.Ce_history_idx = 0;
    es.last_time = NaN;
    es.Ce_delayed_output = 0;
end
