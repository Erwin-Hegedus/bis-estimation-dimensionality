function es = update_effect_site(es, time_k, Cp_k)
%UPDATE_EFFECT_SITE  Advance the effect site and return its delayed output.
%   Ce follows dCe/dt = ke0*(Cp - Ce), discretised exactly over dt. The output
%   is Ce delayed by es.tau_d samples.
    if isnan(es.last_time)
        es.last_time = time_k;
        return;
    end

    dt = time_k - es.last_time;
    if dt <= 0, return; end

    alpha = exp(-es.ke0 * dt);
    es.Ce_current = alpha * es.Ce_current + (1 - alpha) * Cp_k;

    n = numel(es.Ce_history);
    es.Ce_history_idx = mod(es.Ce_history_idx, n) + 1;
    es.Ce_history(es.Ce_history_idx) = es.Ce_current;

    es.Ce_delayed_output = es.Ce_history(mod(es.Ce_history_idx - es.tau_d - 1, n) + 1);
    es.last_time = time_k;
end
