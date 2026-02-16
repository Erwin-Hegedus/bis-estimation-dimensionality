function pk = update_pk_state_online(pk, time_k, rate_k)
    if isnan(pk.last_time)
        pk.last_time = time_k;
        return;
    end
    dt = time_k - pk.last_time;
    if dt <= 0, return; end
    
    input_rate = rate_k * pk.params.input_conc / 3600;
    A = [-(pk.k10 + pk.k12 + pk.k13), pk.k21, pk.k31;
         pk.k12, -pk.k21, 0;
         pk.k13, 0, -pk.k31];
    B = [1; 0; 0];
    
    n_sub = max(1, ceil(dt));
    h = dt / n_sub;
    for s = 1:n_sub
        pk.x = pk.x + h * (A * pk.x + B * input_rate);
        pk.x = max(pk.x, 0);
    end
    pk.Cp = pk.x(1) / pk.params.V1;
    pk.last_time = time_k;
end
