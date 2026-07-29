function state = init_ekf_kscale(cfg)
% Estimates phi = log(k); P, Q and rate_max are fractional.
    state = struct();
    state.initialized   = true;
    state.sample_count  = 0;

    state.C50P_pop = cfg.population_params_van(1);
    state.C50R_pop = cfg.population_params_van(2);
    state.gamma    = cfg.population_params_van(3);
    state.beta     = cfg.population_params_van(4);

    state.k        = 1.0;
    state.P        = 0.20^2;
    state.Q        = cfg.q * state.P;
    state.rate_max = cfg.rate_cap * sqrt(state.P);
    state.R_base   = cfg.R_base;

    state.lb_k     = 0.3;
    state.ub_k     = 3.0;

    state.k_hist   = [];
    state.P_hist   = [];
end
