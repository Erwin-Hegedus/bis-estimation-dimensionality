function state = init_ekf_kscale(cfg)
    state = struct();
    state.initialized   = true;
    state.sample_count  = 0;

    state.C50P_pop = cfg.population_params_van(1);
    state.C50R_pop = cfg.population_params_van(2);
    state.gamma    = cfg.population_params_van(3);
    state.beta     = cfg.population_params_van(4);

    state.k        = 1.0;
    state.P        = 0.20^2;
    state.Q        = (0.002)^2;
    state.R_base   = cfg.R_base;
    state.R_diseq  = cfg.R_disequilibrium_factor;

    state.lb_k     = 0.3;
    state.ub_k     = 3.0;

    state.k_hist   = [];
    state.P_hist   = [];
end
