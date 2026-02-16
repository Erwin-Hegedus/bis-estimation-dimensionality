function state = init_rls_state(cfg, model)
    state = struct();
    state.model = model;
    state.initialized = true;
    state.sample_count = 0;
    
    if strcmpi(model, 'vanluchene')
        p0 = cfg.population_params_van(:);
    else
        p0 = cfg.population_params_gre(:);
    end
    
    state.current_params = p0;
    state.theta_prior = p0;
    
    state.P = diag([1.0, 1.5, 0.4, 0.15].^2);
    state.P_min = [0.01, 0.02, 0.005, 0.002];
    state.P_max = [2.0, 5.0, 0.8, 0.25];
    
    state.Q = zeros(4, 4);
    
    state.FIM = zeros(4, 4);
    state.FIM_forgetting = 0.995;
    state.FIM_condition = 1;
    state.FIM_eigenvalues = zeros(4, 1);
    state.FIM_eigenvectors = eye(4);
    
    state.ident_eigenvalue_ratio = 0.01;
    state.ident_condition_threshold = 100;
    
    state.param_rate_max = [0.003, 0.008, 0.002, 0.001];
    
    state.R_base = 400;
    state.R_disequilibrium_factor = 15;
    
    state.RSE_current = [];
    state.update_count = 0;
    state.projection_count = 0;
end
