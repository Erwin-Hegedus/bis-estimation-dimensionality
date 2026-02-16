function state = init_ekf_2d(cfg)

    state = struct();
    state.initialized = true;
    state.sample_count = 0;
    
    % Population parameters
    state.C50P_pop = cfg.population_params_van(1);
    state.C50R_pop = cfg.population_params_van(2);
    state.gamma    = cfg.population_params_van(3);
    state.beta     = cfg.population_params_van(4);
    
    % Initial parameters
    state.kP = 1.0;
    state.kR = 1.0;
    state.current_params = [state.kP; state.kR];
    state.theta_prior = state.current_params;
    
    % Covariance - same structure as LogLin
    state.P = diag([0.3, 0.3].^2);
    state.P_min = [0.02, 0.02];
    state.P_max = [0.5, 0.5];
    
    % Process noise - same as LogLin
    state.Q = diag([0.002, 0.002].^2);
    
    % Measurement noise
    state.R_base = cfg.R_base;
    state.R_diseq = cfg.R_disequilibrium_factor;
    
    % Rate limits - same as LogLin
    state.param_rate_max = [0.01, 0.01];
    
    % Diagnostics
    state.param_hist = [];
    state.P_hist = [];
end
