function state = init_ekf_loglin3d(cfg)
    state = struct();
    state.current_params = [0.5; 1.0; 0.2];  % [a0, aP, aR]
    state.initialized = true;
    state.sample_count = 0;
    
    state.P = diag([1.0, 0.5, 0.1].^2);
    state.Q = diag([0.005, 0.002, 0.0005].^2);
    
    state.P_min = [0.01, 0.005, 0.001];
    state.P_max = [1.0, 0.5, 0.1];
    
    % FIM for identifiability
    state.FIM = eye(3) * 0.01;
    state.FIM_forgetting = 0.995;
    state.ident_eigenvalue_ratio = 0.01;
    state.FIM_condition = Inf;
    state.projection_count = 0;
    
    % Bounds and rate limits
    state.lb = [-5; 0.1; 0.01];
    state.ub = [5; 4; 1];
    state.param_rate_max = [0.02; 0.01; 0.005];
    
    state.eps_P = 0.05;
    state.eps_R = 0.0001;
    
    state.R_base = 400;
    state.R_diseq = 15;
    
    state.param_hist = [];
    state.P_hist = [];
end
