function state = init_ekf_loglin3d(cfg)
% Estimated in raw parameters, unlike the 1D, 2D and 4D. a0 is a location
% parameter in logit space and may be zero or negative, so a fractional step is
% undefined for it and this filter keeps an additive walk.
    state = struct();
    state.current_params = [0.5; 1.0; 0.2];  % [a0, aP, aR]
    state.initialized = true;
    state.sample_count = 0;
    
    state.P = diag([1.0, 0.5, 0.1].^2);
    state.Q = cfg.q * state.P;
    
    state.P_min = [0.01, 0.005, 0.001];
    state.P_max = [1.0, 0.5, 0.1];
    
    % FIM for identifiability
    state.FIM = zeros(3, 3);
    state.FIM_forgetting = cfg.fim_forgetting;
    state.ident_eigenvalue_ratio = cfg.ident_eigenvalue_ratio;
    state.FIM_condition = Inf;
    state.projection_count = 0;
    
    % Bounds and rate limits
    state.lb = [-5; 0.1; 0.01];
    state.ub = [5; 4; 1];
    state.param_rate_max = cfg.rate_cap * sqrt(diag(state.P));
    
    state.eps_P = 0.05;
    state.eps_R = 0.0001;
    
    state.R_base = cfg.R_base;

    state.param_hist = [];
    state.P_hist = [];
end
