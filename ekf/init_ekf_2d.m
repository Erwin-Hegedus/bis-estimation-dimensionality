function ekf = init_ekf_2d(cfg)
    ekf = struct();
    ekf.current_params = [1.0; 1.0];  % [kP, kR]
    
    % Covariance in NATURAL coordinates (FIM-whitened)
    ekf.P = diag([0.3^2, 0.3^2]);
    ekf.Q = cfg.q * ekf.P;
    
    % FIM accumulation
    ekf.FIM = zeros(2, 2);
    ekf.FIM_forgetting = cfg.fim_forgetting;
    ekf.ident_eigenvalue_ratio = cfg.ident_eigenvalue_ratio;
    
    % Natural gradient parameters
    ekf.lambda_ridge = 0.01;  % Tikhonov regularization
    ekf.use_natural_gradient = true;
    
    % Bounds and rate limits
    ekf.lb = [0.4; 0.4];
    ekf.ub = [2.5; 2.5];
    ekf.rate_max = [0.02; 0.02];
    
    % Diagnostics
    ekf.n_updates = 0;
    ekf.FIM_condition = Inf;
end
