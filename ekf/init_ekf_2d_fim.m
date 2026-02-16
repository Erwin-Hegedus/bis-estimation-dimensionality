function ekf = init_ekf_2d_fim(cfg)
    ekf = struct();
    ekf.current_params = [1.0; 1.0];  % [kP, kR]
    
    % Covariance in NATURAL coordinates (FIM-whitened)
    ekf.P = diag([0.3^2, 0.3^2]);
    ekf.Q = diag([0.001^2, 0.001^2]);  % Small process noise
    
    % FIM accumulation
    ekf.FIM = eye(2) * 0.01;  % Small initialization
    ekf.FIM_forgetting = 0.99;  % Faster forgetting for 2D
    
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
