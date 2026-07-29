function ekf = init_ekf_2d(cfg)
% Estimates phi = log([kP; kR]); P, Q and rate_max are fractional.
    ekf = struct();
    ekf.current_params = [1.0; 1.0];  % [kP, kR]

    ekf.P = diag([0.3^2, 0.3^2]);
    ekf.Q = cfg.q * ekf.P;

    ekf.FIM = zeros(2, 2);
    ekf.FIM_forgetting = cfg.fim_forgetting;
    ekf.ident_eigenvalue_ratio = cfg.ident_eigenvalue_ratio;

    % Bounds and rate limits
    ekf.lb = [0.4; 0.4];
    ekf.ub = [2.5; 2.5];
    ekf.rate_max = cfg.rate_cap * sqrt(diag(ekf.P));

    % Diagnostics
    ekf.n_updates = 0;
    ekf.FIM_condition = Inf;
end
