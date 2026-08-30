function ekf = init_ekf_kgamma(cfg)
%INIT_EKF_KGAMMA  Shared potency plus Hill-slope EKF in log coordinates.
%   The prior spreads are inherited without retuning from the existing 1D
%   k estimator (sd 0.20) and the 4D gamma coordinate (sd 0.40 at gamma=2.8).

    p = cfg.population_params_van(:);

    ekf = struct();
    ekf.current_params = [1.0; p(3)];  % [k; gamma]

    sd_log = [0.20; 0.40 / p(3)];
    ekf.P = diag(sd_log.^2);
    ekf.P0 = ekf.P;
    ekf.Q = cfg.q * ekf.P0;

    ekf.P_min = [0.01^2; (0.005 / p(3))^2];
    ekf.P_max = [0.50^2; (0.80 / p(3))^2];

    ekf.lb = [0.3; cfg.lb(3)];
    ekf.ub = [3.0; cfg.ub(3)];
    ekf.rate_max = cfg.rate_cap * sqrt(diag(ekf.P0));

    ekf.FIM = zeros(2, 2);
    ekf.FIM_cum = zeros(2, 2);
    ekf.FIM_forgetting = cfg.fim_forgetting;
    ekf.ident_eigenvalue_ratio = cfg.ident_eigenvalue_ratio;
    ekf.FIM_condition = Inf;
    ekf.FIM_eigenvalues = zeros(2, 1);

    ekf.R_base = cfg.R_base;
    ekf.sample_count = 0;
    ekf.update_count = 0;
    ekf.projection_count = 0;
    ekf.rate_limit_count = 0;
end
