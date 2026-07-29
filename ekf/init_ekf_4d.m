function state = init_ekf_4d(cfg, model)
% Estimates phi = log([C50P; C50R; gamma; beta]), so the FIM is dimensionless.
% Dividing the covariances and caps by p0 carries the raw tuning over unchanged.
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

    state.P     = diag(([1.0, 1.5, 0.4, 0.15]' ./ p0).^2);
    state.P_min = [0.01, 0.02, 0.005, 0.002]' ./ p0.^2;
    state.P_max = [2.0,  5.0,  0.8,   0.25 ]' ./ p0.^2;

    state.Q = cfg.q * state.P;

    state.FIM = zeros(4, 4);
    state.FIM_cum = zeros(4, 4);
    state.FIM_forgetting = cfg.fim_forgetting;
    state.FIM_condition = 1;
    state.FIM_eigenvalues = zeros(4, 1);

    state.ident_eigenvalue_ratio = cfg.ident_eigenvalue_ratio;

    state.param_rate_max = cfg.rate_cap * sqrt(diag(state.P));

    state.R_base = cfg.R_base;

    state.RSE_current = [];
    state.update_count = 0;
    state.projection_count = 0;
end
