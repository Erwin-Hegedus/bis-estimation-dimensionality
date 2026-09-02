function state = init_ekf_4d(cfg, model)
% Test build. cfg.c50_variant selects the arm; absent or 'shipped' reproduces
% the canonical estimator exactly.
%   'shipped'      projection on the forgetting FIM, full 4-parameter state
%   'fim_cum'      same, but the projection uses the cumulative FIM
%   'q_norank'     shipped, with Q carrying no component along log(C50R/C50P)
%   'freeze_ratio' 3-parameter state psi = log[kappa, gamma, beta],
%                  C50P = kappa*p0(1), C50R = kappa*p0(2)
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
    state.p0 = p0;

    state.variant = 'shipped';
    if isfield(cfg, 'c50_variant'), state.variant = cfg.c50_variant; end

    P_full     = diag(([1.0, 1.5, 0.4, 0.15]' ./ p0).^2);
    Pmin_full  = [0.01, 0.02, 0.005, 0.002]' ./ p0.^2;
    Pmax_full  = [2.0,  5.0,  0.8,   0.25 ]' ./ p0.^2;

    e_ratio  = [1; -1; 0; 0] / sqrt(2);
    state.Pi = eye(4) - e_ratio * e_ratio';

    if strcmp(state.variant, 'freeze_ratio')
        % Reduced basis. The common-scale variance is the mean of the two C50
        % log-variances it replaces, so the marginal prior on each C50 is
        % unchanged to two decimals.
        state.T = [1 0 0; 1 0 0; 0 1 0; 0 0 1];
        d = diag(P_full);
        state.P     = diag([mean(d(1:2)); d(3); d(4)]);
        state.P_min = [mean(Pmin_full(1:2)); Pmin_full(3); Pmin_full(4)];
        state.P_max = [mean(Pmax_full(1:2)); Pmax_full(3); Pmax_full(4)];
        state.kappa_lb = max(cfg.lb(1)/p0(1), cfg.lb(2)/p0(2));
        state.kappa_ub = min(cfg.ub(1)/p0(1), cfg.ub(2)/p0(2));
        state.psi = zeros(3,1);
        nd = 3;
    else
        state.P     = P_full;
        state.P_min = Pmin_full;
        state.P_max = Pmax_full;
        nd = 4;
    end

    state.Q = cfg.q * state.P;
    if strcmp(state.variant, 'q_norank')
        state.Q = state.Pi * state.Q * state.Pi;
    end

    state.FIM = zeros(nd, nd);
    state.FIM_cum = zeros(nd, nd);
    state.FIM_forgetting = cfg.fim_forgetting;
    state.FIM_condition = 1;
    state.FIM_eigenvalues = zeros(nd, 1);

    state.ident_eigenvalue_ratio = cfg.ident_eigenvalue_ratio;

    state.param_rate_max = cfg.rate_cap * sqrt(diag(state.P));

    state.R_base = cfg.R_base;

    state.RSE_current = [];
    state.update_count = 0;
    state.projection_count = 0;
end
