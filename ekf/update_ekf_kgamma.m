function [pred, ekf] = update_ekf_kgamma(ekf, y, CeP, CeR, E0, BISmin, cfg, ...
                                               learning_enabled, Rmult)
%UPDATE_EKF_KGAMMA  Update shared potency k and Hill slope gamma.
%   Both states are estimated in log coordinates. The Jacobian is an exact
%   reparameterization of JACOBIAN_4D, keeping this reduced model on the same
%   numerical path as the other Hill estimators.

    ekf.sample_count = ekf.sample_count + 1;
    k = ekf.current_params(1);
    gamma = ekf.current_params(2);

    pred = predict_bis_kgamma_model(k, gamma, CeP, CeR, E0, BISmin, cfg);
    if ~learning_enabled
        return;
    end

    p = cfg.population_params_van(:);
    theta_eff = [p(1) / k; p(2) / k; gamma; p(4)];
    H4_log = jacobian_4d(theta_eff, CeP, CeR, E0, BISmin, ...
        'vanluchene') .* theta_eff';
    H = [-(H4_log(1) + H4_log(2)), H4_log(3)];

    R = Rmult * ekf.R_base;
    info = (H' * H) / R;
    ekf.FIM = ekf.FIM_forgetting * ekf.FIM + info;
    ekf.FIM_cum = ekf.FIM_cum + info;

    [V, D] = eig(ekf.FIM);
    [eigenvalues, order] = sort(real(diag(D)), 'descend');
    V = V(:, order);
    ekf.FIM_eigenvalues = eigenvalues;
    ekf.FIM_condition = eigenvalues(1) / max(eigenvalues(end), 1e-12);

    threshold = ekf.ident_eigenvalue_ratio * eigenvalues(1);
    identifiable = eigenvalues > threshold;
    n_identifiable = sum(identifiable);
    if n_identifiable == 0
        return;
    elseif n_identifiable < 2
        V_ident = V(:, identifiable);
        P_proj = V_ident * V_ident';
        ekf.projection_count = ekf.projection_count + 1;
    else
        P_proj = eye(2);
    end

    P_pred = ekf.P + ekf.Q;
    innovation = y - pred;
    S = H * P_pred * H' + R;
    if S < 1e-6
        return;
    end
    K = P_pred * H' / S;

    dphi_unclamped = P_proj * (K * innovation);
    if any(abs(dphi_unclamped) > ekf.rate_max)
        ekf.rate_limit_count = ekf.rate_limit_count + 1;
    end
    dphi = max(-ekf.rate_max, min(ekf.rate_max, dphi_unclamped));

    x_new = ekf.current_params .* exp(dphi);
    x_new = max(ekf.lb, min(ekf.ub, x_new));

    I_KH = eye(2) - K * H;
    P_new = I_KH * P_pred * I_KH' + K * R * K';
    P_new = (P_new + P_new') / 2;
    for ii = 1:2
        P_new(ii, ii) = max(ekf.P_min(ii), min(ekf.P_max(ii), P_new(ii, ii)));
    end

    ekf.current_params = x_new;
    ekf.P = P_new;
    ekf.update_count = ekf.update_count + 1;
end
