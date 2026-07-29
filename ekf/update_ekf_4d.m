function [bis_pred, state] = update_ekf_4d(...
    state, bis_obs, CeP, CeR, E0, BISmin, model, cfg, learning_enabled, Rmult)

    state.sample_count = state.sample_count + 1;
    x = state.current_params(:);

    bis_pred = predict_bis_4d(x, CeP, CeR, E0, BISmin, model);

    if ~learning_enabled
        return;
    end

    % Chain rule to log-parameters: dBIS/dlog(theta) = dBIS/dtheta * theta.
    % The FIM accumulated from this is dimensionless.
    H = jacobian_4d(x, CeP, CeR, E0, BISmin, model) .* x(:)';

    R = Rmult * state.R_base;

    info = (H' * H) / R;
    state.FIM = state.FIM_forgetting * state.FIM + info;
    state.FIM_cum = state.FIM_cum + info;

    [V, D] = eig(state.FIM);
    eigenvalues = diag(D);
    [eigenvalues_sorted, sort_idx] = sort(real(eigenvalues), 'descend');
    V_sorted = V(:, sort_idx);

    state.FIM_eigenvalues = eigenvalues_sorted;
    state.FIM_condition = eigenvalues_sorted(1) / max(eigenvalues_sorted(end), 1e-12);

    lambda_max = eigenvalues_sorted(1);
    lambda_threshold = state.ident_eigenvalue_ratio * lambda_max;
    identifiable_mask = eigenvalues_sorted > lambda_threshold;
    n_identifiable = sum(identifiable_mask);

    if n_identifiable > 0 && n_identifiable < 4
        V_ident = V_sorted(:, identifiable_mask);
        P_proj = V_ident * V_ident';
        state.projection_count = state.projection_count + 1;
    elseif n_identifiable == 0
        return;
    else
        P_proj = eye(4);
    end

    P = state.P + state.Q;
    innovation = bis_obs - bis_pred;
    S = H * P * H' + R;
    if S < 1e-6, return; end
    K = P * H' / S;
    max_change = state.param_rate_max(:);
    dphi = max(-max_change, min(max_change, P_proj * (K * innovation)));
    x_new = x .* exp(dphi);
    x_new = max(cfg.lb(:), min(cfg.ub(:), x_new));

    I_KH = eye(4) - K * H;
    P_new = I_KH * P * I_KH' + K * R * K';
    P_new = (P_new + P_new') / 2;

    for ii = 1:4
        P_new(ii,ii) = max(state.P_min(ii), min(state.P_max(ii), P_new(ii,ii)));
    end

    state.current_params = x_new;
    state.P = P_new;
    state.update_count = state.update_count + 1;
    % P is on log-parameters, so its sd is already the relative standard error
    state.RSE_current = sqrt(diag(P_new)) * 100;
end
