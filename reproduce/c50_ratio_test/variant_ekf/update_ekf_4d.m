function [bis_pred, state] = update_ekf_4d(...
    state, bis_obs, CeP, CeR, E0, BISmin, model, cfg, learning_enabled, Rmult)
% Test build. See init_ekf_4d for the arms. With cfg.c50_variant absent this is
% line-for-line the shipped estimator.

    state.sample_count = state.sample_count + 1;
    x = state.current_params(:);

    bis_pred = predict_bis_4d(x, CeP, CeR, E0, BISmin, model);

    if ~learning_enabled
        return;
    end

    reduced = strcmp(state.variant, 'freeze_ratio');
    nd = 4; if reduced, nd = 3; end

    % Chain rule to log-parameters: dBIS/dlog(theta) = dBIS/dtheta * theta.
    H = jacobian_4d(x, CeP, CeR, E0, BISmin, model) .* x(:)';
    if reduced, H = H * state.T; end

    R = Rmult * state.R_base;

    info = (H' * H) / R;
    state.FIM = state.FIM_forgetting * state.FIM + info;
    state.FIM_cum = state.FIM_cum + info;

    if strcmp(state.variant, 'fim_cum')
        F = state.FIM_cum;
    else
        F = state.FIM;
    end

    [V, D] = eig(F);
    eigenvalues = diag(D);
    [eigenvalues_sorted, sort_idx] = sort(real(eigenvalues), 'descend');
    V_sorted = V(:, sort_idx);

    state.FIM_eigenvalues = eigenvalues_sorted;
    state.FIM_condition = eigenvalues_sorted(1) / max(eigenvalues_sorted(end), 1e-12);

    lambda_max = eigenvalues_sorted(1);
    lambda_threshold = state.ident_eigenvalue_ratio * lambda_max;
    identifiable_mask = eigenvalues_sorted > lambda_threshold;
    n_identifiable = sum(identifiable_mask);

    if n_identifiable > 0 && n_identifiable < nd
        V_ident = V_sorted(:, identifiable_mask);
        P_proj = V_ident * V_ident';
        state.projection_count = state.projection_count + 1;
    elseif n_identifiable == 0
        return;
    else
        P_proj = eye(nd);
    end

    P = state.P + state.Q;
    innovation = bis_obs - bis_pred;
    S = H * P * H' + R;
    if S < 1e-6, return; end
    K = P * H' / S;
    max_change = state.param_rate_max(:);
    dphi = max(-max_change, min(max_change, P_proj * (K * innovation)));

    if reduced
        psi_new = state.psi + dphi;
        kappa = min(state.kappa_ub, max(state.kappa_lb, exp(psi_new(1))));
        g     = min(cfg.ub(3), max(cfg.lb(3), state.p0(3) * exp(psi_new(2))));
        b     = min(cfg.ub(4), max(cfg.lb(4), state.p0(4) * exp(psi_new(3))));
        state.psi = [log(kappa); log(g/state.p0(3)); log(b/state.p0(4))];
        x_new = [state.p0(1)*kappa; state.p0(2)*kappa; g; b];
    else
        x_new = x .* exp(dphi);
        x_new = max(cfg.lb(:), min(cfg.ub(:), x_new));
    end

    I_KH = eye(nd) - K * H;
    P_new = I_KH * P * I_KH' + K * R * K';
    P_new = (P_new + P_new') / 2;

    for ii = 1:nd
        P_new(ii,ii) = max(state.P_min(ii), min(state.P_max(ii), P_new(ii,ii)));
    end

    state.current_params = x_new;
    state.P = P_new;
    state.update_count = state.update_count + 1;
    state.RSE_current = sqrt(diag(P_new)) * 100;
end
