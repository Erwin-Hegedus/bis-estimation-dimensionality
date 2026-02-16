function [bis_pred, state] = update_level_rate_rls(...
    state, bis_obs, CeP, CeR, E0, BISmin, model, cfg, data_quality, learning_enabled, ...
    Cp_P, Cp_R, ke0P, ke0R)

    state.sample_count = state.sample_count + 1;
    x = state.current_params(:);
    C50P = x(1); C50R = x(2);
    
    bis_pred = predict_bis_with_endpoints(x, CeP, CeR, E0, BISmin, model);
    
    u_total = CeP / C50P + CeR * 1000 / C50R;
    if ~learning_enabled || data_quality > 2 || u_total < 0.15
        return;
    end
    
    H = compute_jacobian_with_endpoints(x, CeP, CeR, E0, BISmin, model);
    
    diseq_P = Cp_P - CeP;
    diseq_R = Cp_R - CeR;
    diseq_magnitude = abs(diseq_P)/max(C50P,0.1) + abs(diseq_R)*1000/max(C50R,0.1);
    R = state.R_base * (1 + state.R_disequilibrium_factor * diseq_magnitude^2);
    
    state.FIM = state.FIM_forgetting * state.FIM + (H' * H) / R;
    
    [V, D] = eig(state.FIM);
    eigenvalues = diag(D);
    [eigenvalues_sorted, sort_idx] = sort(real(eigenvalues), 'descend');
    V_sorted = V(:, sort_idx);
    
    state.FIM_eigenvalues = eigenvalues_sorted;
    state.FIM_eigenvectors = V_sorted;
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
    
    P = state.P;
    innovation = bis_obs - bis_pred;
    S = H * P * H' + R;
    if S < 1e-6, return; end
    K = P * H' / S;
    dx_raw = K * innovation;
    dx_projected = P_proj * dx_raw;
    max_change = state.param_rate_max(:);
    dx = max(-max_change, min(max_change, dx_projected));
    x_new = x + dx;
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
    state.RSE_current = sqrt(diag(P_new)) ./ max(abs(x_new), 0.1) * 100;
end
