function [bis_pred, state] = update_ekf_loglin3d(state, bis_obs, CeP, CeR, E0, BISmin, cfg, Cp_P, Cp_R)
    state.sample_count = state.sample_count + 1;
    
    a0 = state.current_params(1);
    aP = state.current_params(2);
    aR = state.current_params(3);
    
    % Prediction
    xP = log(CeP + state.eps_P);
    xR = log(CeR + state.eps_R);
    
    Z = a0 + aP * xP + aR * xR;
    Z = max(-8, min(8, Z));
    H_sigmoid = 1 / (1 + exp(-Z));
    
    bis_pred = BISmin + (E0 - BISmin) * (1 - H_sigmoid);
    bis_pred = max(0, min(100, bis_pred));
    
    % Store history
    state.param_hist(end+1, :) = state.current_params';
    state.P_hist(end+1, :) = diag(state.P)';
    
    % Skip if low concentrations
    if CeP < 0.5 && CeR < 0.001
        return;
    end
    
    % === JACOBIAN ===
    dH_dZ = H_sigmoid * (1 - H_sigmoid);
    dBIS_dH = -(E0 - BISmin);
    H_jac = dBIS_dH * dH_dZ * [1, xP, xR];
    
    % === ADAPTIVE R ===
    diseq_P = Cp_P - CeP;
    diseq_R = Cp_R - CeR;
    diseq_mag = abs(diseq_P)/max(CeP + 0.5, 1) + abs(diseq_R)*1000/max(CeR*1000 + 1, 2);
    R = state.R_base * (1 + state.R_diseq * diseq_mag^2);
    
    % === UPDATE FIM ===
    state.FIM = state.FIM_forgetting * state.FIM + (H_jac' * H_jac) / R;
    
    % === EIGENVALUE DECOMPOSITION ===
    [V, D] = eig(state.FIM);
    eigenvalues = diag(D);
    [eigenvalues_sorted, sort_idx] = sort(real(eigenvalues), 'descend');
    V_sorted = V(:, sort_idx);
    
    state.FIM_condition = eigenvalues_sorted(1) / max(eigenvalues_sorted(end), 1e-12);
    
    % Identifiable directions
    lambda_threshold = state.ident_eigenvalue_ratio * eigenvalues_sorted(1);
    identifiable_mask = eigenvalues_sorted > lambda_threshold;
    n_identifiable = sum(identifiable_mask);
    
    % Build projection matrix
    if n_identifiable > 0 && n_identifiable < 3
        V_ident = V_sorted(:, identifiable_mask);
        P_proj = V_ident * V_ident';
        state.projection_count = state.projection_count + 1;
    elseif n_identifiable == 0
        return;
    else
        P_proj = eye(3);
    end
    
    % === EKF UPDATE ===
    innovation = bis_obs - bis_pred;
    
    if abs(innovation) < 0.5
        return;
    end
    
    P_pred = state.P + diag(diag(state.Q));
    S = H_jac * P_pred * H_jac' + R;
    
    if S < 1e-6
        return;
    end
    
    K = P_pred * H_jac' / S;
    dx_raw = K * innovation;
    
    % === PROJECT TO IDENTIFIABLE SUBSPACE ===
    dx_projected = P_proj * dx_raw;
    
    % === RATE LIMITING ===
    dx = max(-state.param_rate_max, min(state.param_rate_max, dx_projected));
    
    % === UPDATE ===
    x_new = state.current_params + dx;
    x_new = max(state.lb, min(state.ub, x_new));
    
    % === COVARIANCE ===
    I_KH = eye(3) - K * H_jac;
    P_new = I_KH * P_pred * I_KH' + K * R * K';
    P_new = (P_new + P_new') / 2;
    
    for ii = 1:3
        P_new(ii,ii) = max(state.P_min(ii), min(state.P_max(ii), P_new(ii,ii)));
    end
    
    state.current_params = x_new;
    state.P = P_new;
end
