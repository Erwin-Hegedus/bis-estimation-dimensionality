function [pred, ekf] = update_ekf_2d_fim(ekf, y, CeP, CeR, E0, BISmin, cfg, ...
                                          data_quality, learning_enabled, CpP, CpR, ke0P, ke0R)
    
    kP = ekf.current_params(1);
    kR = ekf.current_params(2);
    
    % Predict
    pred = predict_bis_2d_model(kP, kR, CeP, CeR, E0, BISmin, cfg);
    
    if ~learning_enabled || data_quality > 2
        return;
    end
    
    % === 1. COMPUTE SENSITIVITY (JACOBIAN) ===
    H = compute_jacobian_2d(kP, kR, CeP, CeR, E0, BISmin, cfg);  % 1x2
    
    % Adaptive measurement noise
    D_P = abs(CpP - CeP) / max(CeP, 0.5);
    D_R = abs(CpR - CeR*1000) / max(CeR*1000, 1);
    R = cfg.R_base * (1 + cfg.R_disequilibrium_factor * (D_P + D_R)^2);
    
    % === 2. UPDATE FIM (Exponentially Weighted) ===
    outer_product = (H' * H) / R;
    ekf.FIM = ekf.FIM_forgetting * ekf.FIM + (1 - ekf.FIM_forgetting) * outer_product;
    ekf.n_updates = ekf.n_updates + 1;
    
    % === 3. EIGENVALUE-BASED PROJECTION (like 4D) ===
    [V, D] = eig(ekf.FIM);
    eigenvalues = diag(D);
    [eigenvalues_sorted, sort_idx] = sort(real(eigenvalues), 'descend');
    V_sorted = V(:, sort_idx);
    
    ekf.FIM_condition = eigenvalues_sorted(1) / max(eigenvalues_sorted(2), 1e-12);
    
    % Identifiability threshold (like 4D)
    lambda_max = eigenvalues_sorted(1);
    lambda_threshold = 0.01 * lambda_max;  % Same ratio as 4D
    identifiable_mask = eigenvalues_sorted > lambda_threshold;
    n_identifiable = sum(identifiable_mask);
    
    % Build projection matrix
    if n_identifiable > 0 && n_identifiable < 2
        V_ident = V_sorted(:, identifiable_mask);
        P_proj = V_ident * V_ident';
    elseif n_identifiable == 0
        return;  % No update if nothing identifiable
    else
        P_proj = eye(2);  % Full update if both identifiable
    end
    
    % === 4. STANDARD EKF UPDATE ===
    innovation = y - pred;
    P_pred = ekf.P + ekf.Q;
    S = H * P_pred * H' + R;
    if S < 1e-6, return; end
    K = P_pred * H' / S;
    
    % Raw update
    delta_theta_raw = K * innovation;
    
    % Project to identifiable subspace (like 4D)
    delta_theta = P_proj * delta_theta_raw;
    
    % === 5. RATE LIMITING ===
    delta_theta = max(-ekf.rate_max, min(ekf.rate_max, delta_theta));
    
    % === 6. PARAMETER UPDATE ===
    ekf.current_params = ekf.current_params + delta_theta;
    ekf.current_params = max(ekf.lb, min(ekf.ub, ekf.current_params));
    
    % === 7. COVARIANCE UPDATE ===
    I_KH = eye(2) - K * H;
    P_new = I_KH * P_pred * I_KH' + K * R * K';
    P_new = (P_new + P_new') / 2;
    P_min = 0.01^2;
    P_max = 0.40^2;
    for i = 1:2
        P_new(i,i) = max(P_min, min(P_max, P_new(i,i)));
    end
    ekf.P = P_new;
end
