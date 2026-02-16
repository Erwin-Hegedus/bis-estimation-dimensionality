function [bis_pred, state] = update_ekf_2d(state, y_eff, CeP, CeR, E0, BISmin, cfg, Cp_P, Cp_R)

    state.sample_count = state.sample_count + 1;
    
    % Current parameters
    kP = state.current_params(1);
    kR = state.current_params(2);
    
    % === 1. PREDICT BIS ===
    bis_pred = compute_bis_2d(state.current_params, CeP, CeR, E0, BISmin, cfg);
    
    % === 2. LEARNING GATE ===
    if CeP < 0.5 && CeR < 0.001
        state.param_hist(end+1, :) = state.current_params';
        state.P_hist(end+1, :) = diag(state.P)';
        return;
    end
    
    % === 3. ADAPTIVE MEASUREMENT NOISE  ===
    diseq_P = Cp_P - CeP;
    diseq_R = Cp_R - CeR;
    diseq_mag = abs(diseq_P)/max(CeP + 0.5, 1) + abs(diseq_R)*1000/max(CeR*1000 + 1, 2);
    R = state.R_base * (1 + state.R_diseq * diseq_mag^2);
    
    % === 4. COMPUTE JACOBIAN (numerical) ===
    delta = 1e-4;
    
    params_pP = state.current_params; params_pP(1) = kP + delta;
    params_pR = state.current_params; params_pR(2) = kR + delta;
    
    bis_pP = compute_bis_2d(params_pP, CeP, CeR, E0, BISmin, cfg);
    bis_pR = compute_bis_2d(params_pR, CeP, CeR, E0, BISmin, cfg);
    
    H_jac = [(bis_pP - bis_pred) / delta, (bis_pR - bis_pred) / delta];
    
    % === 5. INNOVATION ===
    innovation = y_eff - bis_pred;
    
    % Skip tiny innovations 
    if abs(innovation) < 0.5
        state.param_hist(end+1, :) = state.current_params';
        state.P_hist(end+1, :) = diag(state.P)';
        return;
    end
    
    % === 6. EKF UPDATE  ===
    P = state.P + diag(diag(state.Q));
    
    S = H_jac * P * H_jac' + R;
    if S < 1e-6
        state.param_hist(end+1, :) = state.current_params';
        state.P_hist(end+1, :) = diag(state.P)';
        return;
    end
    
    K = P * H_jac' / S;  % 2x1
    
    % Parameter update
    dx = K * innovation;
    
    % === 7. RATE LIMITING ===
    dx = max(-state.param_rate_max(:), min(state.param_rate_max(:), dx));
    
    x_new = state.current_params + dx;
    
    % === 8. PARAMETER BOUNDS ===
    x_new(1) = max(0.4, min(2.5, x_new(1)));  % kP
    x_new(2) = max(0.4, min(2.5, x_new(2)));  % kR
    
    % === 9. COVARIANCE UPDATE (Joseph form) ===
    I_KH = eye(2) - K * H_jac;
    P_new = I_KH * P * I_KH' + K * R * K';
    P_new = (P_new + P_new') / 2;
    
    % Apply bounds
    for i = 1:2
        P_new(i,i) = max(state.P_min(i), min(state.P_max(i), P_new(i,i)));
    end
    
    % === 10. STORE ===
    state.current_params = x_new;
    state.P = P_new;
    
    state.param_hist(end+1, :) = x_new';
    state.P_hist(end+1, :) = diag(P_new)';
end
