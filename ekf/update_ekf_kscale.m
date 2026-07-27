function [bis_pred, state] = update_ekf_kscale(state, bis_obs, CeP, CeR, E0, BISmin, cfg, ...
                                               learning_enabled, Rmult, Cp_P, Cp_R)

    state.sample_count = state.sample_count + 1;

    k  = state.k;
    Pk = state.P;

    bis_pred = predict_bis_1d_internal(k, CeP, CeR, E0, BISmin, state);

    if ~learning_enabled
        state.k_hist(end+1,1) = k;
        state.P_hist(end+1,1) = Pk;
        return;
    end

    dk_step = 0.01;
    k_test = k + dk_step;
    

    if k_test > state.ub_k
        k_test = k - dk_step;
        dk_step = -dk_step;
    end
    
    y_plus = predict_bis_1d_internal(k_test, CeP, CeR, E0, BISmin, state);
    
    Hk = (y_plus - bis_pred) / dk_step;
    
    diseq_P = Cp_P - CeP;
    diseq_R = Cp_R - CeR;
  
    C50P_curr = state.C50P_pop / k;
    C50R_curr = state.C50R_pop / k;
    diseq_mag = abs(diseq_P)/max(C50P_curr, 0.1) + abs(diseq_R)*1000/max(C50R_curr, 0.1);
    
    R_curr = Rmult * state.R_base * (1 + state.R_diseq * diseq_mag^2);
    
    if abs(Hk) < 1e-4
        state.P = min(Pk + state.Q, 0.5^2);
        
        state.k_hist(end+1,1) = k;
        state.P_hist(end+1,1) = state.P;
        return; 
    end

    innovation = bis_obs - bis_pred;

    P_minus = Pk + state.Q;

    S  = Hk * P_minus * Hk' + R_curr;
    Kk = P_minus * Hk' / max(S, 1e-6);

    max_step_per_sample = 0.05;
    dk_clamped = max(-max_step_per_sample, min(max_step_per_sample, Kk * innovation));
    k_new = k + dk_clamped;
    k_new = max(state.lb_k, min(state.ub_k, k_new));

    % Joseph form, scalar case: P = (1 - KH)^2 P + K^2 R
    P_new = (1 - Kk * Hk)^2 * P_minus + Kk^2 * R_curr;
    P_new = max(1e-6, min(0.5^2, P_new));

    state.k = k_new;
    state.P = P_new;
    state.k_hist(end+1,1) = k_new;
    state.P_hist(end+1,1) = P_new;
end
