function [bis_pred, state] = update_ekf_kscale(state, bis_obs, CeP, CeR, E0, BISmin_fixed, cfg, Cp_P, Cp_R)

    state.sample_count = state.sample_count + 1;
    
    k  = state.k;
    Pk = state.P; 

    bis_pred = predict_bis_1d_internal(k, CeP, CeR, E0, BISmin_fixed, state);
    
    if CeP < 0.5 && CeR < 0.1
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
    
    y_plus = predict_bis_1d_internal(k_test, CeP, CeR, E0, BISmin_fixed, state);
    
    Hk = (y_plus - bis_pred) / dk_step;
    
    diseq_P = Cp_P - CeP;
    diseq_R = Cp_R - CeR;
  
    C50P_curr = state.C50P_pop / k;
    C50R_curr = state.C50R_pop / k;
    diseq_mag = abs(diseq_P)/max(C50P_curr, 0.1) + abs(diseq_R)*1000/max(C50R_curr, 0.1);
    
    R_curr = state.R_base * (1 + state.R_diseq * diseq_mag^2);
    
    if abs(Hk) < 1e-4
        state.P = min(Pk + state.Q, 0.5^2);
        
        state.k_hist(end+1,1) = k;
        state.P_hist(end+1,1) = state.P;
        return; 
    end

    innovation = bis_obs - bis_pred;
    
    if abs(innovation) < 0.5
        state.k_hist(end+1,1) = k;
        state.P_hist(end+1,1) = Pk;
        return;
    end
    
    % Időbeli frissítés (Time Update)
    P_minus = Pk + state.Q; 
    
    % Mérési frissítés (Measurement Update)
    S  = Hk * P_minus * Hk' + R_curr;
    Kk = P_minus * Hk' / max(S, 1e-6); % Osztásvédelem
    
    % Nyers paraméter frissítés
    k_raw_change = Kk * innovation;
    
    % --- 8. Rate Limiting és Korlátok (A Szigorú Rész) ---
    % Nem engedünk azonnali ugrást, csak fizikai sebességgel
    max_step_per_sample = 0.05; % Pl. max 0.05 változás mintánként
    
    dk_clamped = max(-max_step_per_sample, min(max_step_per_sample, k_raw_change));
    k_new = k + dk_clamped;
    
    % Abszolút korlátok
    k_new = max(state.lb_k, min(state.ub_k, k_new));
    
    % Kovariancia frissítése (Joseph-form a stabilitásért)
    % P = (I - KH)P(I - KH)' + KRK' 
    % Skalár esetben egyszerűsítve: P = (1 - KH) * P_minus
    P_new = (1 - Kk * Hk) * P_minus;
    
    % Kovariancia "Safe Guard" (Ne legyen túl kicsi, se túl nagy)
    % Túl kicsi P -> A szűrő "elalszik" (nem tanul)
    % Túl nagy P -> A szűrő instabil
    P_new = max(1e-6, min(0.5^2, P_new));
    
    % --- 9. Állapot mentése ---
    state.k = k_new;
    state.P = P_new;
    state.k_hist(end+1,1) = k_new;
    state.P_hist(end+1,1) = P_new;
end
