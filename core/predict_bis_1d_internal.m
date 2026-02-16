function bis = predict_bis_1d_internal(k, CeP, CeR, E0, BISmin, state)

    k_safe = max(k, 0.01);
    C50P_eff = state.C50P_pop / k_safe;
    C50R_eff = state.C50R_pop / k_safe;
    
    gamma = state.gamma;
    beta  = state.beta;
    
    
    uP = CeP / C50P_eff;
    uR = (CeR * 1000) / C50R_eff;
    
    sum_u = uP + uR;
    
    if sum_u < 1e-9
        bis = E0;
        return;
    end
    
    theta = uR / sum_u;
    U50_theta = 1 - beta * theta * (1 - theta);
    
    U = sum_u / U50_theta;
    
    
    if U <= 0
        E = 0;
    else
        U_gamma = U^gamma;
        E = U_gamma / (1 + U_gamma);
    end
    
    
    bis = BISmin + (E0 - BISmin) * (1 - E);
    bis = max(0, min(100, bis));
end
