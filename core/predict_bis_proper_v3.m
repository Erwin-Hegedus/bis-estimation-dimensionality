function y = predict_bis_proper_v3(p, CeP, CeR, E0, model, BISmin_fixed)
    epsv = 1e-12;
    C50P = p(1); C50R = p(2); gamma = p(3); syn_param = p(4);
    CeP_mgL = max(CeP, 0);
    CeR_ngmL = 1000 * max(CeR, 0);
    U_P = CeP_mgL / max(C50P, epsv);
    U_R = CeR_ngmL / max(C50R, epsv);
    
    if strcmpi(model, 'vanluchene')
        denom = U_P + U_R + epsv;
        theta = U_R / denom;
        beta = max(0, syn_param);
        U50_theta = 1 - beta * theta * (1 - theta);
        U_mix = (U_P + U_R) / max(0.01, U50_theta);
    else
        alpha = syn_param;
        U_mix = U_P + U_R + alpha * (U_P * U_R);
    end
    
    Hill_Effect = (U_mix^gamma) / (1 + U_mix^gamma);
    y = BISmin_fixed + (E0 - BISmin_fixed) * (1 - Hill_Effect);
    y = max(0, min(100, y));
end
