function y = predict_bis_with_endpoints(p, CeP, CeR, E0, BISmin, model)
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
    
    H = (U_mix^gamma) / (1 + U_mix^gamma);
    y = BISmin + (E0 - BISmin) * (1 - H);
    y = max(0, min(100, y));
end
