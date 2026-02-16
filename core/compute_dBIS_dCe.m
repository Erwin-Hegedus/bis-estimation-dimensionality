function dBIS_dCe = compute_dBIS_dCe(params, CeP, CeR, model_type, drug)
    C50P = params(1); C50R = params(2);
    gamma = params(3); syn = params(4);
    E0 = 93; BISmin = 32;
    eps = 1e-12;
    
    uP = CeP / max(C50P, eps);
    uR = 1000 * CeR / max(C50R, eps);
    
    if strcmpi(model_type, 'vanluchene')
        denom = uP + uR + eps;
        theta = uR / denom;
        U50 = 1 - syn * theta * (1 - theta);
        U = (uP + uR) / max(U50, 0.01);
        
        if strcmpi(drug, 'propofol')
            dTheta_du = -uR / (denom^2);
            dU50_dTheta = -syn * (1 - 2*theta);
            dU50_du = dU50_dTheta * dTheta_du;
            dU_du = (1/U50) - (uP + uR) * dU50_du / (U50^2);
            du_dCe = 1 / max(C50P, eps);
        else
            dTheta_du = uP / (denom^2);
            dU50_dTheta = -syn * (1 - 2*theta);
            dU50_du = dU50_dTheta * dTheta_du;
            dU_du = (1/U50) - (uP + uR) * dU50_du / (U50^2);
            du_dCe = 1000 / max(C50R, eps);
        end
    else
        U = uP + uR + syn * uP * uR;
        if strcmpi(drug, 'propofol')
            dU_du = 1 + syn * uR;
            du_dCe = 1 / max(C50P, eps);
        else
            dU_du = 1 + syn * uP;
            du_dCe = 1000 / max(C50R, eps);
        end
    end
    
    if U > eps
        Ug = U^gamma;
        dH_dU = gamma * (U^(gamma-1)) / ((1 + Ug)^2);
    else
        dH_dU = 0;
    end
    
    dBIS_dH = -(E0 - BISmin);
    dBIS_dCe = dBIS_dH * dH_dU * dU_du * du_dCe;
end
