function H_jac = compute_jacobian_with_endpoints(x, CeP, CeR, E0, BISmin, model)
    C50P = x(1);
    C50R = x(2);
    gamma = x(3);
    beta = x(4);
    
    eps = 1e-8;
    
    uP = CeP / max(C50P, eps);
    uR = CeR * 1000 / max(C50R, eps);
    
    if strcmpi(model, 'vanluchene')
        U_sum = uP + uR + eps;
        theta = uR / U_sum;
        U50 = 1 - beta * theta * (1 - theta);
        U = U_sum / max(U50, 0.01);
        
        duP_dC50P = -CeP / C50P^2;
        duR_dC50R = -CeR * 1000 / C50R^2;
        
        dtheta_duP = -uR / U_sum^2;
        dtheta_duR = uP / U_sum^2;
        dU50_dtheta = -beta * (1 - 2*theta);
        
        dU_duP = (1/U50) - U_sum * dU50_dtheta * dtheta_duP / U50^2;
        dU_duR = (1/U50) - U_sum * dU50_dtheta * dtheta_duR / U50^2;
        
        dU_dC50P = dU_duP * duP_dC50P;
        dU_dC50R = dU_duR * duR_dC50R;
        dU_dbeta = -U_sum * (-theta * (1-theta)) / U50^2;
    else
        U = uP + uR + beta * uP * uR;
        
        duP_dC50P = -CeP / C50P^2;
        duR_dC50R = -CeR * 1000 / C50R^2;
        
        dU_dC50P = (1 + beta * uR) * duP_dC50P;
        dU_dC50R = (1 + beta * uP) * duR_dC50R;
        dU_dbeta = uP * uR;
    end
    
    if U < eps
        H_jac = [0, 0, 0, 0];
        return;
    end
    
    Ug = U^gamma;
    dH_dU = gamma * U^(gamma-1) / (1 + Ug)^2;
    dH_dgamma = Ug * log(max(U, eps)) / (1 + Ug)^2;
    
    dBIS_dH = -(E0 - BISmin);
    
    H_jac = zeros(1, 4);
    H_jac(1) = dBIS_dH * dH_dU * dU_dC50P;
    H_jac(2) = dBIS_dH * dH_dU * dU_dC50R;
    H_jac(3) = dBIS_dH * dH_dgamma;
    H_jac(4) = dBIS_dH * dH_dU * dU_dbeta;
end
