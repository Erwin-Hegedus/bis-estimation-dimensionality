function bis = compute_bis_2d(theta2d, CeP, CeR, E0, BISmin, cfg)
% COMPUTE_BIS_2D - Vanluchene model with separate kP, kR scaling
%
% Effective C50s:
%   C50P_eff = C50P_pop / kP
%   C50R_eff = C50R_pop / kR
%
% Vanluchene stimulus:
%   uP = CeP * kP / C50P_pop
%   uR = CeR * 1000 * kR / C50R_pop
%   U_sum = uP + uR
%   theta_u = uR / U_sum
%   U50 = 1 - beta * theta_u * (1 - theta_u)
%   U = U_sum / U50
%
% Hill: H = U^gamma / (1 + U^gamma)
% BIS = BISmin + (E0 - BISmin) * (1 - H)

    kP = theta2d(1);
    kR = theta2d(2);
    
    C50P_pop = cfg.population_params_van(1);
    C50R_pop = cfg.population_params_van(2);
    gamma    = cfg.population_params_van(3);
    beta     = cfg.population_params_van(4);
    
    eps = 1e-12;
    
    % Vanluchene stimulus with separate scaling
    uP = CeP * kP / max(C50P_pop, eps);
    uR = CeR * 1000 * kR / max(C50R_pop, eps);
    
    U_sum = uP + uR + eps;
    theta_u = uR / U_sum;
    U50 = 1 - beta * theta_u * (1 - theta_u);
    U = U_sum / max(U50, 0.01);
    
    % Hill
    Ug = U^gamma;
    H = Ug / (1 + Ug);
    
    % BIS
    bis = BISmin + (E0 - BISmin) * (1 - H);
    bis = max(0, min(100, bis));
end
