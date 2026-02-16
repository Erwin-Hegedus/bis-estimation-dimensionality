function bis = predict_bis_proof_internal(kP, kR, CeP, CeR, E0, cfg)
    % Simplified prediction for landscape generation
    C50P_pop = cfg.population_params_van(1);
    C50R_pop = cfg.population_params_van(2);
    gamma    = cfg.population_params_van(3);
    beta     = cfg.population_params_van(4);
    BISmin   = 30; 
    
    uP = (CeP * kP) / max(C50P_pop, 1e-6);
    uR = (CeR * 1000 * kR) / max(C50R_pop, 1e-6);
    
    U_sum = uP + uR;
    if U_sum < 1e-9, bis = E0; return; end
    
    theta = uR / U_sum;
    U50 = 1 - beta * theta * (1 - theta);
    U = U_sum / U50;
    
    E = (U^gamma) / (1 + U^gamma);
    bis = BISmin + (E0 - BISmin) * (1 - E);
    bis = max(0, min(100, bis));
end
