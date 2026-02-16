function pred = predict_bis_2d_model(kP, kR, CeP, CeR, E0, BISmin, cfg)
% Predict BIS using 2D model (separate kP, kR scaling)
    
    % Population parameters
    C50P_pop = cfg.population_params_van(1);
    C50R_pop = cfg.population_params_van(2);
    gamma = cfg.population_params_van(3);
    beta = cfg.population_params_van(4);
    
    % Scaled loads
    uP = (CeP * kP) / C50P_pop;
    uR = (CeR * 1000 * kR) / C50R_pop;
    
    % Bouillon interaction
    epsilon = 1e-9;
    theta_u = uR / (uP + uR + epsilon);
    U50 = 1 - beta * theta_u * (1 - theta_u);
    U = (uP + uR) / U50;
    
    % Hill equation
    H = U^gamma / (1 + U^gamma);
    
    % BIS
    pred = BISmin + (E0 - BISmin) * (1 - H);
    pred = max(0, min(100, pred));
end
