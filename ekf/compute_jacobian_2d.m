function H = compute_jacobian_2d(kP, kR, CeP, CeR, E0, BISmin, cfg)
% Compute Jacobian [∂BIS/∂kP, ∂BIS/∂kR] numerically
    
    delta = 1e-4;
    
    % Baseline prediction
    y0 = predict_bis_2d_model(kP, kR, CeP, CeR, E0, BISmin, cfg);
    
    % Perturb kP
    y_kP = predict_bis_2d_model(kP + delta, kR, CeP, CeR, E0, BISmin, cfg);
    dBIS_dkP = (y_kP - y0) / delta;
    
    % Perturb kR
    y_kR = predict_bis_2d_model(kP, kR + delta, CeP, CeR, E0, BISmin, cfg);
    dBIS_dkR = (y_kR - y0) / delta;
    
    H = [dBIS_dkP, dBIS_dkR];
end
