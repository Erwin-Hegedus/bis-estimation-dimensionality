function pred = predict_bis_kgamma_model(k, gamma, CeP, CeR, E0, BISmin, cfg)
%PREDICT_BIS_KGAMMA_MODEL  Reduced Hill model on the two FIM-supported axes.
%   k is the existing shared potency factor: C50P = C50P_pop/k and
%   C50R = C50R_pop/k. gamma is estimated, while the population potency
%   ratio and Bouillon interaction parameter beta remain fixed.

    p = cfg.population_params_van(:);
    k_safe = max(k, 0.01);
    theta = [p(1) / k_safe; p(2) / k_safe; gamma; p(4)];
    pred = predict_bis_4d(theta, CeP, CeR, E0, BISmin, 'vanluchene');
end
