function cfg = default_cfg()
% Configuration of the published run. Single source of truth: run_analysis.m and
% the reproduction drivers both take cfg from here.

    cfg = struct();

    cfg.ke0P = 0.456;  % Schnider population (DO NOT inflate)
    cfg.ke0R = 0.595;  % Minto, at the reference age of 40
    cfg.q = 3e-5;      % process noise scale; every estimator uses Q = cfg.q * P0
    cfg.BISmin_fixed = 30;
    cfg.E0_fixed = 93;
    cfg.bis_delay = 25;  % samples at 1 Hz, BIS monitor smoothing lag

    cfg.population_params_van = [3.5, 5.0, 2.8, 0.8];
    cfg.population_params_gre = [3.2, 4.5, 2.8, 0.9];

    % Physiological Bounds
    cfg.lb = [1.0,  2.0,  1.0,  0.3];
    cfg.ub = [12.0, 40.0, 6.0, 1.8];

    cfg.rate_cap = 1/15;  % per-sample step limit; cap = cfg.rate_cap * sqrt(diag(P0))

    cfg.fim_forgetting = 0.995;
    cfg.ident_eigenvalue_ratio = 0.01;  % FIM directions below this fraction of lambda_1 are projected out

    cfg.R_base = 400;  % adapted only by the artifact gate: R = Rmult * R_base

    cfg.pk.prop = struct('V1',4.27,'V2',18.9,'V3',238,'CL',1.89,'Q2',1.29,'Q3',0.836,'input_conc',20,'ke0',cfg.ke0P);
    cfg.pk.remi = struct('V1',5.1,'V2',9.82,'V3',5.42,'CL',2.6,'Q2',2.05,'Q3',0.076,'input_conc',0.02,'ke0',cfg.ke0R);

    cfg.target_band = [40,60];
    cfg.min_drug_effect = 0.15;  % minimum CeP/C50P + 1000*CeR/C50R for a sample to inform
    cfg.online = struct('initialization_samples',10);
    cfg.artifact = struct('bis_min',0,'bis_max',100,'bis_jump_thresh',25,'bis_rate_thresh',15);
end
