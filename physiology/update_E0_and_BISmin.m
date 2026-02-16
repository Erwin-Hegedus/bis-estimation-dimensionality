function processor = update_E0_and_BISmin(processor, params, CeP, CeR, y_eff, data_quality, CpP, CpR)
    if data_quality > 2, return; end

    H = compute_hill_only(params, CeP, CeR, processor.ekf_van.model, []);
    H = max(0.01, min(0.99, H));

    diseq_P = CpP - CeP;
    diseq_R = CpR - CeR;
    dBIS_dCeP = compute_dBIS_dCe(params, CeP, CeR, processor.ekf_van.model, 'propofol');
    dBIS_dCeR = compute_dBIS_dCe(params, CeP, CeR, processor.ekf_van.model, 'remi');
    sigma2_model = (dBIS_dCeP * diseq_P)^2 + (dBIS_dCeR * diseq_R)^2;

    R = 50 + sigma2_model;

    if ~isfield(processor, 'endpoint_x')
        processor.endpoint_x = [93; 30];
        processor.endpoint_P = diag([16, 100]);
    end

    x = processor.endpoint_x;
    P = processor.endpoint_P;

    if H < 0.3
        C = [1-H, H];
        Q = diag([1e-4, 0]);
    elseif H > 0.7
        C = [1-H, H];
        Q = diag([0, 1e-4]);
    else
        processor.E0.x = x(1);
        processor.BISmin.x = x(2);
        return;
    end

    P_pred = P + Q;
    y_pred = C * x;
    innov = y_eff - y_pred;
    S = C * P_pred * C' + R;
    K = P_pred * C' / S;
    x_new = x + K * innov;
    P_new = (eye(2) - K * C) * P_pred;

    x_new(1) = max(85, min(100, x_new(1)));
    x_new(2) = max(10, min(35, x_new(2)));

    processor.endpoint_x = x_new;
    processor.endpoint_P = P_new;
    processor.E0.x = x_new(1);
    processor.E0.P = P_new(1,1);
    processor.BISmin.x = x_new(2);
    processor.BISmin.P = P_new(2,2);
end
