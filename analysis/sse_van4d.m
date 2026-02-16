function sse = sse_van4d(theta, CeP, CeR, E0, BISmin, bis_obs, theta_ref, lambda_reg)
    C50P = theta(1);
    C50R = theta(2);
    gamma = theta(3);
    beta  = theta(4);

    if C50P < 1.0 || C50P > 12.0 || ...
       C50R < 2.0 || C50R > 40.0 || ...
       gamma < 1.0 || gamma > 6.0 || ...
       beta  < 0.3 || beta  > 1.8
        sse = 1e9;
        return;
    end

    CeP = CeP(:);
    CeR = CeR(:);
    bis_obs = bis_obs(:);
    N = numel(bis_obs);

    bis_pred = nan(N,1);
    for k = 1:N
        if isnan(CeP(k)) || isnan(CeR(k)) || isnan(bis_obs(k))
            continue;
        end
        bis_pred(k) = predict_bis_with_endpoints(theta(:), CeP(k), CeR(k), E0, BISmin, 'vanluchene');
    end

    diff_val = bis_obs - bis_pred;
    sse_data = sum(diff_val.^2, 'omitnan');
    dtheta = theta(:) - theta_ref(:);
    sse_reg = lambda_reg * sum(dtheta.^2);
    sse = sse_data + sse_reg;
end
