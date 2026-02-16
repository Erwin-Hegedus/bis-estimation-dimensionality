function [mae, sse] = mae_van4d(theta, CeP, CeR, E0, BISmin, bis_obs)
    if theta(1) < 1.0 || theta(1) > 12.0 || ...
       theta(2) < 2.0 || theta(2) > 40.0 || ...
       theta(3) < 1.0 || theta(3) > 6.0  || ...
       theta(4) < 0.3 || theta(4) > 1.8
        mae = Inf;
        sse = Inf;
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
    mae = mean(abs(diff_val), 'omitnan');
    sse = sum(diff_val.^2, 'omitnan');
end
