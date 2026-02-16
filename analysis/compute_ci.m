function ci = compute_ci(data, confidence)
% Compute confidence interval using t-distribution
    data = data(~isnan(data));
    n = length(data);
    m = mean(data);
    s = std(data);
    alpha = 1 - confidence;
    t_crit = tinv(1 - alpha/2, n - 1);
    margin = t_crit * s / sqrt(n);
    ci = [m - margin, m + margin];
end
