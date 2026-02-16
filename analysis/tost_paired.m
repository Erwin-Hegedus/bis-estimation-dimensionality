function [p_tost, is_equivalent] = tost_paired(x, y, delta)
% Two One-Sided Tests (TOST) for equivalence
    diff_data = x - y;
    n = length(diff_data);
    m = mean(diff_data);
    se = std(diff_data) / sqrt(n);
    
    % Test H0: diff >= delta (upper bound)
    t_upper = (m - delta) / se;
    p_upper = tcdf(t_upper, n - 1);
    
    % Test H0: diff <= -delta (lower bound)
    t_lower = (m + delta) / se;
    p_lower = 1 - tcdf(t_lower, n - 1);
    
    % TOST p-value is the maximum of the two
    p_tost = max(p_upper, p_lower);
    
    % Equivalent if both null hypotheses rejected at alpha = 0.05
    is_equivalent = p_tost < 0.05;
end
