function [p_fdr, is_sig] = fdr_correction(p_values, alpha)
% Benjamini-Hochberg FDR correction
    m = length(p_values);
    [p_sorted, sort_idx] = sort(p_values);
    
    % BH critical values
    bh_crit = (1:m)' / m * alpha;
    
    % Find largest k where p(k) <= k/m * alpha
    k_max = find(p_sorted <= bh_crit, 1, 'last');
    
    if isempty(k_max)
        is_sig = false(m, 1);
    else
        is_sig = false(m, 1);
        is_sig(sort_idx(1:k_max)) = true;
    end
    
    % Adjusted p-values
    p_fdr = nan(m, 1);
    p_fdr(sort_idx) = min(1, p_sorted .* m ./ (1:m)');
    
    % Enforce monotonicity
    for i = m-1:-1:1
        p_fdr(sort_idx(i)) = min(p_fdr(sort_idx(i)), p_fdr(sort_idx(i+1)));
    end
end
