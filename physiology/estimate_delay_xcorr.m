function [tau_opt, confidence] = estimate_delay_xcorr(es)
    if es.buffer_fill < es.buffer_size
        idx = 1:es.buffer_fill;
    else
        idx = [es.buffer_idx+1:es.buffer_size, 1:es.buffer_idx];
    end
    
    Ce = es.Ce_buffer(idx);
    BIS = es.BIS_buffer(idx);
    
    dCe = diff(Ce);
    dBIS = diff(BIS);
    
    valid = isfinite(dCe) & isfinite(dBIS);
    dCe = dCe(valid);
    dBIS = dBIS(valid);
    
    if length(dCe) < 60
        tau_opt = NaN;
        confidence = 0;
        return;
    end
    
    dCe = (dCe - mean(dCe)) / (std(dCe) + 1e-9);
    dBIS = (dBIS - mean(dBIS)) / (std(dBIS) + 1e-9);
    
    tau_grid = 0:1:50;
    n_tau = length(tau_grid);
    corr_vals = zeros(n_tau, 1);
    
    for ii = 1:n_tau
        lag = tau_grid(ii);
        n_overlap = length(dCe) - lag;
        if n_overlap < 30, continue; end
        r = corrcoef(dCe(1:n_overlap), dBIS(lag+1:lag+n_overlap));
        corr_vals(ii) = -r(1,2);
    end
    
    [max_corr, ~] = max(corr_vals);
    
    if max_corr < 0.2
        tau_opt = NaN;
        confidence = 0;
        return;
    end
    
    tau_opt = tau_grid(corr_vals == max_corr);
    tau_opt = tau_opt(1);
    
    sorted_corr = sort(corr_vals, 'descend');
    if length(sorted_corr) > 5
        noise_floor = mean(sorted_corr(5:end));
    else
        noise_floor = 0;
    end
    
    prominence = max_corr - noise_floor;
    confidence = min(1, max(0, prominence / 0.3));
    
    if confidence < 0.4
        tau_opt = NaN;
    end
end
