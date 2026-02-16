function E0 = get_initial_E0_robust(bis_samples)
    valid = bis_samples(bis_samples > 0 & bis_samples <= 100);
    if length(valid) >= 5
        E0 = mean(valid(1:min(10, end)));
    else
        E0 = 85;
    end
    E0 = max(80, min(100, E0));
end
