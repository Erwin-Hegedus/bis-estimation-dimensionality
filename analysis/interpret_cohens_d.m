function interp = interpret_cohens_d(d)
% Interpret Cohen's d effect size
    d = abs(d);
    if d < 0.2
        interp = 'negligible';
    elseif d < 0.5
        interp = 'small';
    elseif d < 0.8
        interp = 'medium';
    else
        interp = 'large';
    end
end
