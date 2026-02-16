function [kP, kR] = compute_demographic_ke0_literature(demographics, cfg)
    kP = 0.35;
    kR = 0.50;
    
    if ~isempty(demographics) && isfield(demographics, 'Age')
        age = demographics.Age;
        if isfinite(age) && age > 0
            age_factor = 1.0 - 0.006 * max(0, age - 40);
            age_factor = max(0.65, min(1.0, age_factor));
            kP = kP * age_factor;
            kR = kR * age_factor;
        end
    end
    
    kP = max(0.20, min(0.50, kP));
    kR = max(0.35, min(0.65, kR));
end
