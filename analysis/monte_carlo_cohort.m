function coh = monte_carlo_cohort(results, cfg)
%MONTE_CARLO_COHORT  Cohort-wide equivalent-set fraction and L2 radius.
%   Runs the per-patient Monte Carlo (monte_carlo_equivalence_radius) over the
%   valid cohort (all seven models non-NaN, N=209) and aggregates the
%   equivalent-set fraction and normalized maximum L2 radius into
%   distributions (median / IQR / range). Patient 105 is retained in the
%   output for use as the illustrative example.

    rng(42);                                     % reproducible cohort run

    % Same valid cohort as the MAE tables: all seven models non-NaN.
    m = results.metrics;
    valid = ~isnan(m.pop.MAE) & ~isnan(m.kscale.MAE) & ~isnan(m.m2d.MAE) & ...
            ~isnan(m.m2d_fim.MAE) & ~isnan(m.loglin.MAE) & ~isnan(m.van.MAE) & ...
            ~isnan(m.gre.MAE);
    pids     = find(valid);
    box_diag = norm(cfg.ub(:) - cfg.lb(:));      % max possible L2 distance in the box

    np      = numel(pids);
    frac    = nan(np, 1);                         % fraction of the 500 sets that are equivalent
    rad_max = nan(np, 1);                         % maximum L2 radius, normalized by box_diag
    skipped = false(np, 1);                       % patients with too few valid samples

    for i = 1:np
        try
            e = monte_carlo_equivalence_radius(results, pids(i), cfg);
        catch
            skipped(i) = true;
            continue;
        end
        frac(i) = e.N_eq / e.N_total;             % 0 when no equivalent set was found
        if e.has_eq
            rad_max(i) = e.radius_max / box_diag;
        end
    end

    coh.pids        = pids;
    coh.box_diag    = box_diag;
    coh.frac        = frac;
    coh.radius_norm = rad_max;
    coh.skipped     = skipped;
    coh.frac_stats  = [median(frac,    'omitnan'), quantile(frac,    [.25 .75]), min(frac),                 max(frac)];
    coh.rad_stats   = [median(rad_max, 'omitnan'), quantile(rad_max, [.25 .75]), min(rad_max, [], 'omitnan'), max(rad_max, [], 'omitnan')];
end
