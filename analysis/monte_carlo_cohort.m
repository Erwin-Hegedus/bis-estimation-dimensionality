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

    np        = numel(pids);
    frac      = nan(np, 1);                       % fraction of the 500 sets that are equivalent
    frac_fix  = nan(np, 1);                       % same, under a fixed +1 BIS tolerance
    rad_max   = nan(np, 1);                       % maximum L2 radius, normalized by box_diag
    norm_max  = nan(np, 1);                       % component-normalised, unit-invariant
    norm_p95  = nan(np, 1);                       % 95th pct: max over 500 draws is unstable
    share     = nan(np, 4);                       % per-parameter share of normalised distance
    skipped   = false(np, 1);                     % patients with too few valid samples

    for i = 1:np
        try
            e = monte_carlo_equivalence_radius(results, pids(i), cfg);
        catch
            skipped(i) = true;
            continue;
        end
        frac(i)     = e.N_eq / e.N_total;         % 0 when no equivalent set was found
        frac_fix(i) = e.N_eq_fixed / e.N_total;
        if e.has_eq
            rad_max(i)  = e.radius_max / box_diag;
            norm_max(i) = e.norm_max / e.norm_bound;
            norm_p95(i) = e.norm_p95 / e.norm_bound;
            share(i,:)  = e.share_norm;
        end
    end

    coh.pids        = pids;
    coh.box_diag    = box_diag;
    coh.frac        = frac;
    coh.frac_fixed  = frac_fix;
    coh.radius_norm = rad_max;
    coh.norm_max    = norm_max;
    coh.norm_p95    = norm_p95;
    coh.share_norm  = share;
    coh.skipped     = skipped;

    st = @(x) [median(x,'omitnan'), quantile(x,[.25 .75]), min(x,[],'omitnan'), max(x,[],'omitnan')];
    coh.frac_stats       = st(frac);
    coh.frac_fixed_stats = st(frac_fix);
    coh.rad_stats        = st(rad_max);     % raw, unit-dependent: kept for continuity
    coh.norm_max_stats   = st(norm_max);    % unit-invariant
    coh.norm_p95_stats   = st(norm_p95);
    coh.share_median     = median(share, 1, 'omitnan');
end
