function eq_info = monte_carlo_equivalence_radius(results, pid, cfg)
    r = results.raw(pid);
    bis_obs = r.bis(:);
    CeP     = r.CeP_trajectory(:);
    CeR     = r.CeR_trajectory(:);

    E0_idx     = find(~isnan(r.E0_trajectory),1,'last');
    BISmin_idx = find(~isnan(r.BISmin_trajectory),1,'last');

    if isempty(E0_idx),     E0 = 93; else,     E0 = r.E0_trajectory(E0_idx); end
    if isempty(BISmin_idx), BISmin = cfg.BISmin_fixed; else, BISmin = r.BISmin_trajectory(BISmin_idx); end

    valid = ~isnan(bis_obs) & ~isnan(CeP) & ~isnan(CeR);
    bis_obs = bis_obs(valid);
    CeP     = CeP(valid);
    CeR     = CeR(valid);

    if numel(bis_obs) < 100
        error('Too few valid samples for patient %d', pid);
    end

    if isfield(r,'Xhist_van') && ~isempty(r.Xhist_van)
        last_idx = find(any(~isnan(r.Xhist_van),1),1,'last');
        theta_ref = r.Xhist_van(:, last_idx);
    else
        theta_ref = cfg.population_params_van(:);
    end

    [mae_ref, ~] = mae_van4d(theta_ref, CeP, CeR, E0, BISmin, bis_obs);

    Nsamples   = 500;
    tol_factor = 1.10;
    lb = cfg.lb(:);
    ub = cfg.ub(:);

    rng(42);
    theta_good  = nan(4, Nsamples);
    radius_good = nan(1, Nsamples);
    radius_norm = nan(1, Nsamples);   % component-normalised, unit-invariant
    d2_raw      = zeros(4,1);         % per-parameter share of raw squared distance
    d2_norm     = zeros(4,1);
    cnt = 0;
    cnt_fixed = 0;                    % accepted under a fixed +1 BIS tolerance

    for ii = 1:Nsamples
        theta_rand = lb + (ub - lb) .* rand(4,1);
        [mae_i, ~] = mae_van4d(theta_rand, CeP, CeR, E0, BISmin, bis_obs);

        if mae_i <= mae_ref + 1
            cnt_fixed = cnt_fixed + 1;
        end

        if mae_i <= mae_ref * tol_factor
            cnt = cnt + 1;
            d = theta_rand - theta_ref;
            z = d ./ (ub - lb);
            theta_good(:, cnt) = theta_rand;
            radius_good(cnt)   = norm(d, 2);
            radius_norm(cnt)   = norm(z, 2);
            d2_raw  = d2_raw  + d.^2;
            d2_norm = d2_norm + z.^2;
        end
    end
    theta_good  = theta_good(:, 1:cnt);
    radius_good = radius_good(1:cnt);
    radius_norm = radius_norm(1:cnt);

    if isempty(theta_good)
        eq_info.has_eq      = false;
        eq_info.theta_ref   = theta_ref(:).';
        eq_info.mae_ref     = mae_ref;
        eq_info.radius_max  = 0;
        eq_info.radius_mean = 0;
        eq_info.norm_max    = 0;
        eq_info.norm_p95    = 0;
        eq_info.norm_median = 0;
        eq_info.share_norm  = nan(1,4);
        eq_info.N_eq        = 0;
        eq_info.N_eq_fixed  = cnt_fixed;
        eq_info.N_total     = Nsamples;
        return;
    end

    eq_info.has_eq      = true;
    eq_info.theta_ref   = theta_ref(:).';
    eq_info.mae_ref     = mae_ref;
    eq_info.theta_good  = theta_good.';
    eq_info.radius_all  = radius_good(:);
    eq_info.radius_max  = max(radius_good);
    eq_info.radius_mean = mean(radius_good);

    eq_info.norm_all    = radius_norm(:);
    eq_info.norm_max    = max(radius_norm);
    eq_info.norm_p95    = quantile(radius_norm, 0.95);
    eq_info.norm_median = median(radius_norm);
    eq_info.norm_bound  = norm(ones(4,1));
    eq_info.share_raw   = (d2_raw  / sum(d2_raw)).';
    eq_info.share_norm  = (d2_norm / sum(d2_norm)).';

    eq_info.N_eq        = numel(radius_good);
    eq_info.N_eq_fixed  = cnt_fixed;
    eq_info.N_total     = Nsamples;
end
