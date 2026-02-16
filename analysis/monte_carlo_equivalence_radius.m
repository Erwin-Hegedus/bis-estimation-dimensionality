function eq_info = monte_carlo_equivalence_radius(results, pid, cfg)
    r = results.raw(pid);
    bis_obs = r.bis(:);
    CeP     = r.CeP_trajectory(:);
    CeR     = r.CeR_trajectory(:);

    E0_idx     = find(~isnan(r.E0_trajectory),1,'last');
    BISmin_idx = find(~isnan(r.BISmin_trajectory),1,'last');

    if isempty(E0_idx),     E0 = 93; else,     E0 = r.E0_trajectory(E0_idx); end
    if isempty(BISmin_idx), BISmin = cfg.bismin.init; else, BISmin = r.BISmin_trajectory(BISmin_idx); end

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

    theta_good = [];
    radius_good = [];

    for ii = 1:Nsamples
        theta_rand = lb + (ub - lb) .* rand(4,1);
        [mae_i, ~] = mae_van4d(theta_rand, CeP, CeR, E0, BISmin, bis_obs);

        if mae_i <= mae_ref * tol_factor
            theta_good(:, end+1) = theta_rand;
            radius_good(end+1)   = norm(theta_rand - theta_ref, 2);
        end
    end

    if isempty(theta_good)
        eq_info.has_eq      = false;
        eq_info.theta_ref   = theta_ref(:).';
        eq_info.mae_ref     = mae_ref;
        eq_info.radius_max  = 0;
        eq_info.radius_mean = 0;
        eq_info.N_eq        = 0;
        return;
    end

    eq_info.has_eq      = true;
    eq_info.theta_ref   = theta_ref(:).';
    eq_info.mae_ref     = mae_ref;
    eq_info.theta_good  = theta_good.';
    eq_info.radius_all  = radius_good(:);
    eq_info.radius_max  = max(radius_good);
    eq_info.radius_mean = mean(radius_good);
    eq_info.N_eq        = numel(radius_good);
end
