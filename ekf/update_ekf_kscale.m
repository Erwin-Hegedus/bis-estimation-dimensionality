function [bis_pred, state] = update_ekf_kscale(state, bis_obs, CeP, CeR, E0, BISmin, cfg, ...
                                               learning_enabled, Rmult)

    state.sample_count = state.sample_count + 1;

    k  = state.k;
    Pk = state.P;

    bis_pred = predict_bis_1d_internal(k, CeP, CeR, E0, BISmin, state);

    if ~learning_enabled
        state.k_hist(end+1,1) = k;
        state.P_hist(end+1,1) = Pk;
        return;
    end

    % k scales both population C50s, so d/dlog(k) = -sum of the 4D's first two.
    theta_eff = [state.C50P_pop / k; state.C50R_pop / k; state.gamma; state.beta];
    H4 = jacobian_4d(theta_eff, CeP, CeR, E0, BISmin, 'vanluchene') .* theta_eff';
    Hphi = -(H4(1) + H4(2));

    R_curr = Rmult * state.R_base;

    if abs(Hphi) < 1e-4
        state.P = min(Pk + state.Q, 0.5^2);

        state.k_hist(end+1,1) = k;
        state.P_hist(end+1,1) = state.P;
        return;
    end

    innovation = bis_obs - bis_pred;

    P_minus = Pk + state.Q;

    S  = Hphi * P_minus * Hphi' + R_curr;
    Kk = P_minus * Hphi' / max(S, 1e-6);

    dphi_clamped = max(-state.rate_max, min(state.rate_max, Kk * innovation));
    k_new = k * exp(dphi_clamped);
    k_new = max(state.lb_k, min(state.ub_k, k_new));

    % Joseph form, scalar case: P = (1 - KH)^2 P + K^2 R
    P_new = (1 - Kk * Hphi)^2 * P_minus + Kk^2 * R_curr;
    P_new = max(1e-6, min(0.5^2, P_new));

    state.k = k_new;
    state.P = P_new;
    state.k_hist(end+1,1) = k_new;
    state.P_hist(end+1,1) = P_new;
end
