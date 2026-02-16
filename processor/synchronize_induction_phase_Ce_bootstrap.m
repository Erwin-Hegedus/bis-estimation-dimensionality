function [prop_sync, remi_sync, bis_sync, time_sync, shift_sec, n_eff, status] = ...
    synchronize_induction_phase_Ce_bootstrap(bis, prop, remi, time, ke0_pop)

    bis  = bis(:);
    prop = prop(:);
    remi = remi(:);
    time = time(:);

    N = min([numel(bis), numel(prop), numel(remi), numel(time)]);
    bis  = bis(1:N);
    prop = prop(1:N);
    remi = remi(1:N);
    time = time(1:N);

    shift_sec = 0;
    status = 'Aligned';

  
    dt = median(diff(time));
    if ~isfinite(dt) || dt <= 0 || dt > 5
        error('Invalid time base (must be seconds)');
    end

    n_win = min(N, ceil(15*60/dt));
    if n_win < 120
        prop_sync = prop;
        remi_sync = remi;
        bis_sync  = bis;
        time_sync = time;
        n_eff = N;
        status = 'Window too short';
        return;
    end

    if ~isscalar(ke0_pop) || ke0_pop <= 0 || ke0_pop > 1
        error('ke0_pop must be scalar in [1/s]');
    end

    CeP = zeros(n_win,1);
    alpha = exp(-ke0_pop * dt);

    for k = 2:n_win
        CeP(k) = alpha * CeP(k-1) + (1 - alpha) * prop(k-1);
    end

    if std(CeP) < 1e-6
        prop_sync = prop;
        remi_sync = remi;
        bis_sync  = bis;
        time_sync = time;
        n_eff = N;
        status = 'Flat CeP proxy';
        return;
    end

    
    bis_target = 100 - bis(1:n_win);

    max_lag = ceil(3*60/dt); 
    [c,lags] = xcorr( ...
        bis_target - mean(bis_target), ...
        CeP - mean(CeP), ...
        max_lag);

    [~,idx] = max(c);
    lag = lags(idx);
    shift_sec = lag * dt;

    if lag > 0
        % Inputs too early → delay inputs
        cut = lag;
        prop_sync = prop(1+cut:end);
        remi_sync = remi(1+cut:end);
        bis_sync  = bis(1:end-cut);
        time_sync = time(1:end-cut);
        status = 'Inputs delayed';

    elseif lag < 0
        % Inputs too late → advance inputs
        cut = abs(lag);
        prop_sync = prop(1:end-cut);
        remi_sync = remi(1:end-cut);
        bis_sync  = bis(1+cut:end);
        time_sync = time(1:end-cut);
        status = 'Inputs advanced';

    else
        prop_sync = prop;
        remi_sync = remi;
        bis_sync  = bis;
        time_sync = time;
    end

    
    n_eff = min([numel(prop_sync), numel(remi_sync), ...
                 numel(bis_sync),  numel(time_sync)]);

    prop_sync = prop_sync(1:n_eff);
    remi_sync = remi_sync(1:n_eff);
    bis_sync  = bis_sync(1:n_eff);
    time_sync = time_sync(1:n_eff);
end
