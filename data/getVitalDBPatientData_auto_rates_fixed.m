function [bis, remi_rate, prop_rate, time] = getVitalDBPatientData_auto_rates_fixed(index, matfile)
    persistent T;
    if isempty(T)
        S = load(matfile, 'patientDataFinal');
        T = S.patientDataFinal;
    end
    
    if istable(T)
        nrows = height(T);
    else
        nrows = numel(T);
    end
    
    if index < 1 || index > nrows
        error('Index %d out of range', index);
    end
    
    bis_series = series_of('BIS_BIS', index);
    prop_series = [];
    for name = {'Orchestra_PPF20_RATE', 'PPF20_RATE', 'Propofol_RATE'}
        s = series_of(name{1}, index);
        if ~isempty(s)
            prop_series = s;
            break;
        end
    end
    remi_series = [];
    for name = {'Orchestra_RFTN20_RATE', 'RFTN20_RATE', 'Remifentanil_RATE'}
        s = series_of(name{1}, index);
        if ~isempty(s)
            remi_series = s;
            break;
        end
    end
    
    [time, bis] = tv_of(bis_series, 'BIS');
    [pt, pv] = tv_of(prop_series, 'PROP');
    [rt, rv] = tv_of(remi_series, 'REMI');
    
    [pt, pv] = make_monotonic_unique_hardened(pt, pv);
    [rt, rv] = make_monotonic_unique_hardened(rt, rv);
    
    prop_rate = interp1(pt, pv, time, 'previous', 0);
    remi_rate = interp1(rt, rv, time, 'previous', 0);
    
    v = isfinite(time) & isfinite(bis);
    time = time(v);
    bis = bis(v);
    prop_rate = prop_rate(v);
    remi_rate = remi_rate(v);
    
    function s = series_of(varname, idx)
        s = [];
        if istable(T) && ismember(varname, T.Properties.VariableNames)
            col = T.(varname);
            if iscell(col)
                s = col{idx};
            else
                s = col(idx);
            end
        end
    end
    
    function [t, v2] = tv_of(s, ~)
        if isstruct(s)
            t = double(s.Time);
            v2 = double(s.Value);
        elseif istable(s)
            t = double(s.Time);
            v2 = double(s.Value);
        elseif istimetable(s)
            t = seconds(s.Time - s.Time(1));
            v2 = double(s{:,1});
        else
            t = [];
            v2 = [];
        end
    end
end
