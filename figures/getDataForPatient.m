function [p_ce, r_ce, bis, time] = getDataForPatient(patientDataFinal, idx)
    p_ce = [];
    r_ce = [];
    bis = [];
    time = [];
    
    try
        if idx <= numel(patientDataFinal.Orchestra_PPF20_CE)
            p_struct = patientDataFinal.Orchestra_PPF20_CE(idx);
            if ~isempty(p_struct.Value)
                p_ce = p_struct.Value(:);
                time = p_struct.Time(:);
            end
        end
        
        if idx <= numel(patientDataFinal.Orchestra_RFTN20_CE)
            r_struct = patientDataFinal.Orchestra_RFTN20_CE(idx);
            if ~isempty(r_struct.Value)
                r_ce = r_struct.Value(:);  % FIXED - was r_ce.Value
            end
        end
        
        if idx <= numel(patientDataFinal.BIS_BIS)
            bis_struct = patientDataFinal.BIS_BIS(idx);
            if ~isempty(bis_struct.Value)
                bis = bis_struct.Value(:);
            end
        end
        
        if ~isempty(p_ce) && ~isempty(r_ce) && ~isempty(bis)
            n = min([numel(p_ce), numel(r_ce), numel(bis)]);
            p_ce = p_ce(1:n);
            r_ce = r_ce(1:n);
            bis = bis(1:n);
            if ~isempty(time)
                time = time(1:n);
            else
                time = (0:n-1)';
            end
        end
    catch
    end
end
