function [y_eff, Rmult_ext, art, data_quality] = artifact_gate_simple(art, t, y, cfg)
    y_eff = y;
    Rmult_ext = 1;
    data_quality = 1;
    
    if isnan(art.last_t)
        art.last_t = t;
        art.last_y = y;
        art.last_valid_y = y;
        return;
    end
    
    dt = max(0.1, t - art.last_t);
    dr = (y - art.last_y) / dt;
    
    if y <= 0 || y > 100 || ~isfinite(y) || abs(dr) > 30
        data_quality = 3;
        y_eff = art.last_valid_y;
        Rmult_ext = 1000;
        art.bad_count = art.bad_count + 1;
    else
        art.bad_count = 0;
        if (y < 10 || y > 95) || (abs(dr) > 15)
            data_quality = 2;
            Rmult_ext = 5;
        end
        if abs(dr) > cfg.artifact.bis_rate_thresh
            y_eff = art.last_y + sign(dr) * cfg.artifact.bis_rate_thresh * dt;
            Rmult_ext = max(Rmult_ext, 3);
        end
        if data_quality <= 2 && y >= 10 && y <= 90
            art.last_valid_y = y_eff;
        end
    end
    
    art.last_t = t;
    art.last_y = y;
end
