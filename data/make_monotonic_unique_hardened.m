function [t2, v2] = make_monotonic_unique_hardened(t2, v2)
    t2 = double(t2(:));
    v2 = double(v2(:));
    g = isfinite(t2);
    t2 = t2(g);
    v2 = v2(g);
    
    if isempty(t2)
        t2 = [0; 1e6];
        v2 = [0; 0];
        return;
    end
    
    [t2, ord] = sort(t2);
    v2 = v2(ord);
    [t2, iu] = unique(t2, 'stable');
    v2 = v2(iu);
    
    if numel(t2) < 2
        t2 = [t2; t2(end) + 1e6];
        v2 = [v2; v2(end)];
    end
end
