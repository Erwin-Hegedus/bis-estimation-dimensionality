function cases = pick_case_studies_v3(vec, results)
    maes = vec.MAE_van(:);
    idx = find(isfinite(maes));
    [~, ord] = sort(maes(idx));
    n = numel(idx);
    if n == 0
        cases.idx = [];
        return;
    end
    pos = unique([round(0.25*n), round(0.5*n)-2, round(0.75*n)+1]);
    cases.idx = idx(ord(pos));
    cases.labels = {'25th', 'Median', '75th'};
end
