function plot_case_study_params_grid_v3(cases, results, cfg, fig_dir)
    m = numel(cases.idx);
    if m == 0, return; end
    
    names = {'C50P', 'C50R', 'Gamma', 'Beta/Alpha'};
    f = figure('Color', 'w', 'Position', [80 80 1200 200*m]);
    
    for r = 1:m
        i = cases.idx(r);
        
        % CUT LAST 15% (emergence phase)
        n_total = length(results.raw(i).time);
        n_keep = floor(0.80 * n_total);
        
        XV = results.raw(i).Xhist_van;
        XG = results.raw(i).Xhist_gre;
        t = results.raw(i).time(1:n_keep) / 60;
        
        for p = 1:4
            subplot(m, 4, (r-1)*4 + p);
            hold on;
            grid on;
            if ~isempty(XV)
                plot(t, XV(p, 1:n_keep), 'b-');
                xlim([0 max(t)]);
            end
            if ~isempty(XG)
                plot(t, XG(p, 1:n_keep), 'r--');
                xlim([0 max(t)]);
            end
            if r == 1
                title(names{p});
            end
            if p == 1
                ylabel(sprintf('Case %d', results.patient_id(i)));
            end
        end
    end
    sgtitle('4D Parameter Evolution');
    save_figure(f, fig_dir, 'figure_4d_params');
end
