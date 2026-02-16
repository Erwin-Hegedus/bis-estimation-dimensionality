function plot_case_studies_v3(cases, results, cfg, fig_dir)
    m = numel(cases.idx);
    if m == 0, return; end
    
    f = figure('Color', 'w', 'Position', [50, 50, 1600, 300*m]);
    for k = 1:m
        i = cases.idx(k);
        
        % CUT LAST 15% (emergence phase)
        n_total = length(results.raw(i).time);
        n_keep = floor(0.80 * n_total);
        
        t = results.raw(i).time(1:n_keep) / 60;
        
        subplot(m, 2, 2*k - 1);
        CeP = results.raw(i).CeP_trajectory(1:n_keep);
        CeR = results.raw(i).CeR_trajectory(1:n_keep);
        yyaxis left;
        plot(t, CeP, 'b-', 'LineWidth', 1.5);
        ylabel('Prop (\mug/ml)');
        set(gca, 'YColor', 'b');
        ylim([0 8]);
        xlim([0 max(t)]);
        yyaxis right;
        plot(t, CeR * 1000, 'r-', 'LineWidth', 1.5);
        ylabel('Remi (ng/ml)');
        set(gca, 'YColor', [0.8 0 0]);
        ylim([0 15]);
        xlim([0 max(t)]);
        title(sprintf('Case %d Concentrations', results.patient_id(i)));
        grid on;
        
        subplot(m, 2, 2*k);
        % MEASURED BIS
        plot(t, results.raw(i).bis(1:n_keep), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Measured BIS');
        hold on;
        
        % ENDPOINTS
        if isfield(results.raw(i), 'E0_trajectory') && ~isempty(results.raw(i).E0_trajectory)
            plot(t, results.raw(i).E0_trajectory(1:n_keep), 'Color', [0.6 0.6 0.6], 'LineStyle', '-.', 'LineWidth', 1.2, 'DisplayName', 'E_0');
        end
        if isfield(results.raw(i), 'BISmin_trajectory') && ~isempty(results.raw(i).BISmin_trajectory)
            plot(t, results.raw(i).BISmin_trajectory(1:n_keep), 'Color', [0.4 0.4 0.4], 'LineStyle', '-.', 'LineWidth', 1.2, 'DisplayName', 'BIS_{min}');
        end
        
        % 4D MODEL PREDICTIONS (Van and Greco)
        plot(t, results.raw(i).pred_van(1:n_keep), 'b-', 'LineWidth', 2.0, 'DisplayName', '4D Bouillon');
        plot(t, results.raw(i).pred_gre(1:n_keep), 'r--', 'LineWidth', 2.0, 'DisplayName', '4D Greco');
        
        % TARGET BAND
        yline(cfg.target_band(1), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
        yline(cfg.target_band(2), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
        
        ylim([0 100]);
        xlim([0 max(t)]);
        ylabel('BIS');
        grid on;
        
        % IMPROVED TITLE with percentile and metrics
        mae_van = results.metrics.van.MAE(i);
        mae_gre = results.metrics.gre.MAE(i);
        title(sprintf('%s Percentile - Case %d\nMAE: Bouillon=%.1f | Greco=%.1f', ...
            cases.labels{k}, results.patient_id(i), mae_van, mae_gre), ...
            'FontSize', 10, 'FontWeight', 'bold');
        
        if k == 1
            legend('Location', 'best', 'FontSize', 9);
        end
    end
    sgtitle('4D Model Predictions', 'FontSize', 13, 'FontWeight', 'bold');
    save_figure(f, fig_dir, 'figure_cases_comparison');
end
