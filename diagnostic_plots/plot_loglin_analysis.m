function plot_loglin_analysis(results, patient_id)
    idx = find(results.patient_id == patient_id, 1);
    if isempty(idx)
        warning('Patient %d not found', patient_id);
        return;
    end
    
    r = results.raw(idx);
    t = r.time / 60;
    
    figure('Name', sprintf('LogLin Analysis - Case %d', patient_id), ...
           'Color', 'w', 'Position', [100 100 1200 800]);
    
    subplot(2, 2, [1 2]);
    plot(t, r.bis, 'k-', 'LineWidth', 1.0, 'DisplayName', 'Measured');
    hold on;
    if ~isempty(r.pred_pop)
        plot(t, r.pred_pop, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Pop (0D)');
    end
    if ~isempty(r.pred_kscale)
        plot(t, r.pred_kscale, 'c-', 'LineWidth', 1.5, 'DisplayName', 'K (1D)');
    end
    if isfield(r, 'pred_2d') && ~isempty(r.pred_2d)
        plot(t, r.pred_2d, 'Color', [1 0.5 0], 'LineWidth', 1.5, 'DisplayName', '2D');
    end
    if ~isempty(r.pred_loglin)
        plot(t, r.pred_loglin, 'm-', 'LineWidth', 2.0, 'DisplayName', 'LL (3D)');
    end
    plot(t, r.pred_van, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Van (4D)');
    ylabel('BIS');
    xlabel('Time (min)');
    legend('Location', 'best');
    title(sprintf('Case %d: Model Comparison', patient_id), 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    ylim([0 100]);
    
    mae_pop = mean(abs(r.pred_pop - r.bis), 'omitnan');
    mae_k = mean(abs(r.pred_kscale - r.bis), 'omitnan');
    mae_2d = mean(abs(r.pred_2d - r.bis), 'omitnan');
    mae_ll = mean(abs(r.pred_loglin - r.bis), 'omitnan');
    mae_van = mean(abs(r.pred_van - r.bis), 'omitnan');
    
    text(0.02, 0.98, sprintf('MAE:\nPop: %.2f\nK(1D): %.2f\n2D: %.2f\nLL(3D): %.2f\nVan(4D): %.2f', ...
        mae_pop, mae_k, mae_2d, mae_ll, mae_van), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontSize', 10);
    
    subplot(2, 2, 3);
    if isfield(r, 'Xhist_loglin') && ~isempty(r.Xhist_loglin)
        plot(t, r.Xhist_loglin(1,:), 'r-', 'LineWidth', 2, 'DisplayName', 'a_0');
        hold on;
        plot(t, r.Xhist_loglin(2,:), 'b-', 'LineWidth', 2, 'DisplayName', 'a_P');
        plot(t, r.Xhist_loglin(3,:), 'g-', 'LineWidth', 2, 'DisplayName', 'a_R');
        ylabel('Parameter Value');
        xlabel('Time (min)');
        title('LogLin Parameter Evolution', 'FontWeight', 'bold');
        legend('Location', 'best');
        grid on;
    end
    
    subplot(2, 2, 4);
    if ~isempty(r.k_trajectory) && ~isempty(r.Xhist_van)
        yyaxis left;
        plot(t, r.k_trajectory, 'c-', 'LineWidth', 2);
        hold on;
        if isfield(r, 'kP_trajectory') && ~isempty(r.kP_trajectory)
            plot(t, r.kP_trajectory, 'b--', 'LineWidth', 1.5);
            plot(t, r.kR_trajectory, 'r--', 'LineWidth', 1.5);
        end
        ylabel('k_{scale} / k_P / k_R');
        ylim([0.3 3]);
        
        yyaxis right;
        C50_eff = r.Xhist_van(1,:) + r.Xhist_van(2,:)/1000;
        plot(t, C50_eff, 'b-', 'LineWidth', 2);
        ylabel('C50_P + C50_R/1000');
        
        xlabel('Time (min)');
        title('Potency Comparison', 'FontWeight', 'bold');
        legend('k (1D)', 'kP (2D)', 'kR (2D)', 'Location', 'best');
        grid on;
    end
end
