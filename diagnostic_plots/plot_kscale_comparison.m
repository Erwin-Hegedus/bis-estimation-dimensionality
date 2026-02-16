function plot_kscale_comparison(results, patient_id)
    idx = find(results.patient_id == patient_id, 1);
    if isempty(idx)
        warning('Patient %d not found', patient_id);
        return;
    end
    
    r = results.raw(idx);
    t = r.time / 60;
    
    figure('Name', sprintf('k_scale Analysis - Case %d', patient_id), ...
           'Color', 'w', 'Position', [100 100 1400 900]);
    
    subplot(3, 2, [1 2]);
    plot(t, r.bis, 'k-', 'LineWidth', 1.0, 'DisplayName', 'Measured BIS');
    hold on;
    if ~isempty(r.pred_pop)
        plot(t, r.pred_pop, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Population');
    end
    plot(t, r.pred_van, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Vanluchene (4D)');
    if ~isempty(r.pred_kscale)
        plot(t, r.pred_kscale, 'r-', 'LineWidth', 2.0, 'DisplayName', 'k\_scale (1D)');
    end
    if isfield(r, 'pred_2d') && ~isempty(r.pred_2d)
        plot(t, r.pred_2d, 'Color', [1 0.5 0], 'LineWidth', 2.0, 'DisplayName', '2D');
    end
    ylabel('BIS');
    xlabel('Time (min)');
    legend('Location', 'best');
    title(sprintf('Case %d: 1-Parameter vs 2-Parameter vs 4-Parameter Model Performance', patient_id), ...
          'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    ylim([0 100]);
    
    mae_pop = mean(abs(r.pred_pop - r.bis), 'omitnan');
    mae_van = mean(abs(r.pred_van - r.bis), 'omitnan');
    mae_k = mean(abs(r.pred_kscale - r.bis), 'omitnan');
    mae_2d = mean(abs(r.pred_2d - r.bis), 'omitnan');
    
    text(0.02, 0.98, sprintf('MAE Pop: %.2f\nMAE Van(4D): %.2f\nMAE k_{scale}(1D): %.2f\nMAE 2D: %.2f', ...
        mae_pop, mae_van, mae_k, mae_2d), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontSize', 11);
    
    subplot(3, 2, 3);
    if isfield(r, 'k_trajectory') && ~isempty(r.k_trajectory)
        plot(t, r.k_trajectory, 'r-', 'LineWidth', 2);
        hold on;
        yline(1.0, 'k--', 'Population (k=1)');
        ylabel('k_{scale}');
        xlabel('Time (min)');
        title('Potency Scale Factor Evolution');
        grid on;
        ylim([0.3 3.0]);
        
        C50P_pop = 3.5;
        C50R_pop = 5.0;
        k_final = r.k_trajectory(end);
        text(0.98, 0.02, sprintf('Final k=%.2f\nEquiv C50P=%.1f\nEquiv C50R=%.1f', ...
            k_final, C50P_pop/k_final, C50R_pop/k_final), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'bottom', 'BackgroundColor', 'w', 'EdgeColor', 'k');
    end
    
    subplot(3, 2, 4);
    if isfield(r, 'k_var') && ~isempty(r.k_var)
        t_k = (1:length(r.k_var))' / 60;
        semilogy(t_k, r.k_var, 'r-', 'LineWidth', 2);
        ylabel('Var(k)');
        xlabel('Time (min)');
        title('k_{scale} Variance (Uncertainty)');
        grid on;
    end
    
    subplot(3, 2, 5);
    if ~isempty(r.Xhist_van)
        X_norm = r.Xhist_van ./ r.Xhist_van(:,1);
        plot(t, X_norm(1,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'C50P');
        hold on;
        plot(t, X_norm(2,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'C50R');
        plot(t, X_norm(3,:), 'g-', 'LineWidth', 1.5, 'DisplayName', '\gamma');
        plot(t, X_norm(4,:), 'm-', 'LineWidth', 1.5, 'DisplayName', '\beta');
        yline(1.0, 'k--');
        ylabel('Relative to Initial');
        xlabel('Time (min)');
        title('4D Parameters (Normalized)');
        legend('Location', 'best');
        grid on;
    end
    
    subplot(3, 2, 6);
    err_van = r.pred_van - r.bis;
    err_k = r.pred_kscale - r.bis;
    err_2d = r.pred_2d - r.bis;
    
    plot(t, err_van, 'b-', 'LineWidth', 1.0, 'DisplayName', 'Van(4D) Error');
    hold on;
    plot(t, err_k, 'r-', 'LineWidth', 1.0, 'DisplayName', 'k_{scale}(1D) Error');
    plot(t, err_2d, 'Color', [1 0.5 0], 'LineWidth', 1.0, 'DisplayName', '2D Error');
    yline(0, 'k-');
    ylabel('Prediction Error (BIS)');
    xlabel('Time (min)');
    title('Prediction Errors: 1D vs 2D vs 4D');
    legend('Location', 'best');
    grid on;
    ylim([-30 30]);
    
    sgtitle(sprintf('IDENTIFIABILITY ANALYSIS: Patient %d\nComparing 1D, 2D, and 4D Models', patient_id), ...
            'FontSize', 12, 'FontWeight', 'bold');
end
