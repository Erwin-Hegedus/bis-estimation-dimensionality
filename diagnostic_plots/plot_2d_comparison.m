function plot_2d_comparison(results, patient_id, cfg, fig_dir)

    idx = find(results.patient_id == patient_id, 1);
    if isempty(idx)
        warning('Patient %d not found', patient_id);
        return;
    end
    
    r = results.raw(idx);
    t = r.time / 60;
    
    figure('Name', sprintf('2D Model Analysis - Case %d', patient_id), ...
           'Color', 'w', 'Position', [100 100 1400 900]);
    
    % Panel 1: BIS Predictions
    subplot(3, 2, [1 2]);
    plot(t, r.bis, 'k-', 'LineWidth', 1.0, 'DisplayName', 'Measured');
    hold on;
    if isfield(r, 'pred_pop') && ~isempty(r.pred_pop)
        plot(t, r.pred_pop, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Pop (0D)');
    end
    if isfield(r, 'pred_kscale') && ~isempty(r.pred_kscale)
        plot(t, r.pred_kscale, 'c-', 'LineWidth', 1.5, 'DisplayName', 'K (1D)');
    end
    if isfield(r, 'pred_2d') && ~isempty(r.pred_2d)
        plot(t, r.pred_2d, 'Color', [1 0.5 0], 'LineWidth', 2.0, 'DisplayName', '2D (kP,kR)');
    end
    plot(t, r.pred_van, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Van (4D)');
    ylabel('BIS'); xlabel('Time (min)');
    legend('Location', 'best');
    title(sprintf('Case %d: 1D vs 2D vs 4D', patient_id), 'FontSize', 14, 'FontWeight', 'bold');
    grid on; ylim([0 100]);
    
    % MAE annotations
    mae_pop = mean(abs(r.pred_pop - r.bis), 'omitnan');
    mae_k = mean(abs(r.pred_kscale - r.bis), 'omitnan');
    mae_2d = mean(abs(r.pred_2d - r.bis), 'omitnan');
    mae_van = mean(abs(r.pred_van - r.bis), 'omitnan');
    
    text(0.02, 0.98, sprintf('MAE:\nPop: %.2f\nK(1D): %.2f\n2D: %.2f\nVan(4D): %.2f', ...
        mae_pop, mae_k, mae_2d, mae_van), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontSize', 10);
    
    % Panel 2: CeP, CeR
    subplot(3, 2, 3);
    yyaxis left;
    plot(t, r.CeP_trajectory, 'b-', 'LineWidth', 1.5);
    ylabel('CeP (\mug/ml)'); ylim([0 8]);
    yyaxis right;
    plot(t, r.CeR_trajectory * 1000, 'r-', 'LineWidth', 1.5);
    ylabel('CeR (ng/ml)'); ylim([0 15]);
    title('Effect-site Concentrations'); grid on;
    xlabel('Time (min)');
    
    % Panel 3: k_scale (1D) vs (kP, kR) comparison
    subplot(3, 2, 4);
    hold on;
    if isfield(r, 'k_trajectory') && ~isempty(r.k_trajectory)
        plot(t, r.k_trajectory, 'c-', 'LineWidth', 2, 'DisplayName', 'k (1D)');
    end
    if isfield(r, 'kP_trajectory') && ~isempty(r.kP_trajectory)
        plot(t, r.kP_trajectory, 'b-', 'LineWidth', 2, 'DisplayName', 'k_P (2D)');
        plot(t, r.kR_trajectory, 'r-', 'LineWidth', 2, 'DisplayName', 'k_R (2D)');
    end
    yline(1.0, 'k--', 'Population');
    ylabel('Potency Scale'); ylim([0.3 3]);
    xlabel('Time (min)');
    legend('Location', 'best');
    title('1D vs 2D Parameters'); grid on;
    
    % Panel 4: kP vs kR scatter (identifiability diagnostic)
    subplot(3, 2, 5);
    if isfield(r, 'kP_trajectory') && ~isempty(r.kP_trajectory)
        kP = r.kP_trajectory(~isnan(r.kP_trajectory));
        kR = r.kR_trajectory(~isnan(r.kR_trajectory));
        n_pts = min(length(kP), length(kR));
        if n_pts > 10
            scatter(kP(1:n_pts), kR(1:n_pts), 20, (1:n_pts)', 'filled');
            colormap(gca, jet);
            colorbar('Ticks', [1 n_pts], 'TickLabels', {'Start', 'End'});
            hold on;
            plot([0.3 3], [0.3 3], 'k--', 'LineWidth', 1.5);
            xlabel('k_P'); ylabel('k_R');
            title('kP vs kR Trajectory (color = time)');
            xlim([0.3 3]); ylim([0.3 3]);
            grid on;
            axis equal;
            
            % Compute correlation
            corr_kP_kR = corrcoef(kP(1:n_pts), kR(1:n_pts));
            text(0.98, 0.02, sprintf('Corr: %.2f', corr_kP_kR(1,2)), ...
                'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'bottom', 'BackgroundColor', 'w');
        end
    end
    
    % Panel 5: Variance evolution
    subplot(3, 2, 6);
    hold on;
    if isfield(r, 'k_var') && ~isempty(r.k_var)
        t_k = (1:length(r.k_var))' / 60;
        semilogy(t_k, r.k_var, 'c-', 'LineWidth', 1.5, 'DisplayName', 'Var(k) 1D');
    end
    if isfield(r, 'kP_var') && ~isempty(r.kP_var)
        t_2d = (1:length(r.kP_var))' / 60;
        semilogy(t_2d, r.kP_var, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Var(k_P)');
        semilogy(t_2d, r.kR_var, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Var(k_R)');
    end
    legend('Location', 'best');
    ylabel('Variance'); xlabel('Time (min)');
    title('Parameter Uncertainty Evolution'); grid on;
    
    sgtitle(sprintf('2D Model (kP, kR) Diagnostic - Patient %d\nQuestion: Does separate P/R sensitivity improve fit?', patient_id), ...
        'FontSize', 12);
    
    print(gcf, fullfile(fig_dir, sprintf('2d_analysis_patient_%d.png', patient_id)), '-dpng', '-r150');
end
