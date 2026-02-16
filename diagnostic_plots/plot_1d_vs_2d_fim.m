function plot_1d_vs_2d_fim(results, patient_id, fig_dir)
% PLOT_1D_VS_2D_FIM_REPRESENTATIVE - Clean comparison of 1D vs 2D (FIM)
    
    idx = find(results.patient_id == patient_id, 1);
    if isempty(idx)
        warning('Patient %d not found', patient_id);
        return;
    end
    
    r = results.raw(idx);
    
    % Cut last 15% (emergence)
    n_total = length(r.time);
    n_keep = floor(0.85 * n_total);
    
    t = r.time(1:n_keep) / 60;
    bis = r.bis(1:n_keep);
    pred_1d = r.pred_kscale(1:n_keep);
    pred_2d = r.pred_2d_fim(1:n_keep);
    CeP = r.CeP_trajectory(1:n_keep);
    CeR = r.CeR_trajectory(1:n_keep);
    k_1d = r.k_trajectory(1:n_keep);
    kP_2d = r.kP_fim_trajectory(1:n_keep);
    kR_2d = r.kR_fim_trajectory(1:n_keep);
    
    figure('Name', sprintf('1D vs 2D Comparison - Case %d', patient_id), ...
           'Color', 'w', 'Position', [50 50 1400 900]);
    
    % === Panel 1: BIS Trajectories ===
    subplot(3, 2, [1 2]);
    plot(t, bis, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Measured');
    hold on;
    plot(t, pred_1d, 'c-', 'LineWidth', 2, 'DisplayName', '1D (k)');
    plot(t, pred_2d, '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 2, 'DisplayName', '2D (k_P, k_R)');
    ylabel('BIS', 'FontSize', 12);
    xlabel('Time (min)', 'FontSize', 12);
    title(sprintf('Case %d: 1D vs 2D BIS Prediction', patient_id), 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    ylim([0 100]);
    
    % === Panel 2: Effect-site Concentrations ===
    subplot(3, 2, 3);
    yyaxis left;
    plot(t, CeP, 'b-', 'LineWidth', 1.5);
    ylabel('Ce_P (\mug/ml)', 'FontSize', 11);
    ylim([0 max(CeP)*1.1]);
    yyaxis right;
    plot(t, CeR*1000, 'r-', 'LineWidth', 1.5);
    ylabel('Ce_R (ng/ml)', 'FontSize', 11);
    ylim([0 max(CeR*1000)*1.1]);
    xlabel('Time (min)', 'FontSize', 11);
    title('Effect-site Concentrations', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % === Panel 3: 1D Parameter Evolution ===
    subplot(3, 2, 4);
    plot(t, k_1d, 'c-', 'LineWidth', 2);
    hold on;
    yline(1.0, 'k--', 'LineWidth', 1);
    ylabel('k (potency scale)', 'FontSize', 11);
    xlabel('Time (min)', 'FontSize', 11);
    title('1D Parameter: k', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    ylim([0.4 2.5]);
    
    % === Panel 4: 2D Parameter Evolution ===
    subplot(3, 2, 5);
    plot(t, kP_2d, 'b-', 'LineWidth', 2, 'DisplayName', 'k_P');
    hold on;
    plot(t, kR_2d, 'r-', 'LineWidth', 2, 'DisplayName', 'k_R');
    yline(1.0, 'k--', 'LineWidth', 1);
    ylabel('Potency Scale', 'FontSize', 11);
    xlabel('Time (min)', 'FontSize', 11);
    title('2D Parameters: k_P, k_R', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on;
    ylim([0.4 2.5]);
    
    % === Panel 5: kP vs kR Trajectory ===
    subplot(3, 2, 6);
    valid_idx = ~isnan(kP_2d) & ~isnan(kR_2d);
    kP_v = kP_2d(valid_idx);
    kR_v = kR_2d(valid_idx);
    t_v = t(valid_idx);
    scatter(kP_v, kR_v, 20, t_v, 'filled');
    colormap(gca, jet);
    cb = colorbar;
    cb.Label.String = 'Time (min)';
    hold on;
    plot(kP_v(1), kR_v(1), 'go', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Start');
    plot(kP_v(end), kR_v(end), 'rs', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'End');
    plot([0.4 2.5], [0.4 2.5], 'k--', 'LineWidth', 1);
    xlabel('k_P', 'FontSize', 11);
    ylabel('k_R', 'FontSize', 11);
    title('2D Parameter Trajectory', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on;
    axis equal;
    xlim([0.4 2.5]);
    ylim([0.4 2.5]);
    
    sgtitle(sprintf('1D vs 2D Model Comparison - Patient %d', patient_id), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    print(gcf, fullfile(fig_dir, sprintf('1d_vs_2d_representative_%d.png', patient_id)), '-dpng', '-r150');
    savefig(gcf, fullfile(fig_dir, sprintf('1d_vs_2d_representative_%d.fig', patient_id)));
end
