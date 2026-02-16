function plot_1d_vs_3d(results, patient_id, fig_dir)
% PLOT_1D_VS_3D_REPRESENTATIVE - Clean comparison of 1D vs 3D (LogLin)
    
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
    pred_3d = r.pred_loglin(1:n_keep);
    CeP = r.CeP_trajectory(1:n_keep);
    CeR = r.CeR_trajectory(1:n_keep);
    k_1d = r.k_trajectory(1:n_keep);
    
    % Extract 3D parameters from Xhist_loglin
    if isfield(r, 'Xhist_loglin') && ~isempty(r.Xhist_loglin)
        a0_3d = r.Xhist_loglin(1, 1:n_keep)';
        aP_3d = r.Xhist_loglin(2, 1:n_keep)';
        aR_3d = r.Xhist_loglin(3, 1:n_keep)';
    else
        warning('No 3D parameter history available');
        a0_3d = nan(n_keep, 1);
        aP_3d = nan(n_keep, 1);
        aR_3d = nan(n_keep, 1);
    end
    
    figure('Name', sprintf('1D vs 3D Comparison - Case %d', patient_id), ...
           'Color', 'w', 'Position', [50 50 1500 900]);
    
    % === Panel 1: BIS Trajectories ===
    subplot(3, 2, [1 2]);
    plot(t, bis, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Measured');
    hold on;
    plot(t, pred_1d, 'c-', 'LineWidth', 2, 'DisplayName', '1D (k)');
    plot(t, pred_3d, 'm-', 'LineWidth', 2, 'DisplayName', '3D (a_0, a_P, a_R)');
    ylabel('BIS', 'FontSize', 12);
    xlabel('Time (min)', 'FontSize', 12);
    title(sprintf('Case %d: 1D vs 3D BIS Prediction', patient_id), 'FontSize', 14, 'FontWeight', 'bold');
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
    
    % === Panel 4: 3D Parameter Evolution (2D time plot) ===
    subplot(3, 2, 5);
    yyaxis left;
    plot(t, a0_3d, 'b-', 'LineWidth', 2, 'DisplayName', 'a_0');
    ylabel('a_0 (intercept)', 'FontSize', 11);
    yyaxis right;
    plot(t, aP_3d, '-', 'Color', [0 0.6 0], 'LineWidth', 2, 'DisplayName', 'a_P');
    hold on;
    plot(t, aR_3d, 'r-', 'LineWidth', 2, 'DisplayName', 'a_R');
    ylabel('a_P, a_R (slopes)', 'FontSize', 11);
    xlabel('Time (min)', 'FontSize', 11);
    title('3D Parameters vs Time', 'FontSize', 12, 'FontWeight', 'bold');
    legend({'a_0', 'a_P', 'a_R'}, 'Location', 'best');
    grid on;
    
    % === Panel 5: 3D Parameter Trajectory in 3D Space ===
    subplot(3, 2, 6);
    valid_idx = ~isnan(a0_3d) & ~isnan(aP_3d) & ~isnan(aR_3d);
    a0_v = a0_3d(valid_idx);
    aP_v = aP_3d(valid_idx);
    aR_v = aR_3d(valid_idx);
    t_v = t(valid_idx);
    
    if ~isempty(a0_v) && length(a0_v) > 10
        % 3D scatter with time color
        scatter3(a0_v, aP_v, aR_v, 25, t_v, 'filled');
        colormap(gca, jet);
        cb = colorbar;
        cb.Label.String = 'Time (min)';
        hold on;
        
        % Plot trajectory line
        plot3(a0_v, aP_v, aR_v, 'k-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5]);
        
        % Mark start and end
        plot3(a0_v(1), aP_v(1), aR_v(1), 'go', 'MarkerSize', 14, 'LineWidth', 3);
        plot3(a0_v(end), aP_v(end), aR_v(end), 'rs', 'MarkerSize', 14, 'LineWidth', 3);
        
        xlabel('a_0', 'FontSize', 11);
        ylabel('a_P', 'FontSize', 11);
        zlabel('a_R', 'FontSize', 11);
        title('3D Parameter Trajectory', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        view(45, 30);
    else
        text(0.5, 0.5, 0.5, 'Insufficient 3D data', 'HorizontalAlignment', 'center');
    end
    
    sgtitle(sprintf('1D vs 3D (LogLin) Model Comparison - Patient %d', patient_id), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    print(gcf, fullfile(fig_dir, sprintf('1d_vs_3d_representative_%d.png', patient_id)), '-dpng', '-r150');
    savefig(gcf, fullfile(fig_dir, sprintf('1d_vs_3d_representative_%d.fig', patient_id)));
end
