function generate_figure4_case_study(results, patient_id, cfg, fig_dir)
% FIGURE 4: Representative Case Study of Online Parameter Adaptation
% Multi-panel figure showing BIS, concentrations, and parameter evolution

    fprintf('Generating Figure 4: Case Study (Patient %d)...\n', patient_id);
    
    idx = find(results.patient_id == patient_id, 1);
    if isempty(idx)
        warning('Patient %d not found, using first valid patient', patient_id);
        idx = find(~isnan(results.metrics.van.MAE), 1);
        patient_id = results.patient_id(idx);
    end
    
    r = results.raw(idx);
    
    % CUT LAST 15% (emergence phase)
    n_total = length(r.time);
    n_keep = floor(0.83 * n_total);
    t = r.time(1:n_keep) / 60;
    
    fig = figure('Name', sprintf('Figure 4: Case Study (Patient %d)', patient_id), ...
                 'Color', 'w', 'Position', [50 50 1400 900]);
    
    % Panel (a): BIS predictions
    subplot(3, 2, [1 2]);
    plot(t, r.bis(1:n_keep), 'k-', 'LineWidth', 1.0, 'DisplayName', 'Measured BIS');
    hold on;
    plot(t, r.pred_pop(1:n_keep), 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', ':', 'DisplayName', '0D (Population)');
    plot(t, r.pred_kscale(1:n_keep), 'c-', 'LineWidth', 2.0, 'DisplayName', '1D (k)');
    plot(t, r.pred_2d(1:n_keep), 'b-', 'LineWidth', 1.5, 'DisplayName', '2D');
    
    % Target band
    yline(cfg.target_band(1), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1, 'DisplayName', 'Target Range (40-60)');
    yline(cfg.target_band(2), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1, 'HandleVisibility', 'off');
    
    xlim([0 max(t)]);
    ylabel('BIS', 'FontSize', 12);
    ylim([0 100]);
    legend('Location', 'best', 'FontSize', 10);
    title('(a) BIS Prediction Comparison', 'FontSize', 13, 'FontWeight', 'bold');
    grid on;
   
    
    % Panel (b): Effect-site concentrations
    subplot(3, 2, 3);
    yyaxis left;
    plot(t, r.CeP_trajectory(1:n_keep), 'b-', 'LineWidth', 1.5);
    ylabel('C_{eP} (\mug/mL)', 'FontSize', 11);
    ylim([0 8]);
    set(gca, 'YColor', 'b');
    
    yyaxis right;
    plot(t, r.CeR_trajectory(1:n_keep) * 1000, 'r-', 'LineWidth', 1.5);
    ylabel('C_{eR} (ng/mL)', 'FontSize', 11);
    ylim([0 15]);
    set(gca, 'YColor', [0.8 0 0]);
    
    xlim([0 max(t)]);
    xlabel('Time (min)', 'FontSize', 11);
    title('(b) Effect-Site Concentrations', 'FontSize', 13, 'FontWeight', 'bold');
    legend({'Propofol', 'Remifentanil'}, 'Location', 'best');
    grid on;
    
    % Panel (c): 1D k trajectory
    subplot(3, 2, 4);
    plot(t, r.k_trajectory(1:n_keep), 'c-', 'LineWidth', 2.5);
    hold on;
    yline(1.0, 'k--', 'LineWidth', 1.5);
    xlim([0 max(t)]);
    ylabel('k_{scale}', 'FontSize', 11);
    xlabel('Time (min)', 'FontSize', 11);
    title('(c) 1D Potency Factor Evolution', 'FontSize', 13, 'FontWeight', 'bold');
    ylim([0.3 2.5]);
    grid on;
    
    % Panel (d): 2D kP, kR comparison
    subplot(3, 2, 5);
    plot(t, r.kP_trajectory(1:n_keep), 'b-', 'LineWidth', 2, 'DisplayName', 'k_P');
    hold on;
    plot(t, r.kR_trajectory(1:n_keep), 'r-', 'LineWidth', 2, 'DisplayName', 'k_R');
    plot(t, r.k_trajectory(1:n_keep), 'c--', 'LineWidth', 1.5, 'DisplayName', 'k (1D)');
    yline(1.0, 'k:', 'LineWidth', 1);
    
    xlim([0 max(t)]);
    ylabel('Potency Scale', 'FontSize', 11);
    xlabel('Time (min)', 'FontSize', 11);
    title('(d) 2D vs 1D Parameter Comparison', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'best');
    ylim([0.3 2.5]);
    grid on;
    
    % Panel (e): Variance evolution
    subplot(3, 2, 6);
    if ~isempty(r.k_var)
        n_keep_var = min(n_keep, length(r.k_var));
        t_var = (1:n_keep_var)' / 60;
        semilogy(t_var, r.k_var(1:n_keep_var), 'c-', 'LineWidth', 2, 'DisplayName', 'Var(k) 1D');
        hold on;
    end
    if ~isempty(r.kP_var)
        n_keep_2d = min(n_keep, length(r.kP_var));
        t_2d = (1:n_keep_2d)' / 60;
        semilogy(t_2d, r.kP_var(1:n_keep_2d), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Var(k_P)');
        semilogy(t_2d, r.kR_var(1:n_keep_2d), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Var(k_R)');
    end
    
    xlim([0 max(t)]);
    ylabel('Parameter Variance (log scale)', 'FontSize', 11);
    xlabel('Time (min)', 'FontSize', 11);
    title('(e) Uncertainty Evolution', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on;
    
    save_figure(fig, fig_dir, 'figure4_case_study');
end
