function generate_figure5_monte_carlo(eq_info, cfg, fig_dir)
% FIGURE 5: Visualization of Practical Non-Identifiability (Monte Carlo)
% 2D projections onto C50P-C50R and gamma-beta planes

    fprintf('Generating Figure 5: Monte Carlo Equivalence...\n');
    
    fig = figure('Name', 'Figure 5: Monte Carlo Equivalence', 'Color', 'w', ...
                 'Position', [50 50 1200 500]);
    
    if ~eq_info.has_eq || isempty(eq_info.theta_good)
        % No equivalent sets found - show message
        text(0.5, 0.5, 'No equivalent parameter sets found within tolerance', ...
             'HorizontalAlignment', 'center', 'FontSize', 14);
        save_figure(fig, fig_dir, 'figure5_monte_carlo');
        return;
    end
    
    theta_good = eq_info.theta_good;  % N x 4 matrix
    theta_ref = eq_info.theta_ref;     % 1 x 4 vector
    
    % Panel (a): C50P vs C50R (Potency Plane)
    subplot(1, 2, 1);
    scatter(theta_good(:,1), theta_good(:,2), 26, [0.3 0.5 0.8], 'filled', ...
            'MarkerFaceAlpha', 0.5, 'DisplayName', 'Equivalent sets');
    hold on;
    scatter(theta_ref(1), theta_ref(2), 90, 'r', 'x', 'LineWidth', 1.8, ...
            'DisplayName', 'Reference estimate');
    
    % Add bounds
    rectangle('Position', [cfg.lb(1), cfg.lb(2), cfg.ub(1)-cfg.lb(1), cfg.ub(2)-cfg.lb(2)], ...
              'EdgeColor', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1);
    
    xlabel('C_{50P} (\mug/mL)');
    ylabel('C_{50R} (ng/mL)');
    title('(a) Potency plane');
    % One legend for the figure: both panels show the same two things.
    legend('Location', 'northwest', 'Box', 'off');
    grid on;
    xlim([cfg.lb(1)-0.5, cfg.ub(1)+0.5]);
    ylim([cfg.lb(2)-2, cfg.ub(2)+2]);
    
    % Panel (b): gamma vs beta (Shape Plane)
    subplot(1, 2, 2);
    scatter(theta_good(:,3), theta_good(:,4), 26, [0.3 0.5 0.8], 'filled', ...
            'MarkerFaceAlpha', 0.5);
    hold on;
    scatter(theta_ref(3), theta_ref(4), 90, 'r', 'x', 'LineWidth', 1.8);
    
    % Add bounds
    rectangle('Position', [cfg.lb(3), cfg.lb(4), cfg.ub(3)-cfg.lb(3), cfg.ub(4)-cfg.lb(4)], ...
              'EdgeColor', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1);
    
    xlabel('\gamma (Hill coefficient)');
    ylabel('\beta (interaction)');
    title('(b) Shape plane');
    grid on;
    xlim([cfg.lb(3)-0.2, cfg.ub(3)+0.2]);
    ylim([cfg.lb(4)-0.1, cfg.ub(4)+0.1]);
    
    % No sgtitle: N_eq and MAE_ref belong in the caption, and a two-line
    % superscripted title collided with the panel titles below it.
    save_figure(fig, fig_dir, 'figure5_monte_carlo', 'double', 2.6);
end
