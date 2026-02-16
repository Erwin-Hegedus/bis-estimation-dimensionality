function generate_figure6_2d_drift(results, fig_dir)
% FIGURE 6: Correlated Parameter Drift in the 2D Model
% kP vs kR trajectories colored by time for multiple cases

    fprintf('Generating Figure 6: 2D Parameter Drift...\n');
    
    % Select representative cases
    valid_idx = find(~isnan(results.metrics.m2d.MAE));
    if length(valid_idx) < 4
        warning('Not enough valid cases for Figure 6');
        return;
    end
    
    % Select 4 diverse cases
    n_valid = length(valid_idx);
    case_indices = valid_idx([1, round(n_valid/3), round(2*n_valid/3), n_valid]);
    
    fig = figure('Name', 'Figure 6: 2D Parameter Drift', 'Color', 'w', ...
                 'Position', [50 50 1200 500]);
    
    % Single large panel with all cases
    hold on;
    colors = {'b', 'r', 'g', 'm'};
    
    for ii = 1:length(case_indices)
        idx = case_indices(ii);
        r = results.raw(idx);
        
        if isempty(r.kP_trajectory) || isempty(r.kR_trajectory)
            continue;
        end
        
        kP = r.kP_trajectory(~isnan(r.kP_trajectory));
        kR = r.kR_trajectory(~isnan(r.kR_trajectory));
        n_pts = min(length(kP), length(kR));
        
        if n_pts < 100
            continue;
        end
        
        % Subsample for clarity
        step = max(1, floor(n_pts / 200));
        idx_plot = 1:step:n_pts;
        
        scatter(kP(idx_plot), kR(idx_plot), 20, (1:length(idx_plot))', ...
                'filled', 'MarkerFaceAlpha', 0.6);
    end
    
    % Add diagonal line (kP = kR)
    plot([0.3 2.5], [0.3 2.5], 'k--', 'LineWidth', 2);
    
    % Add population reference
    scatter(1, 1, 200, 'k', 'p', 'filled', 'DisplayName', 'Population');
    
    xlabel('k_P (Propofol potency scale)', 'FontSize', 12);
    ylabel('k_R (Remifentanil potency scale)', 'FontSize', 12);
    title('Figure 6: Correlated Parameter Drift in 2D Model', 'FontSize', 14, 'FontWeight', 'bold');
    
    colormap(jet);
    cb = colorbar;
    cb.Label.String = 'Time progression';
    cb.Label.FontSize = 11;
    
    xlim([0.3 2.5]);
    ylim([0.3 2.5]);
    axis square;
    grid on;
    
    % Add annotation
    annotation('textbox', [0.15 0.72 0.25 0.1], ...
               'String', {'Parameters drift along', 'constant-potency ridge'}, ...
               'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontSize', 10);
    
    save_figure(fig, fig_dir, 'figure6_2d_drift');
end
