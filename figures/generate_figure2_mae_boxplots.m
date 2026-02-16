function generate_figure2_mae_boxplots(results, fig_dir)
% FIGURE 2: Distribution of Prediction Errors Across Model Complexities
% Clean, professional style with subtle color distinctions

    fprintf('Generating Figure 2: MAE Boxplots...\n');
    
    % Get valid data
    valid_mask = ~isnan(results.metrics.van.MAE) & ...
                 ~isnan(results.metrics.gre.MAE) & ...
                 ~isnan(results.metrics.pop.MAE) & ...
                 ~isnan(results.metrics.kscale.MAE) & ...
                 ~isnan(results.metrics.m2d.MAE) & ...
                 ~isnan(results.metrics.loglin.MAE);
    
    mae_pop = results.metrics.pop.MAE(valid_mask);
    mae_k = results.metrics.kscale.MAE(valid_mask);
    mae_2d = results.metrics.m2d.MAE(valid_mask);
    mae_ll = results.metrics.loglin.MAE(valid_mask);
    mae_van = results.metrics.van.MAE(valid_mask);
    mae_gre = results.metrics.gre.MAE(valid_mask);
    N = length(mae_pop);
    
    fig = figure('Name', 'Figure 2: MAE Boxplots', 'Color', 'w', ...
                 'Position', [50 50 1000 500]);
    
    % Combine data for boxplot - use simple labels first, fix later
    data = [mae_pop; mae_k; mae_2d; mae_ll; mae_van; mae_gre];
    groups = [ones(N,1); 2*ones(N,1); 3*ones(N,1); 4*ones(N,1); 5*ones(N,1); 6*ones(N,1)];
    
    % Create boxplot
    boxplot(data, groups, 'Widths', 0.5);
    
    % Professional muted palette - like the example but with subtle variation
    colors = [0.75 0.75 0.75;    % 0D - gray
              0.70 0.80 0.92;    % 1D - light steel blue
              0.85 0.75 0.70;    % 2D - light taupe
              0.80 0.75 0.85;    % 3D - light lavender
              0.70 0.78 0.88;    % 4D Van - light blue
              0.88 0.75 0.75];   % 4D Gre - light rose
    
    edge_colors = [0.50 0.50 0.50;    % 0D
                   0.40 0.55 0.70;    % 1D
                   0.60 0.50 0.45;    % 2D
                   0.55 0.50 0.60;    % 3D
                   0.40 0.50 0.65;    % 4D Van
                   0.65 0.45 0.45];   % 4D Gre
    
    boxes = findobj(gca, 'Tag', 'Box');
    for jj = 1:length(boxes)
        idx = length(boxes) - jj + 1;
        patch(get(boxes(jj), 'XData'), get(boxes(jj), 'YData'), ...
              colors(idx, :), 'FaceAlpha', 0.5, 'EdgeColor', edge_colors(idx, :));
    end
    
    hold on;
    
    % Add mean markers
    means = [mean(mae_pop), mean(mae_k), mean(mae_2d), mean(mae_ll), mean(mae_van), mean(mae_gre)];
    plot(1:6, means, 'o', 'MarkerSize', 6, 'MarkerFaceColor', [0.3 0.3 0.3], ...
         'MarkerEdgeColor', 'none');
    
    % Add mean labels
    for ii = 1:6
        text(ii, means(ii) + 0.7, sprintf('%.2f', means(ii)), ...
             'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);
    end
    
    ylabel('Mean Absolute Error (BIS units)', 'FontSize', 11);
    xlabel('Model Complexity', 'FontSize', 11);
    title('Prediction Error Distribution', 'FontSize', 12, 'FontWeight', 'normal');
    
    % Fix x-axis labels with proper subscripts
    set(gca, 'XTick', 1:6);
    set(gca, 'XTickLabel', {'0D (Pop)', '1D (k)', '2D (k_P, k_R)', '3D (LogLin)', '4D (Bouillon)', '4D (Greco)'});
    set(gca, 'TickLabelInterpreter', 'tex');
    
    % Clean axis styling
    set(gca, 'FontSize', 10, 'Box', 'on', 'LineWidth', 0.5);
    ylim([0 max(data) * 1.1]);
    
    % Subtle reference line at 1D performance
    yline(mean(mae_k), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    
    save_figure(fig, fig_dir, 'figure2_mae_boxplots');
end
