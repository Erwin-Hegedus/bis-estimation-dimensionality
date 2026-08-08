function generate_figure3_dimensionality(results, fig_dir)
% FIGURE 3: In-sample MAE against model dimension, all models at process noise cfg.q.

    fprintf('Generating Figure 3: MAE vs model dimension...\n');
    
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
    
    N = sum(valid_mask);

    fig = figure('Name', 'Figure 3: MAE vs model dimension', 'Color', 'w', ...
                 'Position', [50 50 900 600]);
    
    % Main plot - just Bouillon line for clarity
    dims_main = [0, 1, 2, 3, 4];
    means_main = [mean(mae_pop), mean(mae_k), mean(mae_2d), mean(mae_ll), mean(mae_van)];
    stds_main = [std(mae_pop), std(mae_k), std(mae_2d), std(mae_ll), std(mae_van)];
    
    errorbar(dims_main, means_main, stds_main, 'ko-', 'LineWidth', 1.4, ...
             'MarkerSize', 6, 'MarkerFaceColor', 'b', 'CapSize', 5);
    hold on;
    
    % Add Greco as separate point
    errorbar(4.15, mean(mae_gre), std(mae_gre), 'rs', 'LineWidth', 1.4, ...
             'MarkerSize', 5, 'MarkerFaceColor', 'r', 'CapSize', 4);
    
    % Add reference line at 1D
    yline(mean(mae_k), 'c--', 'LineWidth', 2);
    text(3.5, mean(mae_k) - 0.3, '1D baseline', 'Color', 'c', 'FontWeight', 'bold');
    
    % Annotations
    for ii = 1:5
        text(dims_main(ii), means_main(ii) + stds_main(ii) + 0.5, ...
             sprintf('%.2f', means_main(ii)), ...
             'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    text(4.15, mean(mae_gre) + std(mae_gre) + 0.5, sprintf('%.2f', mean(mae_gre)), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'r');
    
    xlabel('Model dimension (estimated parameters)');
    ylabel('MAE ± SD (BIS units)');
    title(sprintf('In-sample MAE vs model dimension (N=%d)', N));

    legend({'Bouillon 4D', 'Greco 4D', '1D reference'}, ...
           'Location', 'northeast', 'Box', 'off');
    
    xlim([-0.5 4.7]);
    lo = min([means_main - stds_main, mean(mae_gre) - std(mae_gre)]);
    hi = max([means_main + stds_main, mean(mae_gre) + std(mae_gre)]);
    ylim([max(0, lo - 0.5), hi + 1.5]);
    grid on;
    set(gca, 'XTick', 0:4);
    
    d12 = mean(mae_2d) - mean(mae_k);
    [~, best] = min([mean(mae_pop), mean(mae_k), mean(mae_2d), mean(mae_ll), mean(mae_van)]);
    labels = {'population', '1D', '2D', '3D', '4D'};
    [~, p12] = ttest(mae_2d, mae_k);
    annotation('textbox', [0.17 0.14 0.58 0.15], ...
               'String', {sprintf('MAE is not monotone in dimension: %s is lowest.', labels{best}), ...
                          sprintf('1D\\rightarrow2D: %+.2f BIS (paired t, p = %.2f).', d12, p12), ...
                          'All models at the same process noise q.'}, ...
               'BackgroundColor', 'w', 'EdgeColor', 'k');

    save_figure(fig, fig_dir, 'figure3_dimensionality', 'single', 3.2);
end
