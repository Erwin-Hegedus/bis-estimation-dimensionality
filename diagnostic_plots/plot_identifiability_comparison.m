function plot_identifiability_comparison(results, fig_dir)
% PLOT_IDENTIFIABILITY_COMPARISON - Compare 1D, 2D, 3D, 4D models
%
% KEY FIGURE showing model dimension vs MAE

    mae_pop = results.metrics.pop.MAE;
    mae_van = results.metrics.van.MAE;
    mae_gre = results.metrics.gre.MAE;
    mae_k = results.metrics.kscale.MAE;
    mae_2d = results.metrics.m2d.MAE;
    mae_ll = results.metrics.loglin.MAE;
    
    valid = isfinite(mae_pop) & isfinite(mae_van) & isfinite(mae_k) & ...
            isfinite(mae_2d) & isfinite(mae_ll) & isfinite(mae_gre);
    
    if sum(valid) < 3
        warning('Not enough valid cases for comparison');
        return;
    end
    
    mae_pop = mae_pop(valid);
    mae_van = mae_van(valid);
    mae_gre = mae_gre(valid);
    mae_k = mae_k(valid);
    mae_2d = mae_2d(valid);
    mae_ll = mae_ll(valid);
    N = length(mae_pop);
    
    figure('Name', 'Model Dimension vs Prediction Accuracy', 'Color', 'w', 'Position', [50 50 1600 600]);
    
    % === Panel 1: Boxplot comparison ===
    subplot(1, 3, 1);
    data = [mae_pop; mae_k; mae_2d; mae_ll; mae_van; mae_gre];
    groups = [repmat({'Pop (0D)'}, N, 1); ...
              repmat({'K (1D)'}, N, 1); ...
              repmat({'2D'}, N, 1); ...
              repmat({'LL (3D)'}, N, 1); ...
              repmat({'Van (4D)'}, N, 1); ...
              repmat({'Gre (4D)'}, N, 1)];
    boxplot(data, groups, 'Widths', 0.6);
    ylabel('MAE (BIS units)', 'FontSize', 12);
    title('Model Comparison by Dimension', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    means = [mean(mae_pop), mean(mae_k), mean(mae_2d), mean(mae_ll), mean(mae_van), mean(mae_gre)];
    hold on;
    plot(1:6, means, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    for ii = 1:6
        text(ii, means(ii) + 0.5, sprintf('%.2f', means(ii)), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    
    % === Panel 2: Per-patient comparison ===
    subplot(1, 3, 2);
    scatter(mae_van, mae_k, 80, 'c', 'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'K (1D)');
    hold on;
    scatter(mae_van, mae_2d, 80, [1 0.5 0], 'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', '2D');
    scatter(mae_van, mae_ll, 80, 'm', 'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'LL (3D)');
    scatter(mae_van, mae_gre, 80, 'r', 'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'Gre (4D)');
    plot([0 25], [0 25], 'k--', 'LineWidth', 2);
    xlabel('MAE Bouillon (4D)', 'FontSize', 12);
    ylabel('MAE of Other Models', 'FontSize', 12);
    title('Per-Patient: 4D vs Others', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on;
    axis equal;
    xlim([0 max([mae_van; mae_k; mae_2d; mae_ll; mae_gre]) * 1.1]);
    ylim([0 max([mae_van; mae_k; mae_2d; mae_ll; mae_gre]) * 1.1]);
    
    % === Panel 3: MAE vs Model Dimension ===
    subplot(1, 3, 3);
    dims = [0, 1, 2, 3, 4, 4];
    mean_maes = [mean(mae_pop), mean(mae_k), mean(mae_2d), mean(mae_ll), mean(mae_van), mean(mae_gre)];
    std_maes = [std(mae_pop), std(mae_k), std(mae_2d), std(mae_ll), std(mae_van), std(mae_gre)];
    
    % Plot Bouillon trajectory
    errorbar([0 1 2 3 4], mean_maes(1:5), std_maes(1:5), 'ko-', 'LineWidth', 2, 'MarkerSize', 12, 'MarkerFaceColor', 'b');
    hold on;
    % Add Greco as separate point
    errorbar(4.15, mean_maes(6), std_maes(6), 'rs', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    xlabel('Model Dimension (# parameters)', 'FontSize', 12);
    ylabel('Mean MAE ± Std', 'FontSize', 12);
    title('MAE vs Model Complexity', 'FontSize', 14, 'FontWeight', 'bold');
    xlim([-0.5 4.7]);
    grid on;
    
    yline(mean(mae_k), 'c--', 'LineWidth', 1.5);
    text(3.5, mean(mae_k) + 0.3, '1D baseline', 'Color', 'c', 'FontWeight', 'bold');
    
    legend({'Bouillon', 'Greco', '1D reference'}, 'Location', 'best');
    
    sgtitle(sprintf(['IDENTIFIABILITY ANALYSIS (N=%d patients)\n' ...
        'Adding parameters beyond 1D provides negligible improvement'], N), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    print(gcf, fullfile(fig_dir, 'identifiability_comparison.png'), '-dpng', '-r150');
    savefig(gcf, fullfile(fig_dir, 'identifiability_comparison.fig'));
end
