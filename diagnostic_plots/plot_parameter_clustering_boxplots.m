function param_stats = plot_parameter_clustering_boxplots(results, fig_dir)
% PLOT_PARAMETER_CLUSTERING_BOXPLOTS - Combined boxplot for 1D vs 2D
% Shows k, kP, kR in subplots to demonstrate non-identifiability

    N = length(results.patient_id);
    
    % Preallocate
    k_1d = nan(N, 1);
    kP_2d = nan(N, 1);
    kR_2d = nan(N, 1);
    a0_3d = nan(N, 1);
    aP_3d = nan(N, 1);
    aR_3d = nan(N, 1);
    
    % Extract final values for each patient
    for i = 1:N
        r = results.raw(i);
        
        if isfield(r, 'k_trajectory') && ~isempty(r.k_trajectory)
            idx = find(~isnan(r.k_trajectory), 1, 'last');
            if ~isempty(idx), k_1d(i) = r.k_trajectory(idx); end
        end
        
        if isfield(r, 'kP_fim_trajectory') && ~isempty(r.kP_fim_trajectory)
            idx = find(~isnan(r.kP_fim_trajectory), 1, 'last');
            if ~isempty(idx), kP_2d(i) = r.kP_fim_trajectory(idx); end
        end
        if isfield(r, 'kR_fim_trajectory') && ~isempty(r.kR_fim_trajectory)
            idx = find(~isnan(r.kR_fim_trajectory), 1, 'last');
            if ~isempty(idx), kR_2d(i) = r.kR_fim_trajectory(idx); end
        end
        
        if isfield(r, 'Xhist_loglin') && ~isempty(r.Xhist_loglin)
            idx = find(any(~isnan(r.Xhist_loglin), 1), 1, 'last');
            if ~isempty(idx)
                a0_3d(i) = r.Xhist_loglin(1, idx);
                aP_3d(i) = r.Xhist_loglin(2, idx);
                aR_3d(i) = r.Xhist_loglin(3, idx);
            end
        end
    end
    
    % Store results
    param_stats.k_1d = k_1d;
    param_stats.kP_2d = kP_2d;
    param_stats.kR_2d = kR_2d;
    param_stats.a0_3d = a0_3d;
    param_stats.aP_3d = aP_3d;
    param_stats.aR_3d = aR_3d;
    
    %% ==================== MAIN FIGURE: 1D vs 2D ====================
    figure('Name', '1D vs 2D Parameters', 'Color', 'w', 'Position', [50 50 600 350]);
    
    % Get common y-axis limits
    all_data = [k_1d; kP_2d; kR_2d];
    y_min = min(all_data) - 0.1;
    y_max = max(all_data) + 0.1;
    
    % Subplot 1: 1D model
    subplot(1, 2, 1);
    boxplot(k_1d, 'Widths', 0.5);
    hold on;
    yline(1.0, 'k--', 'LineWidth', 1.2);
    ylabel('Potency Scaling Factor', 'FontSize', 11);
    set(gca, 'XTickLabel', {'k'}, 'FontSize', 10);
    title('1D Model', 'FontSize', 12);
    ylim([y_min, y_max]);
    grid on;
    
    % Subplot 2: 2D model
    subplot(1, 2, 2);
    data_2d = [kP_2d, kR_2d];
    boxplot(data_2d, 'Widths', 0.5);
    hold on;
    yline(1.0, 'k--', 'LineWidth', 1.2);
    
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'k_P', 'k_R'}, ...
        'TickLabelInterpreter', 'tex', 'FontSize', 10);
    
    title('2D Model', 'FontSize', 12);
    ylim([y_min, y_max]);
    grid on;
    
    print(gcf, fullfile(fig_dir, 'figure_param_boxplots_1d2d.png'), '-dpng', '-r300');
    savefig(gcf, fullfile(fig_dir, 'figure_param_boxplots_1d2d.fig'));
    
    %% ==================== FIGURE 2: 3D MODEL (separate) ====================
    figure('Name', '3D Model Parameters', 'Color', 'w', 'Position', [150 150 450 350]);
    
    data_3d = [a0_3d, aP_3d, aR_3d];
    boxplot(data_3d, 'Widths', 0.5);
    
    set(gca, 'XTick', [1, 2, 3], 'XTickLabel', {'a_0', 'a_P', 'a_R'}, ...
        'TickLabelInterpreter', 'tex', 'FontSize', 10);
    
    ylabel('Parameter Value', 'FontSize', 11);
    title('3D Log-Linear Model', 'FontSize', 12);
    grid on;
    
    print(gcf, fullfile(fig_dir, 'figure_param_boxplots_3d.png'), '-dpng', '-r300');
    savefig(gcf, fullfile(fig_dir, 'figure_param_boxplots_3d.fig'));
    
    %% ==================== CONSOLE SUMMARY ====================
    fprintf('\n');
    fprintf('================================================================================\n');
    fprintf('                    PARAMETER CLUSTERING SUMMARY                               \n');
    fprintf('================================================================================\n\n');
    
    cv_1d = std(k_1d, 'omitnan') / mean(k_1d, 'omitnan') * 100;
    cv_2d_kP = std(kP_2d, 'omitnan') / mean(kP_2d, 'omitnan') * 100;
    cv_2d_kR = std(kR_2d, 'omitnan') / mean(kR_2d, 'omitnan') * 100;
    cv_3d_a0 = std(a0_3d, 'omitnan') / abs(mean(a0_3d, 'omitnan')) * 100;
    cv_3d_aP = std(aP_3d, 'omitnan') / mean(aP_3d, 'omitnan') * 100;
    cv_3d_aR = std(aR_3d, 'omitnan') / mean(aR_3d, 'omitnan') * 100;
    
    fprintf('%-12s %-10s %-10s %-10s %-10s\n', 'Parameter', 'Mean', 'SD', 'CV(%)', 'N');
    fprintf('%s\n', repmat('-', 1, 55));
    fprintf('%-12s %10.2f %10.2f %10.1f %10d\n', 'k (1D)', mean(k_1d,'omitnan'), std(k_1d,'omitnan'), cv_1d, sum(~isnan(k_1d)));
    fprintf('%-12s %10.2f %10.2f %10.1f %10d\n', 'kP (2D)', mean(kP_2d,'omitnan'), std(kP_2d,'omitnan'), cv_2d_kP, sum(~isnan(kP_2d)));
    fprintf('%-12s %10.2f %10.2f %10.1f %10d\n', 'kR (2D)', mean(kR_2d,'omitnan'), std(kR_2d,'omitnan'), cv_2d_kR, sum(~isnan(kR_2d)));
    fprintf('%-12s %10.2f %10.2f %10.1f %10d\n', 'a0 (3D)', mean(a0_3d,'omitnan'), std(a0_3d,'omitnan'), cv_3d_a0, sum(~isnan(a0_3d)));
    fprintf('%-12s %10.2f %10.2f %10.1f %10d\n', 'aP (3D)', mean(aP_3d,'omitnan'), std(aP_3d,'omitnan'), cv_3d_aP, sum(~isnan(aP_3d)));
    fprintf('%-12s %10.2f %10.2f %10.1f %10d\n', 'aR (3D)', mean(aR_3d,'omitnan'), std(aR_3d,'omitnan'), cv_3d_aR, sum(~isnan(aR_3d)));
    fprintf('\n');
    
    % Store CVs
    param_stats.cv.k_1d = cv_1d;
    param_stats.cv.kP_2d = cv_2d_kP;
    param_stats.cv.kR_2d = cv_2d_kR;
    param_stats.cv.a0_3d = cv_3d_a0;
    param_stats.cv.aP_3d = cv_3d_aP;
    param_stats.cv.aR_3d = cv_3d_aR;
end
