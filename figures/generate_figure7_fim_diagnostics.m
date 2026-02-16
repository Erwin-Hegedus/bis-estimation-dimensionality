function generate_figure7_fim_diagnostics(results, fig_dir)
% FIGURE 7: Diagnostics of Identifiability
% (a) FIM condition number evolution (cohort)
% (b) Parameter variance evolution

    fprintf('Generating Figure 7: FIM Diagnostics...\n');
    
    valid_idx = find(~isnan(results.metrics.van.MAE));
    
    fig = figure('Name', 'Figure 7: FIM Diagnostics', 'Color', 'w', ...
                 'Position', [50 50 1200 500]);
    
    % Panel (a): FIM condition number - cohort ribbon
    subplot(1, 2, 1);
    
    % Collect FIM condition numbers across cases (normalize time to 0-1)
    max_len = 0;
    for ii = 1:length(valid_idx)
        idx = valid_idx(ii);
        if ~isempty(results.raw(idx).FIM_cond_hist)
            max_len = max(max_len, length(results.raw(idx).FIM_cond_hist));
        end
    end
    
    if max_len > 0
        % Interpolate all cases to common time grid
        t_norm = linspace(0, 1, 100);
        fim_matrix = nan(length(valid_idx), 100);
        
        for ii = 1:length(valid_idx)
            idx = valid_idx(ii);
            fim_hist = results.raw(idx).FIM_cond_hist;
            if ~isempty(fim_hist) && sum(~isnan(fim_hist)) > 10
                t_case = linspace(0, 1, length(fim_hist));
                fim_matrix(ii, :) = interp1(t_case, fim_hist, t_norm, 'linear', NaN);
            end
        end
        
        % Compute statistics
        fim_median = median(fim_matrix, 1, 'omitnan');
        fim_q25 = prctile(fim_matrix, 25, 1);
        fim_q75 = prctile(fim_matrix, 75, 1);
        
        % Plot ribbon
        fill([t_norm, fliplr(t_norm)], [fim_q25, fliplr(fim_q75)], ...
             [0.3 0.5 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        hold on;
        semilogy(t_norm, fim_median, 'b-', 'LineWidth', 2);
        yline(100, 'r--', 'LineWidth', 1.5);
        
        xlabel('Normalized Time (0 = start, 1 = end)', 'FontSize', 11);
        ylabel('FIM Condition Number (log scale)', 'FontSize', 11);
        title('(a) FIM Condition Number Evolution', 'FontSize', 13, 'FontWeight', 'bold');
        legend({'IQR', 'Median', 'Ill-conditioned threshold'}, 'Location', 'best');
        set(gca, 'YScale', 'log');
        grid on;
    else
        text(0.5, 0.5, 'FIM data not available', 'HorizontalAlignment', 'center');
    end
    
    % Panel (b): Parameter variance evolution for representative case
    subplot(1, 2, 2);
    
    % Find a case with good P_hist
    rep_idx = [];
    for ii = 1:length(valid_idx)
        idx = valid_idx(ii);
        if ~isempty(results.raw(idx).Phist_van) && ...
           sum(~isnan(results.raw(idx).Phist_van(:))) > 100
            rep_idx = idx;
            break;
        end
    end
    
    if ~isempty(rep_idx)
        r = results.raw(rep_idx);
        t = r.time / 60;
        P_hist = r.Phist_van;
        
        param_names = {'C_{50P}', 'C_{50R}', '\gamma', '\beta'};
        colors = {'b', 'r', 'g', 'm'};
        
        for pp = 1:4
            semilogy(t, P_hist(pp,:), 'Color', colors{pp}, 'LineWidth', 1.5, ...
                     'DisplayName', param_names{pp});
            hold on;
        end
        
        xlabel('Time (min)', 'FontSize', 11);
        ylabel('Parameter Variance (log scale)', 'FontSize', 11);
        title(sprintf('(b) 4D Parameter Variance (Patient %d)', results.patient_id(rep_idx)), ...
              'FontSize', 13, 'FontWeight', 'bold');
        legend('Location', 'best');
        grid on;
    else
        text(0.5, 0.5, 'Variance history not available', 'HorizontalAlignment', 'center');
    end
    
    sgtitle('Figure 7: Identifiability Diagnostics', 'FontSize', 15, 'FontWeight', 'bold');
    
    save_figure(fig, fig_dir, 'figure7_fim_diagnostics');
end
