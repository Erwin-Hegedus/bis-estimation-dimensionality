function plot_fim_diagnostics(results, patient_id)
    idx = find(results.patient_id == patient_id, 1);
    if isempty(idx)
        warning('Patient %d not found', patient_id);
        return;
    end
    
    r = results.raw(idx);
    t = r.time / 60;
    
    figure('Name', sprintf('FIM Diagnostics - Case %d', patient_id), ...
           'Color', 'w', 'Position', [100 100 1200 800]);
    
    subplot(3, 2, 1);
    if isfield(r, 'FIM_cond_hist') && ~isempty(r.FIM_cond_hist)
        semilogy(t, r.FIM_cond_hist, 'b-', 'LineWidth', 1.5);
        hold on;
        yline(100, 'r--', 'Ill-conditioned threshold');
        title('FIM Condition Number (Lower = Better Identifiability)');
        ylabel('cond(FIM)');
        xlabel('Time (min)');
        grid on;
    else
        text(0.5, 0.5, 'FIM history not available', 'HorizontalAlignment', 'center');
    end
    
    subplot(3, 2, 2);
    P_hist = r.Phist_van;
    if ~isempty(P_hist)
        param_names = {'C50P', 'C50R', '\gamma', '\beta'};
        colors = {'b', 'r', 'g', 'm'};
        for ii = 1:4
            semilogy(t, P_hist(ii,:), 'Color', colors{ii}, 'LineWidth', 1.5);
            hold on;
        end
        legend(param_names, 'Location', 'best');
        title('Parameter Variance Evolution (P_{ii})');
        ylabel('Variance');
        xlabel('Time (min)');
        grid on;
    end
    
    param_names = {'C50P (\mug/ml)', 'C50R (ng/ml)', '\gamma', '\beta'};
    XV = r.Xhist_van;
    
    for p = 1:4
        subplot(3, 2, 2 + p);
        if ~isempty(XV)
            plot(t, XV(p,:), 'b-', 'LineWidth', 1.5);
            hold on;
            
            if ~isempty(P_hist)
                std_p = sqrt(P_hist(p,:));
                fill([t(:); flipud(t(:))], ...
                     [XV(p,:)' + 2*std_p'; flipud(XV(p,:)' - 2*std_p')], ...
                     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            end
            
            title(sprintf('%s with 95%% CI', param_names{p}));
            ylabel('Value');
            xlabel('Time (min)');
            grid on;
        end
    end
    
    sgtitle(sprintf('V6.0 Identifiability Diagnostics - Patient %d', patient_id));
end
