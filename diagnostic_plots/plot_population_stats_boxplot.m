function plot_population_stats_boxplot(results)
    N = length(results.patient_id);
    final_params_van = [];
    final_params_gre = [];
    
    for ii = 1:N
        if ~isempty(results.raw(ii).Xhist_van) && size(results.raw(ii).Xhist_van, 2) > 100
            last_vals = mean(results.raw(ii).Xhist_van(:, end-60:end), 2, 'omitnan');
            final_params_van = [final_params_van; last_vals'];
            
            last_vals_g = mean(results.raw(ii).Xhist_gre(:, end-60:end), 2, 'omitnan');
            final_params_gre = [final_params_gre; last_vals_g'];
        end
    end
    
    param_names = {'C50 Prop (\mug/ml)', 'C50 Remi (ng/ml)', 'Gamma', 'Synergy'};
    
    figure('Name', 'Population PD Parameter Distribution (V6.0)', 'Color', 'w', 'Position', [100 100 1000 600]);
    
    for p = 1:4
        subplot(2, 2, p);
        
        data_col = [final_params_van(:, p); final_params_gre(:, p)];
        groups = [repmat({'Bouillon'}, size(final_params_van, 1), 1); ...
                  repmat({'Greco'}, size(final_params_gre, 1), 1)];
        
        boxplot(data_col, groups, 'Widths', 0.5);
        title(param_names{p}, 'FontWeight', 'bold');
        grid on;
        
        med_v = median(final_params_van(:, p));
        med_g = median(final_params_gre(:, p));
        
    end
   
end
