function plot_covariance_evolution(results, patient_id)
    idx = find(results.patient_id == patient_id, 1);
    if isempty(idx)
        warning('Patient not found');
        return;
    end
    
    r = results.raw(idx);
    t = r.time / 60;
    P_hist = r.Phist_van;
    
    if isempty(P_hist) || all(isnan(P_hist(:)))
        warning('No P matrix history found.');
        return;
    end
    
    param_names = {'C50 Prop', 'C50 Remi', 'Gamma', 'Synergy'};
    colors = {'b', 'r', 'g', 'k'};
    
    figure('Name', sprintf('Covariance Evolution - Case %d (V6.0)', patient_id), 'Color', 'w');
    
    for ii = 1:4
        subplot(2, 2, ii);
        data_to_plot = P_hist(ii, :);
        semilogy(t, data_to_plot(:), 'LineWidth', 2, 'Color', colors{ii});
        grid on;
        title(sprintf('%s Variance (P_{%d,%d})', param_names{ii}, ii, ii));
        xlabel('Time (min)');
        ylabel('Variance (\sigma^2)');
        
        final_var = P_hist(ii, end);
        if final_var < 0.05
            subtitle('Status: CONVERGED');
        elseif final_var < 0.5
            subtitle('Status: CONVERGING');
        else
            subtitle('Status: UNCERTAIN');
        end
    end
    
    sgtitle(sprintf('EKF Covariance Evolution (V6.0 Bounded) - Patient %d', patient_id));
end
