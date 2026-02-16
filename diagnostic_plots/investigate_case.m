function investigate_case(results, n)
    pid = n;
    idx = find(results.patient_id == pid, 1);
    
    if isempty(idx)
        fprintf('Patient %d not found in results.\n', pid);
        return;
    end
    
    r = results.raw(idx);
    t_min = r.time / 60;
    
    figure('Name', sprintf('Forensic Analysis: Case %d', pid), 'Color', 'w', 'Position', [50 50 1200 1100]);
    
    subplot(6, 1, 1);
    plot(t_min, r.bis, 'k-', 'LineWidth', 0.8, 'Color', [0.3 0.3 0.3]);
    hold on;
    if isfield(r, 'pred_pop') && ~isempty(r.pred_pop)
        plot(t_min, r.pred_pop, 'g:', 'LineWidth', 1.5);
    end
    if isfield(r, 'pred_kscale') && ~isempty(r.pred_kscale)
        plot(t_min, r.pred_kscale, 'c-', 'LineWidth', 2.0);
    end
    if isfield(r, 'pred_2d') && ~isempty(r.pred_2d)
        plot(t_min, r.pred_2d, 'Color', [1 0.5 0], 'LineWidth', 2.0);  % Orange
    end
    plot(t_min, r.pred_van, 'b-', 'LineWidth', 1.5);
    plot(t_min, r.pred_gre, 'r--', 'LineWidth', 1.5);
    if isfield(r, 'BISmin_trajectory') && ~isempty(r.BISmin_trajectory)
        plot(t_min, r.BISmin_trajectory, 'm--', 'LineWidth', 1.0);
        legend('Meas', 'Pop', 'K(1D)', '2D', 'Van(4D)', 'Gre', 'BISmin', 'Location', 'best');
    else
        legend('Meas', 'Pop', 'K(1D)', '2D', 'Van(4D)', 'Gre', 'Location', 'best');
    end
    title(sprintf('Case %d: 1D vs 2D vs 4D Model Comparison', pid), 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('BIS');
    ylim([0 100]);
    grid on;
    
    subplot(6, 1, 2);
    yyaxis left;
    plot(t_min, r.CeP_trajectory, 'b', 'LineWidth', 1.2);
    ylabel('Prop Ce (\mug/ml)');
    ylim([0 8]);
    set(gca, 'YColor', 'b');
    yyaxis right;
    plot(t_min, r.CeR_trajectory * 1000, 'r', 'LineWidth', 1.2);
    ylabel('Remi Ce (ng/ml)');
    ylim([0 15]);
    set(gca, 'YColor', [0.8 0 0]);
    title('Effect Site Concentrations');
    grid on;
    
    subplot(6, 1, 3);
    if isfield(r, 'ke0_trajectory') && ~isempty(r.ke0_trajectory)
        plot(t_min, r.ke0_trajectory, 'm', 'LineWidth', 1.5);
        yline(0.456, 'k--', 'Schnider Pop');
        title('ke0 Adaptation');
        ylabel('ke0 (/min)');
        grid on;
        ylim([0.05 0.8]);
    end
    
    subplot(6, 1, 4);
    hold on;
    if isfield(r, 'k_trajectory') && ~isempty(r.k_trajectory)
        plot(t_min, r.k_trajectory, 'c-', 'LineWidth', 2, 'DisplayName', 'k (1D)');
    end
    if isfield(r, 'kP_trajectory') && ~isempty(r.kP_trajectory)
        plot(t_min, r.kP_trajectory, 'b-', 'LineWidth', 2, 'DisplayName', 'kP (2D)');
        plot(t_min, r.kR_trajectory, 'r-', 'LineWidth', 2, 'DisplayName', 'kR (2D)');
    end
    yline(1.0, 'k--');
    title('Potency Scale Parameters: 1D vs 2D');
    ylabel('k');
    ylim([0.3 3]);
    legend('Location', 'best');
    grid on;
    
    subplot(6, 1, 5);
    hold on;
    if isfield(r, 'BISmin_trajectory') && ~isempty(r.BISmin_trajectory)
        plot(t_min, r.BISmin_trajectory, 'm', 'LineWidth', 2);
    end
    if isfield(r, 'E0_trajectory') && ~isempty(r.E0_trajectory)
        plot(t_min, r.E0_trajectory, 'c', 'LineWidth', 2);
    end
    yline(32, 'k--', 'BISmin Pop');
    yline(93, 'k:', 'E0 Pop');
    title('E0 and BISmin Estimation');
    ylabel('BIS');
    ylim([0 100]);
    legend('BISmin', 'E0', 'Location', 'east');
    grid on;
    
    param_names = {'C50P', 'C50R', 'Gamma', 'Alpha'};
    pop_vals_van = [3.5, 5.0, 2.8, 0.8];
    XV = r.Xhist_van;
    XG = r.Xhist_gre;
    
    for p = 1:4
        subplot(6, 4, 20 + p);
        hold on;
        grid on;
        if ~isempty(XV)
            plot(t_min, XV(p,:), 'b-', 'LineWidth', 1.2);
        end
        if ~isempty(XG)
            plot(t_min, XG(p,:), 'r--', 'LineWidth', 1.2);
        end
        yline(pop_vals_van(p), 'k:');
        title(param_names{p}, 'FontSize', 10);
        if p == 1
            ylabel('Value');
        end
        if p == 4
            xlabel('Time (min)');
        end
    end
end
