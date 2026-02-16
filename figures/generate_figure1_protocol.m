function generate_figure1_protocol(patientDataFinal, fig_dir)
    fprintf('Generating Figure 1: Protocol characteristics...\n');
    
    N = height(patientDataFinal);
    rho_patients = nan(N,1);
    durations = nan(N,1);
    
    for i = 1:N
        [p_ce, ~, ~, time] = getDataForPatient(patientDataFinal, i);
        if ~isempty(time)
            durations(i) = time(end) / 60;
        end
    end
    
    valid_cases = find(durations >= 15);
    fprintf('  Total cases: %d, Valid (>=15 min): %d\n', N, numel(valid_cases));
    
    representative_patient = [];
    rep_p_ce = [];
    rep_r_ce = [];
    rep_bis = [];
    rep_time = [];
    
    for i = 1:min(50, numel(valid_cases))
        idx = valid_cases(i);
        [p_ce, r_ce, bis, time] = getDataForPatient(patientDataFinal, idx);
        if ~isempty(p_ce) && numel(p_ce) > 900
            representative_patient = idx;
            rep_p_ce = p_ce;
            rep_r_ce = r_ce;
            rep_bis = bis;
            rep_time = time;
            fprintf('  Representative: Patient %d (%.1f min)\n', idx, time(end)/60);
            break;
        end
    end
    
    for ii = 1:numel(valid_cases)
        i = valid_cases(ii);
        [p_ce, r_ce, ~, ~] = getDataForPatient(patientDataFinal, i);
        
        if isempty(p_ce) || numel(p_ce) < 100
            continue;
        end
        
        valid = ~isnan(p_ce) & ~isnan(r_ce) & p_ce > 0.1 & r_ce > 0.01;
        p_filtered = p_ce(valid);
        r_filtered = r_ce(valid);
        
        if numel(p_filtered) < 50
            continue;
        end
        
        R = corrcoef(p_filtered, r_filtered, 'Rows', 'complete');
        rho_patients(i) = R(1,2);
    end
    
    rho_patients = rho_patients(~isnan(rho_patients));
    fprintf('  Correlations: N=%d, Median=%.3f\n', numel(rho_patients), median(rho_patients));
    
    if numel(rho_patients) < 5
        warning('Too few valid correlations.');
        return;
    end
    
    fig = figure('Name', 'Figure 1: Protocol', 'Color', 'w', 'Position', [50 50 1400 600]);
    
    t_min = rep_time / 60;
    
    % TOP LEFT: BIS
    subplot(2, 2, 1);
    plot(t_min, rep_bis, 'k-', 'LineWidth', 1.5);
    hold on;
    yline(40, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    yline(60, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    ylabel('BIS', 'FontSize', 11);
    ylim([0 100]);
    title(sprintf('(a) Representative Case (Patient %d): BIS Response', representative_patient), ...
          'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    box on;
    set(gca, 'XTickLabel', []);  % Remove x-axis labels (shared with plot below)
    
    % BOTTOM LEFT: Concentrations (dual y-axis)
    subplot(2, 2, 3);
    
    yyaxis left;
    plot(t_min, rep_p_ce, 'b-', 'LineWidth', 2.5);
    ylabel('Propofol C_{eP} (\mug/mL)', 'FontSize', 11);
    ylim([0 max(rep_p_ce) * 1.2]);
    set(gca, 'YColor', 'b');
    
    yyaxis right;
    plot(t_min, rep_r_ce, 'r-', 'LineWidth', 2.5);
    ylabel('Remifentanil C_{eR} (ng/mL)', 'FontSize', 11);
    ylim([0 max(rep_r_ce) * 1.2]);
    set(gca, 'YColor', 'r');
    
    xlabel('Time (min)', 'FontSize', 11);
    title('Effect-Site Concentrations', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    box on;
    
    % RIGHT: Correlation (spans both rows)
    subplot(2, 2, [2 4]);
    histogram(rho_patients, 20, 'FaceColor', [0.3 0.5 0.9], 'EdgeColor', 'k', ...
              'FaceAlpha', 0.7, 'Normalization', 'probability');
    hold on;
    
    med = median(rho_patients, 'omitnan');
    xline(med, 'r-', 'LineWidth', 2.5);
    
    text(med, max(ylim)*0.92, sprintf('Median = %.2f', med), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', 'r', ...
         'BackgroundColor', 'w', 'EdgeColor', 'r', 'LineWidth', 1.5);
    
    xlabel('Correlation Coefficient (\rho)', 'FontSize', 11);
    ylabel('Probability', 'FontSize', 11);
    title(sprintf('(b) Correlation Distribution (N=%d)', 209), ...
          'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    box on;
    xlim([-1 1]);
    
    save_figure(fig, fig_dir, 'figure1_protocol');
    fprintf('  Saved: figure1_protocol.png\n');
end
