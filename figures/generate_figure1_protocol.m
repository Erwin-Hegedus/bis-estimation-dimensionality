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
    
    % Two-column float (figure*), so each panel gets about 3.0 x 1.6 in and the
    % descriptive titles and the histogram's y label all have room. tiledlayout
    % rather than subplot: subplot's default margins spend a third of the width
    % on whitespace, which is the width the 10 pt labels need.
    fig = figure('Name', 'Figure 1: Protocol', 'Color', 'w');
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    t_min = rep_time / 60;

    % (a) BIS, top left
    nexttile(tl, 1);
    plot(t_min, rep_bis, 'k-', 'LineWidth', 1.5);
    hold on;
    yline(40, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    yline(60, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    ylabel('BIS');
    ylim([0 100]);
    title(sprintf('(a) Representative case %d: BIS response', representative_patient));
    grid on;
    box on;
    xlim([0 max(t_min)]);
    set(gca, 'XTickLabel', []);  % shared with the concentration panel below

    % (b) Concentrations, bottom left, same time axis
    nexttile(tl, 3);

    yyaxis left;
    plot(t_min, rep_p_ce, 'b-');
    ylabel('Propofol C_{eP} (\mug/mL)');
    ylim([0 max(rep_p_ce) * 1.2]);
    set(gca, 'YColor', 'b');

    yyaxis right;
    plot(t_min, rep_r_ce, 'r-');
    ylabel('Remifentanil C_{eR} (ng/mL)');
    ylim([0 max(rep_r_ce) * 1.2]);
    set(gca, 'YColor', 'r');

    xlabel('Time (min)');
    % No panel letter: the caption folds this and the BIS trace above it into (a).
    title('Effect-site concentrations');
    grid on;
    box on;
    xlim([0 max(t_min)]);

    % (c) Cohort correlation, right column, spanning both rows
    nexttile(tl, 2, [2 1]);
    histogram(rho_patients, 20, 'FaceColor', [0.3 0.5 0.9], 'EdgeColor', 'k', ...
              'FaceAlpha', 0.7, 'Normalization', 'probability');
    hold on;
    
    med = median(rho_patients, 'omitnan');
    xline(med, 'r-');

    % Left of the median line: the distribution piles up against rho = 1, so a
    % centred label sits on top of the tallest bins.
    text(med - 0.04, max(ylim)*0.95, sprintf('Median %.2f', med), ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
         'Color', 'r', 'BackgroundColor', 'w', 'EdgeColor', 'r');

    xlabel('Correlation coefficient \rho');
    ylabel('Probability');
    title(sprintf('(b) C_{eP}-C_{eR} correlation across cohort (N=%d)', 209));
    grid on;
    box on;
    xlim([-1 1]);

    save_figure(fig, fig_dir, 'figure1_protocol', 'double', 3.3);
    fprintf('  Saved: figure1_protocol.png\n');
end
