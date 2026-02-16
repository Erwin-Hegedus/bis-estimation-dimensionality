function generate_combined_identifiability_figure(results, cfg, fig_dir, target_pid)
% FIGURE 7: COMPREHENSIVE IDENTIFIABILITY ANALYSIS
% Merges FIM diagnostics and Error Landscape proof into one 2x2 figure.
% Target: Patient 105 (Representative Case)

    fprintf('Generating Figure 7: Combined Identifiability Analysis...\n');
    
    pid_idx = find([results.patient_id] == target_pid);
    if isempty(pid_idx)
        warning('Patient 105 not found, picking first valid long case.');
        [~, pid_idx] = max([results.demographics.Duration]);
        target_pid = results.patient_id(pid_idx);
    end
    r = results.raw(pid_idx);
    
    % --- 2. CREATE FIGURE ---
    fig = figure('Name', 'Figure 7: Combined Identifiability', 'Color', 'w', ...
        'Position', [50 50 1200 900]); % Nagyobb méret a 2x2-höz

    % =====================================================================
    % PANEL A: COHORT FIM EVOLUTION (WITH PATIENT HIGHLIGHT)
    % =====================================================================
    subplot(2, 2, 1);

    valid_idx = find(~isnan(results.metrics.van.MAE));
    t_norm = linspace(0, 1, 100);
    fim_matrix = nan(length(valid_idx), 100);

    % Collect Cohort Data
    for ii = 1:length(valid_idx)
        idx = valid_idx(ii);
        fim_hist = results.raw(idx).FIM_cond_hist;
        if ~isempty(fim_hist) && sum(~isnan(fim_hist)) > 10
            t_case = linspace(0, 1, length(fim_hist));
            fim_matrix(ii, :) = interp1(t_case, fim_hist, t_norm, 'linear', NaN);
        end
    end

    fim_median = median(fim_matrix, 1, 'omitnan');
    fim_q25 = prctile(fim_matrix, 25, 1);
    fim_q75 = prctile(fim_matrix, 75, 1);

    % Plot Cohort Ribbon
    fill([t_norm, fliplr(t_norm)], [fim_q25, fliplr(fim_q75)], ...
        [0.8 0.85 1.0], 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
    plot(t_norm, fim_median, 'b-', 'LineWidth', 2);

    if ~isempty(r.FIM_cond_hist)
        fim_pat = r.FIM_cond_hist;
        t_pat = linspace(0, 1, length(fim_pat));
        plot(t_pat, fim_pat, 'r-', 'LineWidth', 1.5);
    end

    yline(100, 'k--', 'LineWidth', 1.5);
    set(gca, 'YScale', 'log');
    xlim([0 1]); ylim([1e0 1e9]);
    xlabel('Normalized Time');
    ylabel('Condition Number \kappa (log)');
    title('(a) FIM Condition Number', 'FontSize', 12, 'FontWeight', 'bold');
    legend({'Cohort IQR', 'Cohort Median', 'Representative Patient', ...
        'Ill-cond. threshold'}, 'Location', 'SouthEast', 'FontSize', 8);
    grid on;

    % =====================================================================
    % PANEL B: PARAMETER VARIANCE EVOLUTION (PATIENT 105)
    % =====================================================================
    subplot(2, 2, 2);

    if ~isempty(r.Phist_van)
        t_min = r.time / 60;
        % Normalize variances to start at 1 for comparison
        P_norm = r.Phist_van ./ max(r.Phist_van(:,1), 1e-6); 
        
        param_names = {'C_{50P}', 'C_{50R}', '\gamma', '\beta'};
        colors = {'b', 'r', 'g', 'm'};
        
        for pp = 1:min(4, size(P_norm,1))
            semilogy(t_min, r.Phist_van(pp,:), 'Color', colors{pp}, 'LineWidth', 1.5, ...
                     'DisplayName', param_names{pp}); 
            hold on;
        end
        title('(b) 4D Parameter Variance', 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Time (min)');
        ylabel('Variance \sigma^2 (log)');
        legend('Location', 'NorthEast', 'FontSize', 8);
        grid on; xlim([0 max(t_min)]);
    else
        text(0.5,0.5,'No Variance Data');
    end

    % =====================================================================
    % PREPARE SNAPSHOT FOR PANELS C & D
    % =====================================================================
    % Pick a stable timepoint (e.g., 50% of case or max concentration)
    t_snap_idx = round(length(r.time) * 0.6); 
    CeP_ref = r.CeP_trajectory(t_snap_idx);
    CeR_ref = r.CeR_trajectory(t_snap_idx);
    E0_ref = r.E0_trajectory(t_snap_idx);
    BIS_meas = r.bis(t_snap_idx);
    
    % Ensure non-zero concentrations
    if CeP_ref < 0.5
        [~, t_snap_idx] = max(r.CeP_trajectory);
        CeP_ref = r.CeP_trajectory(t_snap_idx);
        CeR_ref = r.CeR_trajectory(t_snap_idx);
        BIS_meas = r.bis(t_snap_idx);
    end

    % =====================================================================
    % PANEL C: 1D CONVEXITY PROOF
    % =====================================================================
    subplot(2, 2, 3);
    
    k_range = linspace(0.2, 3.0, 100);
    sse_1d = zeros(size(k_range));
    for i = 1:length(k_range)
        bis_pred = predict_bis_proof_internal(k_range(i), k_range(i), CeP_ref, CeR_ref, E0_ref, cfg);
        sse_1d(i) = (bis_pred - BIS_meas)^2;
    end
    
    plot(k_range, sse_1d, 'b-', 'LineWidth', 2); hold on;
    [min_err, min_idx] = min(sse_1d);
    plot(k_range(min_idx), min_err, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    
    xlabel('Parameter k (1D Model)');
    ylabel('Squared Error (SSE)');
    title('(c) 1D Model: Strict Convexity', 'FontSize', 12, 'FontWeight', 'bold');
    grid on; axis tight;
    text(k_range(min_idx) + 0.25, 80, 'Global Minimum', ...
     'Color', 'r', 'FontSize', 9, 'FontWeight', 'bold');

    % =====================================================================
    % PANEL D: 2D FLAT VALLEY PROOF
    % =====================================================================
    subplot(2, 2, 4);
    
    kp_range = linspace(0.4, 2.8, 50);
    kr_range = linspace(0.4, 2.8, 50);
    [KP, KR] = meshgrid(kp_range, kr_range);
    sse_2d = zeros(size(KP));
    
    for i = 1:size(KP,1)
        for j = 1:size(KP,2)
            bis_pred = predict_bis_proof_internal(KP(i,j), KR(i,j), CeP_ref, CeR_ref, E0_ref, cfg);
            sse_2d(i,j) = (bis_pred - BIS_meas)^2;
        end
    end
    
    % Plot Contour
    contourf(KP, KR, log10(sse_2d + 1e-1), 20, 'LineStyle', 'none'); 
    colormap(gca, 'parula');
    c = colorbar; c.Label.String = 'log(SSE)';
    hold on;
    
    % Highlight the valley floor
    min_val_2d = min(sse_2d(:));
    [r_min, c_min] = find(sse_2d < min_val_2d + 1.0); % Tolerance zone
    plot(kp_range(c_min), kr_range(r_min), 'w.', 'MarkerSize', 2, 'DisplayName', 'Valley Floor');
    
    xlabel('Parameter k_P');
    ylabel('Parameter k_R');
    title('(d) 2D Error Landscape', 'FontSize', 12, 'FontWeight', 'bold');
    axis square;
    
    % Add correlation line annotation
    text(1.5, 1.5, 'Flat Valley', 'Color', 'k', 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Rotation', 45, 'FontSize', 10);

    % --- SAVE ---
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
    saveas(fig, fullfile(fig_dir, 'figure7_combined_identifiability.png'));
    saveas(fig, fullfile(fig_dir, 'figure7_combined_identifiability.fig'));
    fprintf('Saved: figure7_combined_identifiability.png\n');
end
