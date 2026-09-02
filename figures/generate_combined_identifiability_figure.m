function generate_combined_identifiability_figure(results, cfg, fig_dir, target_pid)
% FIGURE 7: COMPREHENSIVE IDENTIFIABILITY ANALYSIS
% Merges FIM diagnostics and Error Landscape proof into one 2x2 figure.
% Target: the patient passed in target_pid (Representative Case)

    fprintf('Generating Figure 7: Combined Identifiability Analysis...\n');
    
    pid_idx = find([results.patient_id] == target_pid);
    if isempty(pid_idx)
        warning('Patient %d not in results, picking longest valid case.', target_pid);
        [~, pid_idx] = max([results.demographics.Duration]);
        target_pid = results.patient_id(pid_idx);
    end
    r = results.raw(pid_idx);
    
    % --- 2. CREATE FIGURE ---
    fig = figure('Name', 'Figure 7: Combined Identifiability', 'Color', 'w', ...
        'Position', [50 50 1200 900]);

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

    yline(100, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    set(gca, 'YScale', 'log');
    % Upper limit follows the data: both the cohort ribbon and the highlighted
    % case run past 1e11, and a fixed 1e9 ceiling clipped the red trace.
    y_top = 10^ceil(log10(max([fim_q75(:); r.FIM_cond_hist(:)], [], 'omitnan')));
    xlim([0 1]); ylim([1e0 max(1e9, y_top)]);
    xlabel('Normalized time');
    ylabel('Condition number \kappa (log)');
    title('(a) FIM condition number');
    legend({'IQR', 'Median', sprintf('Patient %d', target_pid)}, ...
        'Location', 'SouthEast', 'Box', 'on', 'Color', 'w', 'NumColumns', 1);
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
        title('(b) 4D parameter variance');
        xlabel('Time (min)');
        ylabel('Variance \sigma^2 (log)');
        % A decade of headroom so the legend row clears the traces. The beta
        % variance is nearly flat across the top of the panel, so every in-axes
        % legend position landed on a curve.
        yl = ylim; ylim([yl(1) yl(2)*20]);
        legend('Location', 'North', 'Box', 'off', 'NumColumns', 2);
        grid on; xlim([0 max(t_min)]);
    else
        text(0.5,0.5,'No Variance Data');
    end

    % =====================================================================
    % TRAJECTORY FOR PANELS C & D
    % =====================================================================
    % Landscapes are evaluated over the whole case. A single sample gives an
    % equal-output curve for any two parameters and says nothing about
    % identifiability from the trajectory.
    ok = ~isnan(r.bis) & ~isnan(r.CeP_trajectory) & ~isnan(r.CeR_trajectory) & ...
         ~isnan(r.E0_trajectory) & r.CeP_trajectory > 0.5;
    idx = find(ok);
    stride = max(1, floor(numel(idx)/600));
    idx = idx(1:stride:end);

    bis_t = r.bis(idx);
    CeP_t = r.CeP_trajectory(idx);
    CeR_t = r.CeR_trajectory(idx);
    E0_t  = r.E0_trajectory(idx);
    n_t   = numel(idx);

    % =====================================================================
    % PANEL C: 1D CONVEXITY PROOF
    % =====================================================================
    subplot(2, 2, 3);
    
    k_range = linspace(0.2, 6.0, 150);
    mae_1d = zeros(size(k_range));
    for i = 1:length(k_range)
        p = predict_bis_2d_model(k_range(i), k_range(i), CeP_t, CeR_t, E0_t, cfg.BISmin_fixed, cfg);
        mae_1d(i) = mean(abs(p - bis_t));
    end

    plot(k_range, mae_1d, 'b-', 'LineWidth', 2); hold on;
    [min_err, min_idx] = min(mae_1d);
    plot(k_range(min_idx), min_err, 'rx', 'MarkerSize', 9, 'LineWidth', 1.6);

    xlabel('Parameter k (1D model)');
    ylabel('MAE (BIS)');
    title('(c) 1D error landscape');
    grid on; axis tight;

    % =====================================================================
    % PANEL D: 2D (k, gamma) LANDSCAPE
    % =====================================================================
    subplot(2, 2, 4);

    % k spans the same range as panel (c) so the two landscapes can be read
    % against each other; this case's minimum lies above the estimator's
    % k bound and would sit off the edge of a bounds-width sweep.
    k_range_2d  = linspace(0.2, 6.0, 50);
    ga_range_2d = linspace(1.0, 6.0, 50);
    [KK, GG] = meshgrid(k_range_2d, ga_range_2d);
    mae_2d = zeros(size(KK));

    % predict_bis_kgamma_model is scalar-valued, so the trajectory is walked
    % sample by sample rather than passed in as a vector.
    p = zeros(n_t, 1);
    for i = 1:size(KK,1)
        for j = 1:size(KK,2)
            for s = 1:n_t
                p(s) = predict_bis_kgamma_model(KK(i,j), GG(i,j), ...
                    CeP_t(s), CeR_t(s), E0_t(s), cfg.BISmin_fixed, cfg);
            end
            mae_2d(i,j) = mean(abs(p - bis_t(:)));
        end
    end

    contourf(KK, GG, mae_2d, 20, 'LineStyle', 'none', 'HandleVisibility', 'off');
    colormap(gca, 'parula');
    c = colorbar; c.Label.String = 'MAE (BIS)';
    hold on;

    % Sets within 1 BIS of the best fit: the practically indistinguishable region
    min_val_2d = min(mae_2d(:));
    [r_min, c_min] = find(mae_2d < min_val_2d + 1.0);
    plot(k_range_2d(c_min), ga_range_2d(r_min), 'w.', 'MarkerSize', 2, ...
         'DisplayName', 'within 1 BIS');

    xlabel('Parameter k');
    ylabel('Hill coefficient \gamma');
    title('(d) 2D error landscape');
    % No axis square: it shrank the axes away from the colorbar and forced the
    % tick labels to tilt. The tile is already close to square.
    legend('Location', 'NorthEast', 'Box', 'off');

    % --- SAVE ---
    % No bare apply_figure_style(fig) here: it sets the 'styled' guard, which
    % made the sized call below a silent no-op and shipped this figure at the
    % fallback geometry with (a) and (b) overprinting each other.
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
    export_figure_png(fig, fullfile(fig_dir, 'figure7_combined_identifiability.png'), 'double', 4.0);
    savefig(fig, fullfile(fig_dir, 'figure7_combined_identifiability.fig'));
    fprintf('Saved: figure7_combined_identifiability.png\n');
end
