function generate_figure4_case_study(results, patient_id, cfg, fig_dir)
% FIGURE 4: Representative Case Study of Online Parameter Adaptation
% Multi-panel figure showing BIS, concentrations, and parameter evolution
%
% Two-column float (figure*) laid out with tiledlayout. In a single column each
% panel was 1.6 x 1.5 in, so a 10 pt label was a seventh of the panel height and
% the titles had to be abbreviated to fit; at \textwidth the same 10 pt text sits
% on a 3.3 in panel and the figure reads as an ordinary plot. tiledlayout rather
% than subplot because subplot's default margins spend a third of the width on
% whitespace.

    fprintf('Generating Figure 4: Case Study (Patient %d)...\n', patient_id);

    idx = find(results.patient_id == patient_id, 1);
    if isempty(idx)
        warning('Patient %d not found, using first valid patient', patient_id);
        idx = find(~isnan(results.metrics.van.MAE), 1);
        patient_id = results.patient_id(idx);
    end

    r = results.raw(idx);

    % CUT LAST 15% (emergence phase)
    n_total = length(r.time);
    n_keep = floor(0.83 * n_total);
    t = r.time(1:n_keep) / 60;
    tmax = max(t);

    fig = figure('Name', sprintf('Figure 4: Case Study (Patient %d)', patient_id), ...
                 'Color', 'w');
    tl = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Panel (a): BIS predictions, spanning the full width
    ax = nexttile(tl, 1, [1 2]);
    hM = plot(t, r.bis(1:n_keep), 'k-', 'LineWidth', 1.0); hold on;
    h0 = plot(t, r.pred_pop(1:n_keep), 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
    h1 = plot(t, r.pred_kscale(1:n_keep), 'c-');
    h2 = plot(t, r.pred_kgamma(1:n_keep), 'b-');

    % Target band: drawn, not in the legend. Two dashed grey lines at 40 and 60
    % need one clause of caption, not a fifth of the legend row.
    yline(cfg.target_band(1), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
    yline(cfg.target_band(2), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);

    xlim([0 tmax]);
    ylabel('BIS');
    ylim([0 100]);
    title('(a) BIS prediction comparison');
    grid on;
    set(ax, 'XTickLabel', []);

    % One legend for the figure, on its own strip above the layout, so it never
    % has to share a panel with the data.
    lg = legend([hM h0 h1 h2], ...
                {'Measured BIS', '0D (population)', '1D (k)', '2D (k, \gamma)'}, ...
                'Orientation', 'horizontal', 'Box', 'off');
    lg.Layout.Tile = 'north';

    % Panel (b): Effect-site concentrations. The axis colours identify the two
    % drugs, so the legend that used to sit here is redundant.
    nexttile(tl);
    yyaxis left;
    plot(t, r.CeP_trajectory(1:n_keep), 'b-');
    ylabel('C_{eP} (\mug/mL)');
    ylim([0 8]);
    set(gca, 'YColor', 'b');

    yyaxis right;
    plot(t, r.CeR_trajectory(1:n_keep) * 1000, 'r-');
    ylabel('C_{eR} (ng/mL)');
    ylim([0 15]);
    set(gca, 'YColor', [0.8 0 0]);

    xlim([0 tmax]);
    title('(b) Effect-site concentrations');
    grid on;

    % Panel (c): 1D k trajectory
    nexttile(tl);
    plot(t, r.k_trajectory(1:n_keep), 'c-');
    hold on;
    yline(1.0, 'k--', 'LineWidth', 1.0);
    xlim([0 tmax]);
    ylabel('k_{scale}');
    title('(c) 1D potency factor evolution');
    ylim([0.3 2.5]);
    grid on;

    % Panel (d): the two 2D coordinates against the 1D estimate. Potency and
    % steepness are different quantities, so gamma carries its own axis.
    nexttile(tl);
    yyaxis left;
    plot(t, r.kg_k_trajectory(1:n_keep), 'b-', 'DisplayName', 'k (2D)');
    hold on;
    plot(t, r.k_trajectory(1:n_keep), 'c--', 'DisplayName', 'k (1D)');
    yline(1.0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
    ylabel('Potency scale');
    ylim([0.3 2.0]);

    yyaxis right;
    plot(t, r.kg_gamma_trajectory(1:n_keep), 'r-', 'DisplayName', '\gamma (2D)');
    yline(cfg.population_params_van(3), 'r:', 'LineWidth', 1, 'HandleVisibility', 'off');
    ylabel('Hill coefficient \gamma');
    ylim([0.8 3.2]);

    xlim([0 tmax]);
    xlabel('Time (min)');
    title('(d) 2D vs 1D parameter comparison');
    legend('Location', 'best', 'Orientation', 'horizontal', 'Box', 'on');
    grid on;

    % Panel (e): Variance evolution
    nexttile(tl);
    if ~isempty(r.k_var)
        n_keep_var = min(n_keep, length(r.k_var));
        t_var = (1:n_keep_var)' / 60;
        semilogy(t_var, r.k_var(1:n_keep_var), 'c-', 'DisplayName', 'Var(k), 1D');
        hold on;
    end
    if ~isempty(r.kg_k_var)
        n_keep_2d = min(n_keep, length(r.kg_k_var));
        t_2d = (1:n_keep_2d)' / 60;
        semilogy(t_2d, r.kg_k_var(1:n_keep_2d), 'b-', 'DisplayName', 'Var(log k), 2D');
        semilogy(t_2d, r.kg_gamma_var(1:n_keep_2d), 'r-', 'DisplayName', 'Var(log \gamma), 2D');
    end

    xlim([0 tmax]);
    yl = ylim; ylim([yl(1) yl(2)*2]);
    ylabel('Parameter variance (log)');
    xlabel('Time (min)');
    title('(e) Uncertainty evolution');
    grid on;

    save_figure(fig, fig_dir, 'figure4_case_study', 'double', 4.5);
end
