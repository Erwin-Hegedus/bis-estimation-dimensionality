function plot_case_studies(cases, results, cfg, fig_dir)
%PLOT_CASE_STUDIES  Concentrations and 4D BIS predictions for m selected cases.
%
%   Two-column float (figure*) on tiledlayout. On subplot the five-entry legend
%   was wider than the panel it sat in and spilled across the neighbouring one;
%   here it gets its own strip above the layout, and the time label is shared.

    m = numel(cases.idx);
    if m == 0, return; end

    f = figure('Color', 'w');
    tl = tiledlayout(f, m, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    h = gobjects(0);
    for k = 1:m
        i = cases.idx(k);

        % CUT LAST 15% (emergence phase)
        n_total = length(results.raw(i).time);
        n_keep = floor(0.80 * n_total);

        t = results.raw(i).time(1:n_keep) / 60;

        ax = nexttile(tl, 2*k - 1);
        CeP = results.raw(i).CeP_trajectory(1:n_keep);
        CeR = results.raw(i).CeR_trajectory(1:n_keep);
        yyaxis(ax, 'left');
        plot(ax, t, CeP, 'b-', 'LineWidth', 1.5);
        ylabel(ax, 'C_{eP} (\mug/mL)');
        set(ax, 'YColor', 'b');
        ylim(ax, [0 8]);
        xlim(ax, [0 max(t)]);
        yyaxis(ax, 'right');
        plot(ax, t, CeR * 1000, 'r-', 'LineWidth', 1.5);
        ylabel(ax, 'C_{eR} (ng/mL)');
        set(ax, 'YColor', [0.8 0 0]);
        ylim(ax, [0 15]);
        xlim(ax, [0 max(t)]);
        title(ax, sprintf('Case %d: effect-site concentrations', results.patient_id(i)));
        grid(ax, 'on');
        if k < m, set(ax, 'XTickLabel', []); end

        ax = nexttile(tl, 2*k);
        hold(ax, 'on');
        hb = plot(ax, t, results.raw(i).bis(1:n_keep), 'k-', 'LineWidth', 1.5);

        % ENDPOINTS
        hE = gobjects(0); hB = gobjects(0);
        if isfield(results.raw(i), 'E0_trajectory') && ~isempty(results.raw(i).E0_trajectory)
            hE = plot(ax, t, results.raw(i).E0_trajectory(1:n_keep), ...
                'Color', [0.6 0.6 0.6], 'LineStyle', '-.', 'LineWidth', 1.2);
        end
        if isfield(results.raw(i), 'BISmin_trajectory') && ~isempty(results.raw(i).BISmin_trajectory)
            hB = plot(ax, t, results.raw(i).BISmin_trajectory(1:n_keep), ...
                'Color', [0.4 0.4 0.4], 'LineStyle', '-.', 'LineWidth', 1.2);
        end

        % 4D MODEL PREDICTIONS (Van and Greco)
        hV = plot(ax, t, results.raw(i).pred_van(1:n_keep), 'b-', 'LineWidth', 2.0);
        hG = plot(ax, t, results.raw(i).pred_gre(1:n_keep), 'r--', 'LineWidth', 2.0);

        % TARGET BAND
        yline(ax, cfg.target_band(1), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
        yline(ax, cfg.target_band(2), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);

        ylim(ax, [0 100]);
        xlim(ax, [0 max(t)]);
        ylabel(ax, 'BIS');
        grid(ax, 'on');
        if k < m, set(ax, 'XTickLabel', []); end

        mae_van = results.metrics.van.MAE(i);
        mae_gre = results.metrics.gre.MAE(i);
        title(ax, {sprintf('%s percentile', cases.labels{k}), ...
                   sprintf('MAE Bouillon %.1f, Greco %.1f', mae_van, mae_gre)});

        if k == 1, h = [hb hE hB hV hG]; end
    end

    xlabel(tl, 'Time (min)');

    h = h(isgraphics(h));
    if ~isempty(h)
        lg = legend(h, {'Measured BIS', 'E_0', 'BIS_{min}', '4D Bouillon', '4D Greco'}, ...
                    'Orientation', 'horizontal', 'Box', 'off');
        lg.Layout.Tile = 'north';
    end

    save_figure(f, fig_dir, 'figure_cases_comparison', 'double', 4.5);
end
