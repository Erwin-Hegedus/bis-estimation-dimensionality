function plot_case_study_params(cases, results, cfg, fig_dir)
%PLOT_CASE_STUDY_PARAMS  4D parameter trajectories, Bouillon vs Greco.
%
%   Twelve panels, three cases x four parameters, in a two-column (figure*)
%   slot at \textwidth. Four columns of panels only work at full width: in a
%   single column they are 0.87 in each and the tick labels collide.
%
%   Redundancy is stripped so the panels can be small and still legible: x tick
%   labels and the time label appear on the bottom row only, the legend appears
%   once, and the row label carries the case id.

    m = numel(cases.idx);
    if m == 0, return; end

    names = {'C_{50P} (\mug/mL)', 'C_{50R} (ng/mL)', '\gamma', '\beta / \alpha'};

    f = figure('Color', 'w');
    tl = tiledlayout(f, m, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

    hV = []; hG = [];
    for r = 1:m
        i = cases.idx(r);

        % CUT LAST 15% (emergence phase)
        n_total = length(results.raw(i).time);
        n_keep  = floor(0.80 * n_total);

        XV = results.raw(i).Xhist_van;
        XG = results.raw(i).Xhist_gre;
        t  = results.raw(i).time(1:n_keep) / 60;

        for p = 1:4
            ax = nexttile(tl, (r-1)*4 + p);
            hold(ax, 'on'); grid(ax, 'on');
            if ~isempty(XV)
                hV = plot(ax, t, XV(p, 1:n_keep), 'b-');
            end
            if ~isempty(XG)
                hG = plot(ax, t, XG(p, 1:n_keep), 'r--');
            end
            if ~isempty(t), xlim(ax, [0 max(t)]); end

            if r == 1, title(ax, names{p}); end
            if p == 1, ylabel(ax, sprintf('Case %d', results.patient_id(i))); end

            % x tick labels on the bottom row only
            if r < m, set(ax, 'XTickLabel', []); end
        end
    end

    xlabel(tl, 'Time (min)');

    if ~isempty(hV) && ~isempty(hG)
        lg = legend([hV hG], {'Bouillon', 'Greco'}, 'Orientation', 'horizontal');
        lg.Layout.Tile = 'north';
        lg.Box = 'off';
    end

    % Explicit: two-column slot, 3.3 in tall. The paper is at its page limit, so
    % this is a budget decision, not an aesthetic one.
    save_figure(f, fig_dir, 'figure_4d_params', 'double', 3.3);
end
