function apply_figure_style(fig)
%APPLY_FIGURE_STYLE  Journal styling: 20 pt text, 22 pt bold titles, legends
%   at 16, 2 pt lines and axes, no outer padding.

    FS       = 20;
    FS_TITLE = 22;
    FS_LEG   = FS - 4;
    LW       = 2.0;

    if isappdata(fig, 'styled'), return; end
    setappdata(fig, 'styled', true);

    ax = findall(fig, 'Type', 'axes');
    nr = 1; nc = 1;
    if ~isempty(ax)
        pos = cell2mat(get(ax, {'Position'}));
        if size(pos,1) > 1
            nc = numel(uniquetol(pos(:,1), 0.02));
            nr = numel(uniquetol(pos(:,2), 0.02));
        end
    end
    set(fig, 'Color', 'w', 'Position', [50 50 max(1000, 700*nc), max(800, 550*nr)]);

    tl = findall(fig, 'Type', 'tiledlayout');
    for k = 1:numel(tl)
        tl(k).TileSpacing = 'compact';
        tl(k).Padding     = 'tight';
    end

    for t = {'line','errorbar','stair','constantline','functionline'}
        set(findall(fig, 'Type', t{1}), 'LineWidth', LW);
    end

    for k = 1:numel(ax)
        a = ax(k);
        set(a, 'FontSize', FS, 'LineWidth', LW);
        if ~isempty(a.Title),  set(a.Title,  'FontSize', FS_TITLE, 'FontWeight', 'bold'); end
        if ~isempty(a.XLabel), set(a.XLabel, 'FontSize', FS); end
        if ~isempty(a.YLabel), set(a.YLabel, 'FontSize', FS); end
        if ~isempty(a.ZLabel), set(a.ZLabel, 'FontSize', FS); end
        try
            a.LooseInset = a.TightInset;
        catch
        end
    end

    lg = findall(fig, 'Type', 'legend');
    set(lg, 'FontSize', FS_LEG);

    cb = findall(fig, 'Type', 'colorbar');
    for k = 1:numel(cb)
        set(cb(k), 'FontSize', FS_LEG, 'LineWidth', LW);
        if ~isempty(cb(k).Label), set(cb(k).Label, 'FontSize', FS); end
    end

    sg = findall(fig, 'Type', 'text', '-and', 'Tag', 'suptitle');
    set(sg, 'FontSize', FS_TITLE, 'FontWeight', 'bold');
end
