function apply_figure_style(fig, kind, height_in)
%APPLY_FIGURE_STYLE  Journal styling at TRUE PRINTED SIZE for the JBHI two-column
%   layout.
%
%   apply_figure_style(fig)                  width chosen from the subplot grid
%   apply_figure_style(fig, 'single')        one column,  3.487 in
%   apply_figure_style(fig, 'double')        two columns, 7.140 in
%   apply_figure_style(fig, kind, height)    override height, inches
%
%   Why this changed (2026-08-07). The previous version set 20 pt fonts and then
%   sized the canvas at max(1000, 700*ncol) x max(800, 550*nrow) PIXELS, i.e.
%   10-31 inches wide. LaTeX then scaled that down to \columnwidth = 3.49 in, so
%   every font was divided by 4-9x and landed on the page at 2-5 pt:
%
%       figure2_mae_boxplots  13.8in -> 3.49in  = 0.25x  -> 20 pt printed as 5.1
%       figure7_combined_ident 16.0in -> 3.14in = 0.20x  -> 20 pt printed as 3.9
%       figure_4d_params      31.1in -> 3.49in  = 0.11x  -> 20 pt printed as 2.2
%
%   Bigger numbers in the .m file made the printed text SMALLER, because a
%   larger canvas meant a harder downscale. The fix is to build the figure at
%   the size it is printed at, so LaTeX scales by 1.0 and a 9 pt font is 9 pt.
%
%   Geometry is measured from ieeecolor.cls on US letter: paperwidth 8.5in,
%   side margins 0.680in, \columnsep 1pc -> \textwidth 7.140in,
%   \columnwidth 3.487in. Font sizes below are true points at final size;
%   IEEEtran captions are 8 pt, so 8 pt ticks match the caption and 9 pt labels
%   read one step above it. Do not go below 7 pt.
%
%   Revision 1 (2026-08-07). Reviewer 1: "I would enlarge figures
%   significantly; in particular, labels are really hardly readable." Two things
%   were wrong, and only one of them was in this file.
%
%     1. Every float in manuscript_jbhi_revised.tex is a single-column
%        \begin{figure} at \linewidth = 3.487 in, but the wide figures were
%        being built at TEXTWIDTH and then squeezed to 0.58x on the page. The
%        default is now 'single'; 'double' is only correct for a figure* float,
%        and nothing in this manuscript is one.
%     2. Body text is 9 pt in this class. Ticks now match it and labels sit one
%        step above, so nothing in a panel is smaller than the running text.
%
%   Panels get taller with the fonts: a 10 pt label on a 1.2 in panel is not
%   readable either, it just overlaps.

    % ---- geometry, measured from ieeecolor.cls -------------------------
    COLUMNWIDTH = 3.487;    % in
    TEXTWIDTH   = 7.140;    % in

    % ---- fonts, TRUE POINTS at final printed size ----------------------
    FS_TICK   = 9;
    FS_LABEL  = 10;
    FS_TITLE  = 10;
    FS_LEG    = 9;
    FS_SUPER  = 11;

    % ---- strokes -------------------------------------------------------
    LW_DATA = 1.4;
    LW_AXIS = 0.7;

    if isappdata(fig, 'styled'), return; end
    setappdata(fig, 'styled', true);

    % ---- count the subplot grid ---------------------------------------
    ax = findall(fig, 'Type', 'axes');
    nr = 1; nc = 1;
    if ~isempty(ax)
        pos = cell2mat(get(ax, {'Position'}));
        if size(pos,1) > 1
            nc = numel(uniquetol(pos(:,1), 0.02));
            nr = numel(uniquetol(pos(:,2), 0.02));
        end
    end

    % ---- pick the printed width ---------------------------------------
    % Single column is the default because every float in the manuscript is a
    % single-column one. Guessing 'double' from the panel count only produced
    % figures that LaTeX then halved.
    if nargin < 2 || isempty(kind), kind = 'single'; end
    switch lower(kind)
        case {'single','column','1col'}, W = COLUMNWIDTH;
        case {'double','full','2col','wide'}, W = TEXTWIDTH;
        otherwise, error('apply_figure_style:kind', 'kind must be single or double');
    end

    % A panel narrower than about 1.4 in cannot carry a 9 pt tick label and a
    % 10 pt axis label at the same time. Say so rather than emitting a figure
    % whose text silently overlaps.
    if W / max(nc,1) < 1.4
        warning('apply_figure_style:narrow', ...
            ['%d panels across %.2f in leaves %.2f in each; labels will ' ...
             'collide. Restructure the grid.'], nc, W, W/max(nc,1));
    end

    % ---- pick the printed height --------------------------------------
    % Panels are near square. The 4:3 aspect used before came from the days of
    % wide canvases; in one column a squat panel just wastes the tick labels it
    % is too short to separate.
    if nargin >= 3 && ~isempty(height_in)
        H = height_in;
    else
        H = (W / max(nc,1)) * 0.95 * nr + 0.45;   % +0.45 for the sgtitle/xlabel band
        H = min(max(H, 1.8), 8.5);
    end

    % ---- size the canvas in INCHES, not pixels ------------------------
    set(fig, 'Color', 'w', 'Units', 'inches');
    set(fig, 'Position', [1 1 W H]);
    set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0 0 W H], 'PaperSize', [W H]);

    tl = findall(fig, 'Type', 'tiledlayout');
    for k = 1:numel(tl)
        tl(k).TileSpacing = 'compact';
        tl(k).Padding     = 'compact';
    end

    for t = {'line','errorbar','stair','constantline','functionline'}
        set(findall(fig, 'Type', t{1}), 'LineWidth', LW_DATA);
    end

    % Free-standing annotations first: findall('Type','text') also returns the
    % axis labels and titles, so doing this AFTER the loop below would silently
    % clobber them back down to the legend size.
    set(findall(fig, 'Type', 'text'), 'FontSize', FS_LEG);

    owned = gobjects(0);
    for k = 1:numel(ax)
        a = ax(k);
        set(a, 'FontSize', FS_TICK, 'LineWidth', LW_AXIS);
        set(a, 'GridAlpha', 0.12, 'Box', 'on', 'TickDir', 'out', 'TickLength', [0.015 0.015]);
        if ~isempty(a.Title),  set(a.Title,  'FontSize', FS_TITLE, 'FontWeight', 'bold'); end
        if ~isempty(a.XLabel), set(a.XLabel, 'FontSize', FS_LABEL); end
        if ~isempty(a.YLabel), set(a.YLabel, 'FontSize', FS_LABEL); end
        if ~isempty(a.ZLabel), set(a.ZLabel, 'FontSize', FS_LABEL); end
        owned = [owned; a.Title; a.XLabel; a.YLabel]; %#ok<AGROW>

        % Thin the ticks rather than restating them: 4 columns of panels leaves
        % no room for six labelled ticks, and forcing linspace() over ylim
        % produces unround numbers. Keep MATLAB's choice, drop every other one.
        %
        % How many survive depends on how wide the panel actually is. Left at a
        % fixed count, MATLAB rotates the labels it cannot fit, and a column of
        % 45-degree numbers is the thing the reviewer called unreadable.
        try
            w_in = a.Position(3) * W;
            h_in = a.Position(4) * H;
        catch
            w_in = W; h_in = H;
        end
        thin_ticks(a, 'Y', tick_budget(h_in, 0.42));
        thin_ticks(a, 'X', tick_budget(w_in, 0.62));
    end

    % tiledlayout's own xlabel/ylabel/title are not owned by any axes, so they
    % were caught by the blanket set above. Put them back at label size.
    for k = 1:numel(tl)
        for pr = {'XLabel','YLabel','Title'}
            h = tl(k).(pr{1});
            if ~isempty(h) && ~ismember(h, owned), set(h, 'FontSize', FS_LABEL); end
        end
    end

    % NB: the old version set a.LooseInset = a.TightInset here. That fights
    % exportgraphics' own cropping and clipped rotated y-labels on the tall
    % grids; exportgraphics crops on its own.

    set(findall(fig, 'Type', 'legend'), 'FontSize', FS_LEG);

    cb = findall(fig, 'Type', 'colorbar');
    for k = 1:numel(cb)
        set(cb(k), 'FontSize', FS_LEG, 'LineWidth', LW_AXIS);
        if ~isempty(cb(k).Label), set(cb(k).Label, 'FontSize', FS_LABEL); end
    end

    sg = findall(fig, 'Type', 'text', '-and', 'Tag', 'suptitle');
    set(sg, 'FontSize', FS_SUPER, 'FontWeight', 'bold');
end


function n = tick_budget(extent_in, per_label_in)
%TICK_BUDGET  How many labelled ticks fit along an axis of this length.
%   per_label_in is the width a numeric label plus its breathing room needs at
%   9 pt: about 0.62 in horizontally, 0.42 in for the line spacing vertically.
    n = max(2, min(5, floor(extent_in / per_label_in)));
end


function thin_ticks(a, axname, maxn)
%THIN_TICKS  Keep MATLAB's tick values, drop every other one until <= maxn.
%   Preserves round numbers, unlike forcing linspace() over the limits.
%
%   Only touches axes where BOTH the ticks and their labels are automatic.
%   Thinning a categorical axis or one with hand-written labels deletes
%   categories -- that is what silently dropped three of the six model names
%   from figure2_mae_boxplots on the first pass.
    try
        if strcmpi(get(a, [axname 'Scale']), 'log'), return; end
        if ~strcmpi(get(a, [axname 'TickMode']),      'auto'), return; end
        if ~strcmpi(get(a, [axname 'TickLabelMode']), 'auto'), return; end
        ruler = get(a, [axname 'Axis']);
        if ~isa(ruler, 'matlab.graphics.axis.decorator.NumericRuler'), return; end
        tv = get(a, [axname 'Tick']);
        while numel(tv) > maxn
            tv = tv(1:2:end);
        end
        set(a, [axname 'Tick'], tv);
        % Now that they fit, stop MATLAB from tilting them. Auto-rotation is
        % what it does when the labels are too dense; the thinning above is the
        % fix, the rotation was the symptom.
        set(a, [axname 'TickLabelRotation'], 0);
    catch
    end
end
