function save_figure(fig, fig_dir, filename, kind, height_in)
% SAVE_FIGURE - Save as .fig, 600 dpi .png and vector .pdf at true printed size.
%
%   save_figure(fig, dir, name)                 width from the subplot grid
%   save_figure(fig, dir, name, 'single')       one column,  3.487 in
%   save_figure(fig, dir, name, 'double', 3.2)  two columns, explicit height
%
%   Pass 'kind' whenever you know where the figure goes in the paper. The
%   automatic choice only looks at how many panels sit side by side, which is a
%   guess; stating it is better.

    if nargin < 4, kind = []; end
    if nargin < 5, height_in = []; end

    apply_figure_style(fig, kind, height_in);
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
    savefig(fig, fullfile(fig_dir, [filename '.fig']));
    export_figure_png(fig, fullfile(fig_dir, [filename '.png']), kind, height_in);
    fprintf('  Saved: %s.fig, %s.png, %s.pdf\n', filename, filename, filename);
end
