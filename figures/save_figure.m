function save_figure(fig, fig_dir, filename)
% SAVE_FIGURE - Save figure as both .fig and .png at journal styling
    apply_figure_style(fig);
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
    savefig(fig, fullfile(fig_dir, [filename '.fig']));
    export_figure_png(fig, fullfile(fig_dir, [filename '.png']));
    fprintf('  Saved: %s.fig and %s.png\n', filename, filename);
end
