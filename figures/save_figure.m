function save_figure(fig, fig_dir, filename)
% SAVE_FIGURE - Save figure as both .fig and .png
    savefig(fig, fullfile(fig_dir, [filename '.fig']));
    print(fig, fullfile(fig_dir, [filename '.png']), '-dpng', '-r300');
    fprintf('  Saved: %s.fig and %s.png\n', filename, filename);
end
