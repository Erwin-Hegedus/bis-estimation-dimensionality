function export_figure_png(fig, png_path)
%EXPORT_FIGURE_PNG  Styled export at 300 dpi, tight.

    apply_figure_style(fig);
    d = fileparts(png_path);
    if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end

    exportgraphics(fig, png_path, 'Resolution', 300, ...
        'ContentType', 'image', 'BackgroundColor', 'white');
end
