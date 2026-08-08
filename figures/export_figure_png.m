function export_figure_png(fig, png_path, kind, height_in)
%EXPORT_FIGURE_PNG  Export at true printed size: 600 dpi PNG plus a vector PDF.
%
%   The figure is styled to the IEEE two-column geometry first, so the file is
%   written at exactly the width it occupies on the page and LaTeX scales it by
%   1.0. Include it WITHOUT a width override, or with width=\columnwidth
%   (single) / width=\textwidth inside figure* (double) -- both are no-ops at
%   the native size, which is the point.
%
%   A vector PDF is written alongside the PNG. pdflatex resolves a bare
%   \includegraphics{name} to .pdf before .png, so the vector copy is picked up
%   automatically: sharper line art, smaller file, no resampling. Delete the PDF
%   if you specifically want the raster.

    if nargin < 3, kind = []; end
    if nargin < 4, height_in = []; end

    apply_figure_style(fig, kind, height_in);

    d = fileparts(png_path);
    if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end

    exportgraphics(fig, png_path, 'Resolution', 600, ...
        'ContentType', 'image', 'BackgroundColor', 'white');

    % exportgraphics crops to the drawn content, so the file is narrower than
    % the canvas. \includegraphics[width=\linewidth] stretches that crop back to
    % the full column, which magnifies every font by the crop ratio. Report the
    % ratio: it is the difference between the point sizes set in
    % apply_figure_style and the point sizes that reach the page.
    [p, n] = fileparts(png_path);
    info = imfinfo(png_path);
    built_w = get(fig, 'Position') * [0 0 1 0]';
    fprintf('    %-38s %5.2f in built -> %5.2f in exported (fonts x%.2f on page)\n', ...
        n, built_w, info.Width/600, built_w / (info.Width/600));

    pdf_path = fullfile(p, [n '.pdf']);
    try
        exportgraphics(fig, pdf_path, 'ContentType', 'vector', ...
            'BackgroundColor', 'white');
    catch err
        warning('export_figure_png:pdf', ...
            'Vector PDF failed for %s (%s); PNG written.', n, err.message);
    end
end
