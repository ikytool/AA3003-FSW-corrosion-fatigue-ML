function export_vector(fig, name, width_mm, height_mm)
%EXPORT_VECTOR Export a figure as a high-quality vector PDF to ../figures.
%   export_vector(fig, 'fig_name', width_mm, height_mm)
%   Widths for International Journal of Fatigue:
%     single column 90 mm, 1.5 column 140 mm, full width 190 mm.

if nargin < 3, width_mm  = 90; end
if nargin < 4, height_mm = 70; end

set(fig, 'Units', 'centimeters');
pos = get(fig, 'Position');
set(fig, 'Position', [pos(1) pos(2) width_mm/10 height_mm/10]);

outdir = fullfile(fileparts(mfilename('fullpath')), '..', 'figures');
if ~exist(outdir, 'dir'), mkdir(outdir); end

outfile = fullfile(outdir, [char(name) '.pdf']);
exportgraphics(fig, outfile, 'ContentType', 'vector');
fprintf('Exported %s\n', outfile);
end
