%S0_SETUP_MONTAGE One-figure montage of the experimental setup photos.
%   Assembles the equipment photographs of the two source papers
%   (authors' own published figures, re-used with citation) into a
%   single 2x3 figure: (a) FSW process, (b) FSW tool, (c) fatigue
%   machine with crack monitoring, (d) SENT specimen in the grips,
%   (e) pre-corrosion samples with temperature logger, (f) thermostatic
%   bath. Sources: JMRT Fig. 9 (page08.png) and Colloids Surf. A
%   Figs. 1 and 4 (colloids_p03/p04.png).
%   Produces figures/fig_setup.pdf (180 mm wide).

clear; close all;
fig_defaults();
here    = fileparts(mfilename('fullpath'));
pagedir = fullfile(here, '..', 'data', 'source_pages');

jm = imread(fullfile(pagedir, 'page08.png'));      % JMRT page 8 (Fig. 9)
j3 = imread(fullfile(pagedir, 'ismail_p03.png'));  % JMRT page 3 (Fig. 3)
c3 = imread(fullfile(pagedir, 'colloids_p03.png'));% Colloids page 3 (Fig. 1)

% crops [y1 y2 x1 x2] in page-pixel coordinates (5000 px wide renders)
P = { ...
  c3, [ 480 1720 1090 2860], '(a) FSW process';
  c3, [ 480 1680 3170 3900], '(b) FSW tool';
  jm, [4640 5860  320 1500], '(c) fatigue machine';
  jm, [4640 5860 1520 2400], '(d) SENT specimen'};

fig = figure('Visible', 'off');
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:size(P, 1)
    nexttile;
    r = P{i, 2};
    imshow(P{i, 1}(r(1):r(2), r(3):r(4), :));
    title(P{i, 3}, 'FontWeight', 'normal', 'FontSize', 10, ...
        'FontName', 'Times New Roman');
end
% panel (e): JMRT Fig. 3 (immersion setup), portrait image over 2 tiles
nexttile([1 2]);
imshow(j3(2680:4990, 2580:4650, :));
title('(e) pre-corrosion immersion setup', 'FontWeight', 'normal', ...
    'FontSize', 10, 'FontName', 'Times New Roman');

set(fig, 'Units', 'centimeters');
pos = get(fig, 'Position');
set(fig, 'Position', [pos(1) pos(2) 18 10.5]);
outdir = fullfile(here, '..', 'figures');
exportgraphics(fig, fullfile(outdir, 'fig_setup.pdf'), ...
    'ContentType', 'image', 'Resolution', 300);
fprintf('Exported figures/fig_setup.pdf\n');