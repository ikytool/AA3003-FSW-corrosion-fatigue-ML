%DIGITIZE_ALL Digitize Figs 10/12/13/14 of the JMRT source paper.
%   Input:  ../../data/source_pages/page{08,09,10}.png (5000 px wide
%           renders of the publisher PDF pages).
%   Output: ../../data/digitized/fig*.csv  (per curve)
%           diagnostic overlay PNGs in ./diagnostics for visual QA.
%
%   Uses base MATLAB only (no Image Processing Toolbox): morphology via
%   movmax/movmin, connected components via graph/conncomp.
%
%   Method: axes-frame detection -> major-tick detection (ticks at the
%   frame corners are dropped; C.xticks/C.yticks list INTERIOR majors
%   only) -> per-series colour segmentation (black series get the union
%   of colour masks subtracted, so coloured markers' dark outlines do
%   not leak in) -> marker centroids (scatter figs) or column-median
%   line tracing (a-N figs) -> axis transform. The QA gate against the
%   published Paris fits is applied by s1_paris_baseline.m.

clear; close all;
here    = fileparts(mfilename('fullpath'));
pagedir = fullfile(here, '..', '..', 'data', 'source_pages');
outdir  = fullfile(here, '..', '..', 'data', 'digitized');
diagdir = fullfile(here, 'diagnostics');
if ~exist(outdir, 'dir'),  mkdir(outdir);  end
if ~exist(diagdir, 'dir'), mkdir(diagdir); end

cfg = figure_configs();

for f = 1:numel(cfg)
    C = cfg(f);
    fprintf('\n=== %s ===\n', C.name);
    img  = imread(fullfile(pagedir, C.page));
    crop = img(C.crop(3):C.crop(4), C.crop(1):C.crop(2), :);

    [fr, cal] = calibrate_axes(crop, C);

    % union of all colour (hue-mode) masks, for outline subtraction
    humask = false(size(crop, 1), size(crop, 2));
    for s = 1:numel(C.series)
        if strcmp(C.series(s).mode, 'hue')
            humask = humask | colour_mask(crop, C.series(s));
        end
    end
    humask_d = dilate(humask, 3);

    overlay = crop;
    for s = 1:numel(C.series)
        SS = C.series(s);
        mask = colour_mask(crop, SS);
        if ~strcmp(SS.mode, 'hue')
            mask = mask & ~humask_d;
        end
        mask = restrict_to_frame(mask, fr, C.excl);
        if C.scatterfig
            ar = C.arange;
            if strcmp(SS.mode, 'gray'), ar(1) = 160; end  % reject gridline dashes
            [px, py] = marker_centroids(mask, ar);
        else
            [px, py] = trace_line(mask, fr);
        end
        [xd, yd] = px2data(px, py, cal);
        keep = xd >= cal.xlim(1) & xd <= cal.xlim(2) & ...
               yd >= cal.ylim(1) & yd <= cal.ylim(2);
        xd = xd(keep); yd = yd(keep); px = px(keep); py = py(keep);
        [xd, i] = sort(xd); yd = yd(i); px = px(i); py = py(i);

        T = table(xd(:), yd(:), 'VariableNames', C.cols);
        writetable(T, fullfile(outdir, [SS.file '.csv']));
        fprintf('  %-18s %3d points\n', SS.file, numel(xd));

        overlay = draw_points(overlay, px, py);
    end

    k = ceil(size(overlay, 2) / 1300);
    imwrite(overlay(1:k:end, 1:k:end, :), ...
        fullfile(diagdir, [C.name '_overlay.png']));
end
fprintf('\nDiagnostics in matlab/digitize/diagnostics - inspect visually.\n');

% ======================================================================
function cfg = figure_configs()
    % xticks/yticks: INTERIOR major ticks only (corner ticks dropped)

    c12.name = 'fig12'; c12.page = 'page09.png';
    c12.crop = [2600 4700 350 1900];
    c12.scatterfig = true;  c12.arange = [70 8000];
    c12.cols = {'dK', 'dadN'};
    c12.xlog = true;  c12.ylog = true;
    c12.xticks = [4 6 8 10 12 14 16];
    c12.yticks = [1e-6 1e-7 1e-8];   % TOP to BOTTOM
    c12.xlim = [2.85 18.6]; c12.ylim = [8e-10 4e-6];
    c12.excl = [0.00 0.33 0.00 0.32;
                0.47 1.00 0.60 1.00];
    c12.series = [ ...
        srs('fig12_bm',        'dark',  [],        [0    0.35], [0    0.45]);
        srs('fig12_fsw',       'hue',   [265 315], [0.25 1],    [0.25 1.00]);
        srs('fig12_fswc_25c',  'hue',   [205 255], [0.35 1],    [0.25 1.00]);
        srs('fig12_fswc_55c',  'hue',   [ 95 165], [0.30 1],    [0.25 1.00]);
        srs('fig12_fswc_85c',  'hue',   [345  15], [0.40 1],    [0.30 1.00])];

    c14.name = 'fig14'; c14.page = 'page10.png';
    c14.crop = [2600 4700 350 1900];
    c14.scatterfig = true;  c14.arange = [70 8000];
    c14.cols = {'dK', 'dadN'};
    c14.xlog = true;  c14.ylog = true;
    c14.xticks = [4 6 8 10 12 14 16];
    c14.yticks = [1e-6 1e-7 1e-8];   % TOP to BOTTOM
    c14.xlim = [2.85 16.2]; c14.ylim = [8e-10 4e-6];
    c14.excl = [0.00 0.40 0.00 0.38;
                0.45 1.00 0.55 1.00];
    c14.series = [ ...
        srs('fig14_fswc_25c',  'hue',   [ 95 165], [0.30 1],    [0.15 0.80]);
        srs('fig14_fswi_25c',  'hue',   [168 205], [0.35 1],    [0.45 1.00]);
        srs('fig14_fswc_55c',  'hue',   [210 255], [0.35 1],    [0.25 1.00]);
        srs('fig14_fswi_55c',  'gray',  [],        [0    0.22], [0.30 0.80]);
        srs('fig14_fswc_85c',  'hue',   [345  15], [0.40 1],    [0.30 1.00]);
        srs('fig14_fswi_85c',  'hue',   [295 345], [0.30 1],    [0.35 1.00])];

    c10.name = 'fig10'; c10.page = 'page08.png';
    c10.crop = [2500 4700 4550 5990];
    c10.scatterfig = false; c10.arange = [120 8000];
    c10.cols = {'N', 'a'};
    c10.xlog = false; c10.ylog = false;
    c10.xticks = [0 5e5 1e6 1.5e6];
    c10.yticks = [30 20 10];         % TOP to BOTTOM
    c10.xlim = [-2e4 2.05e6]; c10.ylim = [0 36];
    c10.excl = [0.00 0.34 0.00 0.30];
    c10.series = [ ...
        srs('fig10_bm',        'dark',  [],        [0    0.35], [0    0.45]);
        srs('fig10_fsw',       'hue',   [265 320], [0.30 1],    [0.30 1.00]);
        srs('fig10_fswc_25c',  'hue',   [215 255], [0.35 1],    [0.25 1.00]);
        srs('fig10_fswc_55c',  'hue',   [ 95 165], [0.30 1],    [0.25 1.00]);
        srs('fig10_fswc_85c',  'hue',   [345  15], [0.40 1],    [0.30 1.00])];

    c13.name = 'fig13'; c13.page = 'page10.png';
    c13.crop = [250 2450 4200 5990];
    c13.scatterfig = false; c13.arange = [120 8000];
    c13.cols = {'N', 'a'};
    c13.xlog = false; c13.ylog = false;
    c13.xticks = [0 5e5 1e6 1.5e6];
    c13.yticks = [30 20 10];         % TOP to BOTTOM
    c13.xlim = [-2e4 1.75e6]; c13.ylim = [0 36];
    c13.excl = [0.00 0.34 0.00 0.40];
    c13.series = [ ...
        srs('fig13_fswc_25c',  'hue',   [135 185], [0.30 1],    [0.75 1.00]);
        srs('fig13_fswi_25c',  'hue',   [ 95 150], [0.30 1],    [0.15 0.74]);
        srs('fig13_fswc_55c',  'hue',   [215 260], [0.35 1],    [0.15 0.85]);
        srs('fig13_fswi_55c',  'gray',  [],        [0    0.22], [0.30 0.90]);
        srs('fig13_fswc_85c',  'hue',   [345  15], [0.40 1],    [0.30 1.00]);
        srs('fig13_fswi_85c',  'hue',   [315 345], [0.35 1],    [0.35 1.00])];

    cfg = [c12 c14 c10 c13];
end

function s = srs(file, mode, hue, sat, val)
    s.file = file; s.mode = mode; s.hue = hue; s.sat = sat; s.val = val;
end

% ======================================================================
function [fr, cal] = calibrate_axes(crop, C)
    g = double(rgb2gray(crop));
    dark = g < 120;
    [H, W] = size(dark);

    rowfrac = sum(dark, 2) / W;
    colfrac = sum(dark, 1) / H;
    rows = find(rowfrac > 0.55);
    cols = find(colfrac > 0.55);
    assert(~isempty(rows) && ~isempty(cols), 'frame not found');
    fr.y1 = rows(1); fr.y2 = rows(end);
    fr.x1 = cols(1); fr.x2 = cols(end);
    fprintf('  frame px: x [%d %d], y [%d %d]\n', fr.x1, fr.x2, fr.y1, fr.y2);

    cap = 45;   % max tick stroke length considered (px)

    % x axis: tick strokes rise from the bottom frame line into the plot
    runs = zeros(1, W);
    for x = 1:W
        r = 0;
        for y = fr.y2-3:-1:max(1, fr.y2-3-cap)
            if dark(y, x), r = r + 1; else, break; end
        end
        runs(x) = r;
    end
    xt = major_ticks(runs, fr.x1, fr.x2, 'x');

    % y axis: tick strokes extend right from the left frame line
    runsy = zeros(1, H);
    for y = 1:H
        r = 0;
        for x = fr.x1+3:min(W, fr.x1+3+cap)
            if dark(y, x), r = r + 1; else, break; end
        end
        runsy(y) = r;
    end
    yt = major_ticks(runsy, fr.y1, fr.y2, 'y');

    cal.xlog = C.xlog; cal.ylog = C.ylog;
    cal.xlim = C.xlim; cal.ylim = C.ylim;
    cal.xfit = tickfit(xt, C.xticks, C.xlog, 'x');
    cal.yfit = tickfit(yt, C.yticks, C.ylog, 'y');
end

function pk = major_ticks(runs, lo, hi, label)
    % group contiguous tick-stroke columns; measure each group's length;
    % keep long ("major") ticks; drop groups near the frame corners
    on = runs >= 8;
    d = diff([0 on 0]);
    a = find(d == 1); b = find(d == -1) - 1;
    pos = zeros(1, numel(a)); len = zeros(1, numel(a));
    for i = 1:numel(a)
        idx = a(i):b(i);
        pos(i) = sum(idx .* runs(idx)) / sum(runs(idx));
        len(i) = max(runs(idx));
    end
    ok = pos > lo + 10 & pos < hi - 10;
    pos = pos(ok); len = len(ok);
    if isempty(pos), pk = []; return; end
    thr = 0.60 * max(len);
    pk = pos(len >= thr);
    fprintf('  %s ticks: %d groups, lengths %s -> %d majors at %s\n', ...
        label, numel(len), mat2str(round(len)), numel(pk), mat2str(round(pk)));
end

function p = tickfit(pix, vals, uselog, label)
    v = vals; if uselog, v = log10(vals); end
    assert(numel(pix) >= 2, '%s axis: too few major ticks', label);
    if numel(pix) == numel(v)
        p = polyfit(pix, v, 1);
    else
        error('%s axis: found %d majors, expected %d - adjust config', ...
            label, numel(pix), numel(v));
    end
    % report fit residuals in px for sanity
    resid = (v - polyval(p, pix)) / p(1);
    fprintf('  %s calibration residuals (px): %s\n', label, mat2str(round(resid, 1)));
end

function [x, y] = px2data(px, py, cal)
    x = polyval(cal.xfit, px);
    y = polyval(cal.yfit, py);
    if cal.xlog, x = 10.^x; end
    if cal.ylog, y = 10.^y; end
end

% ======================================================================
function mask = colour_mask(crop, S)
    hsv = rgb2hsv(crop);
    Hh = hsv(:,:,1) * 360; Ss = hsv(:,:,2); Vv = hsv(:,:,3);
    switch S.mode
        case 'hue'
            if S.hue(1) <= S.hue(2)
                hm = Hh >= S.hue(1) & Hh <= S.hue(2);
            else
                hm = Hh >= S.hue(1) | Hh <= S.hue(2);
            end
            mask = hm & Ss >= S.sat(1) & Ss <= S.sat(2) & ...
                        Vv >= S.val(1) & Vv <= S.val(2);
        case 'dark'
            mask = Ss <= S.sat(2) & Vv <= S.val(2);
        case 'gray'
            mask = Ss <= S.sat(2) & Vv >= S.val(1) & Vv <= S.val(2);
    end
end

function mask = restrict_to_frame(mask, fr, excl)
    m = false(size(mask));
    pad = 30;   % must exceed the longest tick stroke (26 px)
    m(fr.y1+pad:fr.y2-pad, fr.x1+pad:fr.x2-pad) = ...
        mask(fr.y1+pad:fr.y2-pad, fr.x1+pad:fr.x2-pad);
    w = fr.x2 - fr.x1; h = fr.y2 - fr.y1;
    for e = 1:size(excl, 1)
        xa = fr.x1 + round(excl(e,1) * w); xb = fr.x1 + round(excl(e,2) * w);
        ya = fr.y1 + round(excl(e,3) * h); yb = fr.y1 + round(excl(e,4) * h);
        m(max(1,ya):yb, max(1,xa):xb) = false;
    end
    mask = m;
end

% ---- base-MATLAB morphology -----------------------------------------
function m = dilate(m, k)
    m = movmax(movmax(double(m), k, 1), k, 2) > 0;
end
function m = erode(m, k)
    m = movmin(movmin(double(m), k, 1), k, 2) > 0;
end
function m = closeop(m, k)
    m = erode(dilate(m, k), k);
end

function [px, py] = marker_centroids(mask, arange)
    mask = closeop(mask, 5);
    [lbl, n] = concomp(mask);
    px = []; py = [];
    if n == 0, return; end
    [yy, xx] = find(mask);
    id = lbl(sub2ind(size(mask), yy, xx));
    for i = 1:n
        sel = id == i;
        A = nnz(sel);
        w = max(xx(sel)) - min(xx(sel)) + 1;
        h = max(yy(sel)) - min(yy(sel)) + 1;
        aspect = w / h;
        % markers are roughly round; gridline dashes are thin and long
        if A >= arange(1) && A <= arange(2) && aspect >= 0.45 && aspect <= 2.2
            px(end+1) = mean(xx(sel)); %#ok<AGROW>
            py(end+1) = mean(yy(sel)); %#ok<AGROW>
        end
    end
    px = px(:); py = py(:);
end

function [lbl, n] = concomp(mask)
    % 4-connected components via graph/conncomp (no toolbox needed)
    idx = find(mask);
    lbl = zeros(size(mask));
    if isempty(idx), n = 0; return; end
    [H, ~] = size(mask);
    map = zeros(size(mask)); map(idx) = 1:numel(idx);
    % right neighbours
    r = idx(mask(min(idx + H, numel(mask))) & (idx + H <= numel(mask)));
    e1 = [map(r), map(r + H)];
    % down neighbours (avoid wrapping across column ends)
    hasdown = mod(idx, H) ~= 0 & idx + 1 <= numel(mask);
    dn = idx(hasdown & mask(min(idx + 1, numel(mask))));
    e2 = [map(dn), map(dn + 1)];
    G = graph([e1(:,1); e2(:,1)], [e1(:,2); e2(:,2)], [], numel(idx));
    bins = conncomp(G);
    lbl(idx) = bins;
    n = max(bins);
end

function [px, py] = trace_line(mask, fr)
    % column scan for the shallow part + row scan for near-vertical
    % failure tails (row samples only where the curve is locally steep)
    mask = closeop(mask, 3);
    px = []; py = [];
    for x = fr.x1+12:6:fr.x2-12
        ys = find(mask(:, x));
        if ~isempty(ys)
            px(end+1) = x;          %#ok<AGROW>
            py(end+1) = median(ys); %#ok<AGROW>
        end
    end
    for y = fr.y1+12:6:fr.y2-12
        xs = find(mask(y, :));
        if ~isempty(xs) && (max(xs) - min(xs)) <= 60
            px(end+1) = median(xs); %#ok<AGROW>
            py(end+1) = y;          %#ok<AGROW>
        end
    end
    px = px(:); py = py(:);
end

function img = draw_points(img, px, py)
    for i = 1:numel(px)
        cx = round(px(i)); cy = round(py(i));
        xs = max(1, cx-11):min(size(img,2), cx+11);
        ys = max(1, cy-2):min(size(img,1), cy+2);
        img(ys, xs, 1) = 255; img(ys, xs, 2) = 140; img(ys, xs, 3) = 0;
        ys2 = max(1, cy-11):min(size(img,1), cy+11);
        xs2 = max(1, cx-2):min(size(img,2), cx+2);
        img(ys2, xs2, 1) = 255; img(ys2, xs2, 2) = 140; img(ys2, xs2, 3) = 0;
    end
end
