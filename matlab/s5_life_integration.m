%S5_LIFE_INTEGRATION Unit recovery, E647 validity, closure checks and
%   conditional life scenarios.
%
%   Stage A - geometry and K-solution. SENT specimen, from the
%   laboratory record of the source campaign: width w = 60 mm,
%   thickness B = 1.8 mm, initial crack a0 = 4 mm, Pmax = 3.5 kN,
%   R = 0.1. K = P/(B*sqrt(w)) * f(a/w) with the ASTM E647 SENT form
%   (Eqs. 1-2 of the source paper). Consistency: dK(a0) = 3.755
%   matches the stated range start of 3.75 MPa*sqrt(m).
%
%   Stage B - recovery of the physical rate unit. The published rate
%   axis label is ambiguous; the factor u between plotted units and
%   mm/cycle is RECOVERED, not assumed: every a-N curve is
%   differentiated with the E647 seven-point incremental polynomial,
%   giving physical rates in mm/cycle, which are compared with the
%   plotted rate curve of the same specimen at the same dK. The median
%   factor across all eight specimens (with spread) converts plotted
%   units to mm/cycle. No condition is consumed by calibration.
%
%   Stage C - ASTM E647 validity. The remaining-ligament criterion
%   w - a >= (4/pi) (Kmax/sigma_ys)^2 is evaluated with the lowest
%   zone elastic limit of the joint (Table 5 of the source paper);
%   integration never leaves the verified domain.
%
%   Stage D - internal closure checks (integrate each refit Paris law
%   back to its own a-N curve; NOT an independent validation) and
%   conditional interpolation scenarios at 45/65 degC (PCHIP rate
%   interpolation; model-form sensitivity is quantified in s6).
%
%   Outputs: figures/fig_aN_validation.pdf, fig_predictions_45_65.pdf,
%            data/derived/geometry_calibration.csv, unit_recovery.csv,
%            e647_validity.csv, life_validation.csv,
%            life_predictions_45_65.csv

clear; close all;
fig_defaults();
S = source_data();

here    = fileparts(mfilename('fullpath'));
datadir = fullfile(here, '..', 'data', 'digitized');
dervdir = fullfile(here, '..', 'data', 'derived');

B_mm   = 1.8;
w_mm   = 60;                  % laboratory record (author-provided)
dP_N   = 0.9 * 3500;          % R = 0.1
sig_ys = 95.3;                % MPa, lowest zone elastic limit (HAZa)
fgeo   = @(aw) sqrt(aw) .* (1.99 - 0.41*aw + 18.7*aw.^2 - 38.48*aw.^3 + 53.85*aw.^4);
dK_fun = @(a, w) dP_N ./ (B_mm .* sqrt(w)) .* fgeo(a ./ w) ./ sqrt(1000); % MPa*sqrt(m)

fprintf('Geometry: w = %.0f mm, B = %.1f mm; dK(a0=4mm) = %.3f MPa*sqrt(m)\n', ...
    w_mm, B_mm, dK_fun(4, w_mm));

% ---------------- Stage B: recover the physical rate unit ------------
pairs = { ...
 'BM',        'fig10_bm.csv',       'fig12_bm.csv';
 'FSW',       'fig10_fsw.csv',      'fig12_fsw.csv';
 'FSW.C 25C', 'fig10_fswc_25c.csv', 'fig12_fswc_25c.csv';
 'FSW.C 55C', 'fig10_fswc_55c.csv', 'fig12_fswc_55c.csv';
 'FSW.C 85C', 'fig10_fswc_85c.csv', 'fig12_fswc_85c.csv';
 'FSW.I 25C', 'fig13_fswi_25c.csv', 'fig14_fswi_25c.csv';
 'FSW.I 55C', 'fig13_fswi_55c.csv', 'fig14_fswi_55c.csv';
 'FSW.I 85C', 'fig13_fswi_85c.csv', 'fig14_fswi_85c.csv'};

urows = [];
for i = 1:size(pairs, 1)
    M = sortrows(readtable(fullfile(datadir, pairs{i,2})), 'N');
    Rt = readtable(fullfile(datadir, pairs{i,3}));
    pR = polyfit(log10(Rt.dK), log10(Rt.dadN), 1);   % plotted-rate Paris fit
    % E647 seven-point incremental polynomial on the a-N curve
    % (duplicate N readings from dense curves are merged first)
    Nu = unique(M.N);
    au = zeros(size(Nu));
    for q = 1:numel(Nu)
        au(q) = mean(M.a(M.N == Nu(q)));
    end
    [av, vv] = deal([]);
    np = numel(au);
    for k = 4:np-3
        idx = k-3:k+3;
        sN = (Nu(idx) - Nu(k)) / 1e5;                 % centred and scaled
        if max(sN) - min(sN) <= 0, continue; end
        c = polyfit(sN, au(idx), 2);
        v = c(2) / 1e5;                               % da/dN, mm/cycle
        if v > 0
            av = [av; au(k)]; vv = [vv; v]; %#ok<AGROW>
        end
    end
    dKv = dK_fun(av, w_mm);
    ok = dKv >= min(Rt.dK) & dKv <= max(Rt.dK);
    fac = vv(ok) ./ 10.^polyval(pR, log10(dKv(ok)));
    urows = [urows; {pairs{i,1}, numel(fac), median(fac), ...
        prctile(fac,25), prctile(fac,75)}]; %#ok<AGROW>
end
Tu = cell2table(urows, 'VariableNames', ...
    {'condition','n_pts','u_median','u_q25','u_q75'});
u = median(Tu.u_median);
fprintf('Recovered unit factor u = %.0f (plotted -> mm/cycle);\n', u);
fprintf('  per-curve medians %.0f to %.0f\n', ...
    min(Tu.u_median), max(Tu.u_median));
writetable(Tu, fullfile(dervdir, 'unit_recovery.csv'));
disp(Tu);

% ---------------- Stage C: ASTM E647 validity domain -----------------
agrid = (4:0.05:28)';
dKg   = dK_fun(agrid, w_mm);
Kmax  = dKg / (1 - S.test.R);
req   = (4/pi) * (Kmax / sig_ys).^2 * 1000;    % mm
lig   = w_mm - agrid;
valid = lig >= req;
iv = find(~valid, 1);
if isempty(iv), a_valid = agrid(end); else, a_valid = agrid(iv - 1); end
dK_valid = dK_fun(a_valid, w_mm);
fprintf('E647 ligament criterion (sigma_ys = %.1f MPa): valid to a = %.2f mm (dK = %.2f)\n', ...
    sig_ys, a_valid, dK_valid);
ii = unique([find(agrid==4), find(abs(agrid-12)<0.026,1), ...
      find(abs(agrid-18)<0.026,1), find(abs(agrid-a_valid)<0.026,1), ...
      find(abs(agrid-24)<0.026,1)]);
Te = table(agrid(ii), dKg(ii), Kmax(ii), agrid(ii)/w_mm, lig(ii), req(ii), valid(ii), ...
    'VariableNames', {'a_mm','dK','Kmax','a_over_w','ligament_mm','required_mm','E647_ok'});
writetable(Te, fullfile(dervdir, 'e647_validity.csv'));
disp(Te);

af_int = min(a_valid, 24);     % integration cap: E647-valid and inside data

writetable(table(w_mm, u, dK_fun(4, w_mm), a_valid, dK_valid, ...
    'VariableNames', {'w_mm','unit_factor','dK_a0','a_valid_mm','dK_valid'}), ...
    fullfile(dervdir, 'geometry_calibration.csv'));

% ---------------- Stage D1: internal closure checks ------------------
Tq   = readtable(fullfile(dervdir, 'paris_refit.csv'));
pick = @(lbl) Tq(strcmp(Tq.label, lbl), :);
conds = { ...
 'BM',        'fig10_bm.csv';       'FSW',       'fig10_fsw.csv';
 'FSW.C 25C', 'fig10_fswc_25c.csv'; 'FSW.C 55C', 'fig10_fswc_55c.csv';
 'FSW.C 85C', 'fig10_fswc_85c.csv';
 'FSW.I 25C', 'fig13_fswi_25c.csv'; 'FSW.I 55C', 'fig13_fswi_55c.csv';
 'FSW.I 85C', 'fig13_fswi_85c.csv'};

fig = figure('Visible', 'off');
tiledlayout(2, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
rows = [];
for i = 1:size(conds, 1)
    P = pick(conds{i,1});
    M = sortrows(readtable(fullfile(datadir, conds{i,2})), 'N');
    a0 = max(4.0, min(M.a));
    af = min(af_int, max(M.a));
    [aa, NN] = integrate_life(P.C * u, P.m, a0, af, w_mm, dK_fun);
    [au, iu] = unique(M.a);
    Nf_meas = interp1(au, M.N(iu), af);
    Nf_pred = max(NN);
    err = 100 * (Nf_pred - Nf_meas) / Nf_meas;
    rows = [rows; {conds{i,1}, af, Nf_meas, Nf_pred, err}]; %#ok<AGROW>

    nexttile; hold on;
    plot(M.N(M.a <= af) / 1e5, M.a(M.a <= af), '.', 'MarkerSize', 4, ...
        'DisplayName', 'measured');
    plot(NN / 1e5, aa, '-', 'DisplayName', 'integrated');
    xlabel('N (10^5 cycles)'); ylabel('a (mm)');
    title(sprintf('%s (%+.0f%%)', conds{i,1}, err), 'FontWeight', 'normal');
    if i == 1, legend('Location', 'northwest'); end
end
export_vector(fig, 'fig_aN_validation', 180, 104);

Tv = cell2table(rows, 'VariableNames', ...
    {'condition','a_end_mm','N_measured','N_integrated','err_pct'});
writetable(Tv, fullfile(dervdir, 'life_validation.csv'));
disp(Tv);
fprintf('Mean |closure error| = %.1f%%\n', mean(abs(Tv.err_pct)));

% ---------------- Stage D2: 45/65 degC interpolation scenarios -------
% Bootstrap (C, m) per anchor temperature from the deduplicated curves,
% PCHIP rate interpolation in 1/T, integration inside the valid domain.
rng(1);
nboot = 300;
targets = [45 65];
fig = figure('Visible', 'off');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
prows = [];
a0s = 4.6;
for panel = 1:2
    inh = panel - 1;
    if inh == 0, fam = 'fig12_fswc'; else, fam = 'fig14_fswi'; end
    D = cell(1, 3); temps = [25 55 85];
    for t = 1:3
        D{t} = readtable(fullfile(datadir, sprintf('%s_%dc.csv', fam, temps(t))));
    end
    nexttile; hold on;
    cols = get(gca, 'ColorOrder');
    for it = 1:numel(targets)
        Tt = targets(it);
        Nf = zeros(nboot, 1);
        agrid2 = linspace(a0s, af_int, 200);
        curves = zeros(nboot, 200);
        for b = 1:nboot
            Cf = zeros(1,3); mf = zeros(1,3);
            for t = 1:3
                d = D{t}; k = randi(height(d), height(d), 1);
                p = polyfit(log10(d.dK(k)), log10(d.dadN(k)), 1);
                mf(t) = p(1); Cf(t) = p(2);
            end
            [aa, NN] = integrate_life_v(@(dk) ...
                u * rate_at_T(dk, Cf, mf, S.TK(temps), S.TK(Tt)), ...
                a0s, af_int, w_mm, dK_fun);
            Nf(b) = max(NN);
            curves(b, :) = interp1(aa, NN, agrid2, 'linear', 'extrap');
        end
        med = median(curves, 1); lo = prctile(curves, 2.5, 1); hi = prctile(curves, 97.5, 1);
        cc = cols(it, :);
        fill([lo, fliplr(hi)] / 1e5, [agrid2, fliplr(agrid2)], cc, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(med / 1e5, agrid2, '-', 'Color', cc, ...
            'DisplayName', sprintf('%d ^{\\circ}C scenario', Tt));
        prows = [prows; {inh, Tt, median(Nf), prctile(Nf, 2.5), prctile(Nf, 97.5)}]; %#ok<AGROW>
    end
    % measured neighbours for context (truncated at the same af)
    nb = {'fig10_fswc_25c.csv', 'fig13_fswi_25c.csv';
          'fig10_fswc_55c.csv', 'fig13_fswi_55c.csv';
          'fig10_fswc_85c.csv', 'fig13_fswi_85c.csv'};
    mkc = [0.4 0.4 0.4];
    for t = 1:3
        M = sortrows(readtable(fullfile(datadir, nb{t, inh+1})), 'N');
        M = M(M.a <= af_int, :);
        vis = 'off'; if t == 1, vis = 'on'; end
        plot(M.N / 1e5, M.a, ':', 'Color', mkc, 'LineWidth', 0.8, ...
            'HandleVisibility', vis, ...
            'DisplayName', 'measured 25/55/85 ^{\circ}C');
    end
    xlabel('N (10^5 cycles)'); ylabel('a (mm)');
    if inh == 0, title('(a) FSW.C, no inhibitor', 'FontWeight', 'normal');
    else,        title('(b) FSW.I, 10% ethylene glycol', 'FontWeight', 'normal'); end
    legend('Location', 'north');
end
export_vector(fig, 'fig_predictions_45_65', 180, 81);

Tp = cell2table(prows, 'VariableNames', {'inhibited','T_C','Nf_median','Nf_lo95','Nf_hi95'});
writetable(Tp, fullfile(dervdir, 'life_predictions_45_65.csv'));
disp(Tp);

% ======================================================================
function [aa, NN] = integrate_life(C_mm, m, a0, af, w, dK_fun)
    aa = linspace(a0, af, 400)';
    v  = max(C_mm .* dK_fun(aa, w) .^ m, 1e-12);   % mm/cycle
    NN = cumtrapz(aa, 1 ./ v);
end

function [aa, NN] = integrate_life_v(vfun, a0, af, w, dK_fun)
    aa = linspace(a0, af, 400)';
    v  = max(vfun(dK_fun(aa, w)), 1e-12);          % mm/cycle
    NN = cumtrapz(aa, 1 ./ v);
end

function v = rate_at_T(dk, Cf, mf, TKanchors, TKq)
    % per-dK shape-preserving interpolation of ln v in 1/T through the
    % three anchor temperatures (exact at anchors, no overshoot between
    % them; the 25-55 degC plateau + sharp 85 degC rise is strongly
    % non-Arrhenius, so a single straight line would overshoot)
    [x, i] = sort(1 ./ TKanchors);
    v = zeros(size(dk));
    for k = 1:numel(dk)
        lv = log(10.^Cf .* dk(k) .^ mf);
        lv = lv(i);
        v(k) = exp(pchip(x, lv, 1 / TKq));
    end
end
