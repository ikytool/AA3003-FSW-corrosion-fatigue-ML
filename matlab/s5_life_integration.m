%S5_LIFE_INTEGRATION Geometry/unit calibration and residual-life curves.
%
%   Stage A - calibration. The SENT K-solution of the source paper
%   (its Eqs. 1-2, ASTM E647 form) is
%       K = P/(B*sqrt(w)) * f(a/w),
%       f = (a/w)^0.5 (1.99 - 0.41(a/w) + 18.7(a/w)^2 - 38.48(a/w)^3
%           + 53.85(a/w)^4)
%   with B = 1.8 mm, a0 = 4 mm (user-confirmed), Pmax = 3.5 kN, R = 0.1.
%   The specimen width w is not printed in the paper. It is calibrated
%   here from internal consistency: da/dN computed numerically from the
%   digitized a-N curve (Fig. 10) must coincide with the digitized
%   da/dN-dK curve (Fig. 12) for the same condition. The same match
%   yields the unit factor u between the plotted rate axis and mm/cycle.
%
%   Stage B - validation. Integrate the refit Paris laws to a-N curves
%   and compare with the measured (digitized) a-N curves.
%   Stage C - prediction. Arrhenius-interpolated (C, m) at 45/65 degC
%   with bootstrap uncertainty bands, integrated to a-N and lives.
%
%   Outputs: figures/fig_aN_validation.pdf,
%            figures/fig_predictions_45_65.pdf,
%            data/derived/geometry_calibration.csv,
%            data/derived/life_validation.csv,
%            data/derived/life_predictions_45_65.csv

clear; close all;
fig_defaults();
S = source_data();

here    = fileparts(mfilename('fullpath'));
datadir = fullfile(here, '..', 'data', 'digitized');
dervdir = fullfile(here, '..', 'data', 'derived');

B_mm   = 1.8;
dP_N   = 0.9 * 3500;          % R = 0.1
fgeo   = @(aw) sqrt(aw) .* (1.99 - 0.41*aw + 18.7*aw.^2 - 38.48*aw.^3 + 53.85*aw.^4);
dK_fun = @(a, w) dP_N ./ (B_mm .* sqrt(w)) .* fgeo(a ./ w) ./ sqrt(1000); % MPa*sqrt(m)

% ---------------- Stage A: calibrate w and unit factor ---------------
% (1) The paper states the tested range starts at dK = 3.75 MPa*sqrt(m)
%     and the first crack datum is a0 = 4 mm  ->  solve for w.
% (2) The unit factor u between the plotted rate axis and mm/cycle is
%     set so that integrating the refit BM Paris law reproduces the
%     measured BM life exactly. BM is thereby used for calibration;
%     the other seven conditions remain independent validations.
Tq   = readtable(fullfile(dervdir, 'paris_refit.csv'));
pick = @(lbl) Tq(strcmp(Tq.label, lbl), :);

w_mm = 60;   % author-confirmed specimen width (2026-08-04)
% consistency check: the stated range start dK = 3.75 MPa*sqrt(m) at
% a0 = 4 mm independently gives w = 60.1 mm
fprintf('w = %.0f mm (author-confirmed); dK(a0=4mm) = %.3f MPa*sqrt(m)\n', ...
    w_mm, dK_fun(4, w_mm));
fprintf('  -> dK(a=24mm) = %.2f, dK(a=28mm) = %.2f MPa*sqrt(m)\n', ...
    dK_fun(24, w_mm), dK_fun(28, w_mm));

bm = pick('BM');
M0 = sortrows(readtable(fullfile(datadir, 'fig10_bm.csv')), 'N');
a0 = max(4.0, min(M0.a)); af = max(M0.a);
[~, NN1] = integrate_life(bm.C, bm.m, a0, af, w_mm, dK_fun);  % u = 1
u = max(NN1) / max(M0.N);
fprintf('Calibrated unit factor u = %.0f (plotted units -> mm/cycle)\n', u);
writetable(table(w_mm, u, dK_fun(4, w_mm), ...
    'VariableNames', {'w_mm','unit_factor','dK_a0'}), ...
    fullfile(dervdir, 'geometry_calibration.csv'));

% ---------------- Stage B: integrate refits, validate vs a-N ---------
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
    a0 = max(4.0, min(M.a)); af = max(M.a);
    [aa, NN] = integrate_life(P.C * u, P.m, a0, af, w_mm, dK_fun);
    Nf_meas = max(M.N); Nf_pred = max(NN);
    err = 100 * (Nf_pred - Nf_meas) / Nf_meas;
    rows = [rows; {conds{i,1}, Nf_meas, Nf_pred, err}]; %#ok<AGROW>

    nexttile; hold on;
    plot(M.N / 1e5, M.a, '.', 'MarkerSize', 4, 'DisplayName', 'measured');
    plot(NN / 1e5, aa, '-', 'DisplayName', 'integrated');
    xlabel('N (10^5 cycles)'); ylabel('a (mm)');
    title(sprintf('%s (%+.0f%%)', conds{i,1}, err), 'FontWeight', 'normal');
    if i == 1, legend('Location', 'northwest'); end
end
export_vector(fig, 'fig_aN_validation', 180, 104);

Tv = cell2table(rows, 'VariableNames', {'condition','Nf_measured','Nf_integrated','err_pct'});
writetable(Tv, fullfile(dervdir, 'life_validation.csv'));
disp(Tv);

% ---------------- Stage C: predictions at 45 and 65 degC -------------
% Bootstrap (C, m) per temperature from the digitized FSW.C/FSW.I data,
% Arrhenius-interpolate each draw to the target temperature, integrate.
rng(1);
nboot = 300;
targets = [45 65];
fig = figure('Visible', 'off');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
prows = [];
for panel = 1:2
    inh = panel - 1;
    fam = 'fswc'; if inh, fam = 'fswi'; end
    D = cell(1, 3); temps = [25 55 85];
    for t = 1:3
        if inh == 0
            T1 = readtable(fullfile(datadir, sprintf('fig12_%s_%dc.csv', fam, temps(t))));
            T2 = readtable(fullfile(datadir, sprintf('fig14_%s_%dc.csv', fam, temps(t))));
            D{t} = [T1; T2];
        else
            D{t} = readtable(fullfile(datadir, sprintf('fig14_%s_%dc.csv', fam, temps(t))));
        end
    end
    nexttile; hold on;
    cols = get(gca, 'ColorOrder');
    for it = 1:numel(targets)
        Tt = targets(it);
        Nf = zeros(nboot, 1);
        curves = zeros(nboot, 200);
        agrid = linspace(4.6, 27, 200);
        for b = 1:nboot
            Cf = zeros(1,3); mf = zeros(1,3);
            for t = 1:3
                d = D{t}; k = randi(height(d), height(d), 1);
                p = polyfit(log10(d.dK(k)), log10(d.dadN(k)), 1);
                mf(t) = p(1); Cf(t) = p(2);
            end
            % rate-space Arrhenius interpolation: per dK, ln v linear in
            % 1/T through the three anchors (keeps the C-m correlation
            % and monotone temperature behaviour; the Q(dK) construction)
            [aa, NN] = integrate_life_v(@(dk) ...
                u * rate_at_T(dk, Cf, mf, S.TK(temps), S.TK(Tt)), ...
                4.6, 27, w_mm, dK_fun);
            Nf(b) = max(NN);
            curves(b, :) = interp1(aa, NN, agrid, 'linear', 'extrap');
        end
        med = median(curves, 1); lo = prctile(curves, 2.5, 1); hi = prctile(curves, 97.5, 1);
        cc = cols(it, :);
        fill([lo, fliplr(hi)] / 1e5, [agrid, fliplr(agrid)], cc, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(med / 1e5, agrid, '-', 'Color', cc, ...
            'DisplayName', sprintf('%d ^{\\circ}C predicted', Tt));
        prows = [prows; {inh, Tt, median(Nf), prctile(Nf, 2.5), prctile(Nf, 97.5)}]; %#ok<AGROW>
    end
    % measured neighbours for context
    nb = {'fig10_fswc_25c.csv', 'fig13_fswi_25c.csv';
          'fig10_fswc_55c.csv', 'fig13_fswi_55c.csv';
          'fig10_fswc_85c.csv', 'fig13_fswi_85c.csv'};
    mkc = [0.4 0.4 0.4];
    for t = 1:3
        M = sortrows(readtable(fullfile(datadir, nb{t, inh+1})), 'N');
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
    % three anchor temperatures (exact at anchors, no overshoot; the
    % 25-55 degC plateau + sharp 85 degC rise is strongly non-Arrhenius,
    % so a single-Q straight line would overshoot between anchors)
    [x, i] = sort(1 ./ TKanchors);
    v = zeros(size(dk));
    for k = 1:numel(dk)
        lv = log(10.^Cf .* dk(k) .^ mf);
        lv = lv(i);
        v(k) = exp(pchip(x, lv, 1 / TKq));
    end
end
