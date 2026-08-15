%S6_PROBABILISTIC Conditional life scenarios with model-form envelope.
%   Replaces the former design-map/B10-B50 outputs. Because only three
%   pre-corrosion temperatures exist, the temperature dependence of the
%   crack-growth rate is NOT identifiable from the data alone; the
%   scenario spread across defensible interpolation forms is therefore
%   reported together with the parameter (bootstrap) uncertainty.
%
%   Forms, per inhibitor family, all anchored on the per-temperature
%   Paris fits of the deduplicated curves:
%     pchip   - shape-preserving interpolation of ln v in 1/T (reference)
%     arrh    - single straight line in 1/T (least squares, 3 anchors)
%     linT    - straight line in T (least squares, 3 anchors)
%     wcorr   - ln v linear in ln Wcorr(T), the measured gravimetric
%               corrosion rate of the companion study (solution Y)
%
%   Outputs: figures/fig_design_maps.pdf (scenario map),
%            data/derived/design_map.csv (life vs T per form),
%            data/derived/scenario_45_65.csv (45/65 degC summary)

clear; close all;
fig_defaults();
S = source_data();

here    = fileparts(mfilename('fullpath'));
datadir = fullfile(here, '..', 'data', 'digitized');
dervdir = fullfile(here, '..', 'data', 'derived');

G = readtable(fullfile(dervdir, 'geometry_calibration.csv'));
w_mm = G.w_mm; u = G.unit_factor;
af   = min(G.a_valid_mm, 24);
a0   = 4.6;
B_mm = 1.8; dP_N = 0.9 * 3500;
fgeo   = @(aw) sqrt(aw) .* (1.99 - 0.41*aw + 18.7*aw.^2 - 38.48*aw.^3 + 53.85*aw.^4);
dK_fun = @(a, w) dP_N ./ (B_mm .* sqrt(w)) .* fgeo(a ./ w) ./ sqrt(1000);

% gravimetric corrosion rate of BM in solution Y vs temperature
% (recovered unit g.m^-2.h^-1; only ratios matter here)
lnW_of_T = @(TC) interp1(S.grav.T_C, log(S.grav.W_BM_Y), TC, 'pchip');

temps = [25 55 85];
forms = {'pchip', 'arrh', 'linT', 'wcorr'};
formlab = {'PCHIP in 1/T', 'Arrhenius (1/T)', 'linear in T', 'corrosion-rate covariate'};
Tgrid = 25:5:85;
rng(1);
nboot = 200;

fig = figure('Visible', 'off');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
rows = []; srows = [];
for f = 1:2
    inh = f - 1;
    if inh == 0, fam = 'fig12_fswc'; else, fam = 'fig14_fswi'; end
    D = cell(1, 3);
    Cf0 = zeros(1,3); mf0 = zeros(1,3);
    for t = 1:3
        D{t} = readtable(fullfile(datadir, sprintf('%s_%dc.csv', fam, temps(t))));
        p = polyfit(log10(D{t}.dK), log10(D{t}.dadN), 1);
        mf0(t) = p(1); Cf0(t) = p(2);
    end

    % point-estimate life vs T for each interpolation form
    NfF = zeros(numel(forms), numel(Tgrid));
    for im = 1:numel(forms)
        for it = 1:numel(Tgrid)
            NfF(im, it) = life_at_T(forms{im}, Tgrid(it), Cf0, mf0, temps, ...
                S, lnW_of_T, u, a0, af, w_mm, dK_fun);
        end
        rows = [rows; table(repmat(inh, numel(Tgrid), 1), Tgrid', ...
            repmat(string(forms{im}), numel(Tgrid), 1), NfF(im,:)', ...
            'VariableNames', {'inhibited','T_C','form','Nf'})]; %#ok<AGROW>
    end

    % parameter (bootstrap) percentiles for the reference form
    NfB = zeros(nboot, numel(Tgrid));
    for b = 1:nboot
        Cf = zeros(1,3); mf = zeros(1,3);
        for t = 1:3
            d = D{t}; k = randi(height(d), height(d), 1);
            p = polyfit(log10(d.dK(k)), log10(d.dadN(k)), 1);
            mf(t) = p(1); Cf(t) = p(2);
        end
        for it = 1:numel(Tgrid)
            NfB(b, it) = life_at_T('pchip', Tgrid(it), Cf, mf, temps, ...
                S, lnW_of_T, u, a0, af, w_mm, dK_fun);
        end
    end
    lo = prctile(NfB, 2.5); hi = prctile(NfB, 97.5); med = median(NfB);

    % 45/65 summary: parameter percentiles (pchip) + model-form envelope
    for Tt = [45 65]
        it = find(Tgrid == Tt);
        env = [min(NfF(:,it)), max(NfF(:,it))];
        srows = [srows; table(inh, Tt, med(it), lo(it), hi(it), ...
            env(1), env(2), 'VariableNames', ...
            {'inhibited','T_C','Nf_pchip_median','Nf_param_lo','Nf_param_hi', ...
             'Nf_form_min','Nf_form_max'})]; %#ok<AGROW>
    end

    nexttile; hold on;
    cols = get(gca, 'ColorOrder');
    fill([Tgrid, fliplr(Tgrid)], [lo, fliplr(hi)] / 1e5, cols(1,:), ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
        'DisplayName', 'parameter band (PCHIP)');
    sty = {'-o', '--s', ':^', '-.d'};
    for im = 1:numel(forms)
        plot(Tgrid, NfF(im,:) / 1e5, sty{im}, 'MarkerSize', 4, ...
            'MarkerFaceColor', 'auto', 'DisplayName', formlab{im});
    end
    % measured anchor lives, truncated at the same crack interval
    nb = {'fig10_fswc_25c.csv','fig10_fswc_55c.csv','fig10_fswc_85c.csv';
          'fig13_fswi_25c.csv','fig13_fswi_55c.csv','fig13_fswi_85c.csv'};
    Nme = zeros(1,3);
    for t = 1:3
        M = sortrows(readtable(fullfile(datadir, nb{inh+1, t})), 'N');
        [au, iu] = unique(M.a);
        Nme(t) = interp1(au, M.N(iu), min(af, max(M.a))) - ...
                 interp1(au, M.N(iu), max(a0, min(M.a)));
    end
    plot(temps, Nme / 1e5, 'kp', 'MarkerSize', 9, 'MarkerFaceColor', 'k', ...
        'DisplayName', 'measured anchors');
    xlabel('Pre-corrosion temperature (^{\circ}C)');
    ylabel('Integrated life (10^5 cycles)');
    if inh == 0, title('(a) FSW.C, no inhibitor', 'FontWeight', 'normal');
    else,        title('(b) FSW.I, 10% ethylene glycol', 'FontWeight', 'normal'); end
    legend('Location', 'southwest', 'FontSize', 7);
end
export_vector(fig, 'fig_design_maps', 180, 85);

writetable(rows, fullfile(dervdir, 'design_map.csv'));
writetable(srows, fullfile(dervdir, 'scenario_45_65.csv'));
disp(srows);
gain45 = srows.Nf_pchip_median(srows.inhibited==1 & srows.T_C==45) / ...
         srows.Nf_pchip_median(srows.inhibited==0 & srows.T_C==45) - 1;
gain65 = srows.Nf_pchip_median(srows.inhibited==1 & srows.T_C==65) / ...
         srows.Nf_pchip_median(srows.inhibited==0 & srows.T_C==65) - 1;
fprintf('Scenario inhibitor gain (PCHIP form): %+.0f%% at 45C, %+.0f%% at 65C\n', ...
    100*gain45, 100*gain65);

% ======================================================================
function Nf = life_at_T(form, TC, Cf, mf, temps, S, lnW_of_T, u, a0, af, w, dK_fun)
    TK  = S.TK(temps);
    TKq = S.TK(TC);
    aa = linspace(a0, af, 400)';
    dk = dK_fun(aa, w);
    lv = zeros(numel(aa), 1);
    for k = 1:numel(aa)
        lva = log(10.^Cf .* dk(k) .^ mf);      % ln v at the three anchors
        switch form
            case 'pchip'
                [x, i] = sort(1 ./ TK);
                lv(k) = pchip(x, lva(i), 1 / TKq);
            case 'arrh'
                p = polyfit(1 ./ TK, lva, 1);
                lv(k) = polyval(p, 1 / TKq);
            case 'linT'
                p = polyfit(temps, lva, 1);
                lv(k) = polyval(p, TC);
            case 'wcorr'
                x = lnW_of_T(temps);
                p = polyfit(x, lva, 1);
                lv(k) = polyval(p, lnW_of_T(TC));
        end
    end
    v = max(u * exp(lv), 1e-12);               % mm/cycle
    NN = cumtrapz(aa, 1 ./ v);
    Nf = max(NN);
end
