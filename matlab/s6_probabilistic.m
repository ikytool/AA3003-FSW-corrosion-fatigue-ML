%S6_PROBABILISTIC Monte Carlo life distributions and design maps.
%   Bootstrap (C, m) draws per condition -> Arrhenius interpolation in
%   temperature -> life integration -> life-vs-temperature design map
%   with 95% bands (a), and probability-of-failure curves (b).
%   Requires s1 (refits) and s5 (geometry calibration) to have run.
%
%   Outputs: figures/fig_design_maps.pdf,
%            data/derived/design_map.csv (B10/B50 lives vs temperature)

clear; close all;
fig_defaults();
S = source_data();

here    = fileparts(mfilename('fullpath'));
datadir = fullfile(here, '..', 'data', 'digitized');
dervdir = fullfile(here, '..', 'data', 'derived');

G = readtable(fullfile(dervdir, 'geometry_calibration.csv'));
w_mm = G.w_mm; u = G.unit_factor;
B_mm = 1.8; dP_N = 0.9 * 3500;
fgeo   = @(aw) sqrt(aw) .* (1.99 - 0.41*aw + 18.7*aw.^2 - 38.48*aw.^3 + 53.85*aw.^4);
dK_fun = @(a, w) dP_N ./ (B_mm .* sqrt(w)) .* fgeo(a ./ w) ./ sqrt(1000);

temps = [25 55 85];
rng(1);
nboot = 400;
Tgrid = 25:5:85;
a0 = 4.6; af = 27;

fig = figure('Visible', 'off');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
cols = [0 114 178; 213 94 0] / 255;
rows = [];
NfT = cell(1, 2);
for f = 1:2
    inh = f - 1;
    D = cell(1, 3);
    for t = 1:3
        if inh == 0
            T1 = readtable(fullfile(datadir, sprintf('fig12_fswc_%dc.csv', temps(t))));
            T2 = readtable(fullfile(datadir, sprintf('fig14_fswc_%dc.csv', temps(t))));
            D{t} = [T1; T2];
        else
            D{t} = readtable(fullfile(datadir, sprintf('fig14_fswi_%dc.csv', temps(t))));
        end
    end
    Nf = zeros(nboot, numel(Tgrid));
    for b = 1:nboot
        Cf = zeros(1,3); mf = zeros(1,3);
        for t = 1:3
            d = D{t}; k = randi(height(d), height(d), 1);
            p = polyfit(log10(d.dK(k)), log10(d.dadN(k)), 1);
            mf(t) = p(1); Cf(t) = p(2);
        end
        for it = 1:numel(Tgrid)
            [~, NN] = integrate_life_v(@(dk) ...
                u * rate_at_T(dk, Cf, mf, S.TK(temps), S.TK(Tgrid(it))), ...
                a0, af, w_mm, dK_fun);
            Nf(b, it) = max(NN);
        end
    end
    NfT{f} = Nf;
    med = median(Nf); lo = prctile(Nf, 2.5); hi = prctile(Nf, 97.5);
    b10 = prctile(Nf, 10); b50 = med;
    rows = [rows; [repmat(inh, numel(Tgrid), 1), Tgrid(:), b10(:), b50(:), lo(:), hi(:)]]; %#ok<AGROW>

    nexttile(1); hold on;
    fill([Tgrid, fliplr(Tgrid)], [lo, fliplr(hi)] / 1e5, cols(f,:), ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    nm = 'FSW.C (no inhibitor)'; if inh, nm = 'FSW.I (10% EG)'; end
    plot(Tgrid, med / 1e5, '-o', 'Color', cols(f,:), 'MarkerFaceColor', cols(f,:), ...
        'MarkerSize', 4, 'DisplayName', nm);
end
nexttile(1);
xlabel('Pre-corrosion temperature (^{\circ}C)');
ylabel('Median life (10^5 cycles)');
legend('Location', 'southwest');
title('(a)', 'FontWeight', 'normal');

% ---- Panel (b): probability of failure vs N at 45 and 65 degC -------
nexttile(2); hold on;
sty = {'-', '--'};
for f = 1:2
    for it = 1:2
        Tt = [45 65];
        idx = find(Tgrid == Tt(it));
        Nfs = sort(NfT{f}(:, idx));
        pof = (1:nboot)' / (nboot + 1);
        nm = sprintf('%s, %d ^{\\circ}C', ternstr(f==1, 'FSW.C', 'FSW.I'), Tt(it));
        plot(Nfs / 1e5, pof * 100, sty{it}, 'Color', cols(f,:), 'DisplayName', nm);
    end
end
yline(10, ':', 'B10', 'HandleVisibility', 'off', 'FontSize', 8, ...
    'LabelHorizontalAlignment', 'left');
yline(50, ':', 'B50', 'HandleVisibility', 'off', 'FontSize', 8, ...
    'LabelHorizontalAlignment', 'left');
xlabel('N (10^5 cycles)');
ylabel('Probability of failure (%)');
legend('Location', 'southeast');
title('(b)', 'FontWeight', 'normal');

export_vector(fig, 'fig_design_maps', 180, 76);

Td = array2table(rows, 'VariableNames', {'inhibited','T_C','B10','B50','lo95','hi95'});
writetable(Td, fullfile(dervdir, 'design_map.csv'));
fprintf('B50 life at 65C: FSW.C %.2e, FSW.I %.2e (gain %.0f%%)\n', ...
    Td.B50(Td.inhibited==0 & Td.T_C==65), Td.B50(Td.inhibited==1 & Td.T_C==65), ...
    100*(Td.B50(Td.inhibited==1 & Td.T_C==65)/Td.B50(Td.inhibited==0 & Td.T_C==65)-1));

function [aa, NN] = integrate_life_v(vfun, a0, af, w, dK_fun)
    aa = linspace(a0, af, 400)';
    v  = max(vfun(dK_fun(aa, w)), 1e-12);
    NN = cumtrapz(aa, 1 ./ v);
end

function v = rate_at_T(dk, Cf, mf, TKanchors, TKq)
    % shape-preserving interpolation of ln v in 1/T (see s5)
    [x, i] = sort(1 ./ TKanchors);
    v = zeros(size(dk));
    for k = 1:numel(dk)
        lv = log(10.^Cf .* dk(k) .^ mf);
        lv = lv(i);
        v(k) = exp(pchip(x, lv, 1 / TKq));
    end
end
function s = ternstr(c, a, b)
    if c, s = a; else, s = b; end
end
