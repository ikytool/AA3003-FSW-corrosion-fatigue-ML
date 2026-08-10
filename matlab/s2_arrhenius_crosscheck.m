%S2_ARRHENIUS_CROSSCHECK The two-study fusion figure.
%   Produces figures/fig_arrhenius_crosscheck.pdf:
%   (a) Arrhenius plots (normalised to 25 degC) of the crack growth rate
%       v* = C(T)*dK*^m(T) at fixed dK*, next to the gravimetric
%       corrosion rate in the same medium (solution Y).
%   (b) Apparent activation energy Q(dK) of crack growth versus the
%       corrosion activation energy Ea re-derived from the tabulated
%       corrosion rates, for both the uninhibited (FSW.C) and inhibited
%       (FSW.I) pre-corroded joints.
%   Uses the REFITTED Paris parameters from the digitized data
%   (data/derived/paris_refit.csv, produced by s1) when available, and
%   falls back to the published constants otherwise.

clear; close all;
fig_defaults();
S = source_data();
R = S.Rgas;

here = fileparts(mfilename('fullpath'));
refitfile = fullfile(here, '..', 'data', 'derived', 'paris_refit.csv');
TC = [25 55 85];
if exist(refitfile, 'file')
    Tq = readtable(refitfile);
    pick = @(lbl) Tq(strcmp(Tq.label, lbl), :);
    Cc = [pick('FSW.C 25C').C,  pick('FSW.C 55C').C,  pick('FSW.C 85C').C];
    mc = [pick('FSW.C 25C').m,  pick('FSW.C 55C').m,  pick('FSW.C 85C').m];
    Ci = [pick('FSW.I 25C').C,  pick('FSW.I 55C').C,  pick('FSW.I 85C').C];
    mi = [pick('FSW.I 25C').m,  pick('FSW.I 55C').m,  pick('FSW.I 85C').m];
    src = 'refitted digitized data';
else
    iT = ~isnan(S.paris.T_C);
    Cc = S.paris.C(iT); mc = S.paris.m(iT); Ci = []; mi = [];
    src = 'published fit constants';
end
fprintf('Using %s\n', src);

TK   = S.TK(TC);
invT = 1000 ./ TK;

fig = figure('Visible', 'off');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% ---- Panel (a): normalised Arrhenius lines --------------------------
nexttile; hold on;
dKstar = [5 8 15];
mk = {'o', 's', '^'};
for j = 1:numel(dKstar)
    v = Cc .* dKstar(j) .^ mc;
    plot(invT, log(v ./ v(1)), ['-' mk{j}], 'MarkerFaceColor', 'auto', ...
        'DisplayName', sprintf('\\DeltaK = %g', dKstar(j)));
end
Wg = S.grav.W_BM_Y;
plot(1000 ./ S.TK(S.grav.T_C), log(Wg ./ Wg(1)), '--d', 'Color', [0 0 0], ...
    'MarkerFaceColor', 'k', 'DisplayName', 'corrosion (gravim.)');
xlabel('1000/T (K^{-1})');
ylabel('ln(rate / rate at 25 ^{\circ}C)');
legend('Location', 'northeast');
title('(a)', 'FontWeight', 'normal');

% ---- Panel (b): Q(dK) for FSW.C and FSW.I vs corrosion Ea -----------
dKgrid = linspace(4, 16, 60);
Qc = arrayfun(@(dk) -polyfitslope(1./TK, log(Cc .* dk.^mc)) * R / 1000, dKgrid);
pg = polyfit(1 ./ S.TK(S.grav.T_C), log(Wg), 1);
Ea = -pg(1) * R / 1000;
fprintf('Corrosion Ea (solution Y, BM, re-derived): %.1f kJ/mol\n', Ea);
fprintf('Q(dK) FSW.C: %.1f to %.1f kJ/mol over dK = 4-16\n', min(Qc), max(Qc));

nexttile; hold on;
plot(dKgrid, Qc, '-', 'DisplayName', 'Q(\DeltaK), FSW.C');
if ~isempty(Ci)
    Qi = arrayfun(@(dk) -polyfitslope(1./TK, log(Ci .* dk.^mi)) * R / 1000, dKgrid);
    plot(dKgrid, Qi, '-.', 'DisplayName', 'Q(\DeltaK), FSW.I (10% EG)');
    fprintf('Q(dK) FSW.I: %.1f to %.1f kJ/mol over dK = 4-16\n', min(Qi), max(Qi));
end
yline(Ea, '--k', 'DisplayName', ...
    sprintf('E_a corrosion (re-derived) = %.1f kJ/mol', Ea));
xlabel('\DeltaK (MPa\cdotm^{1/2})');
ylabel('Apparent activation energy (kJ/mol)');
ylim([0 42]);
legend('Location', 'west');
title('(b)', 'FontWeight', 'normal');

export_vector(fig, 'fig_arrhenius_crosscheck', 180, 76);

function s = polyfitslope(x, y)
    p = polyfit(x, y, 1);
    s = p(1);
end
