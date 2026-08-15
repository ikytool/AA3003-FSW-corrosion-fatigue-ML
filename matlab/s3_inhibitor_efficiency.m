%S3_INHIBITOR_EFFICIENCY Inhibitor efficiency across the two studies.
%   Produces figures/fig_inhibitor_efficiency.pdf:
%   (a) Gravimetric inhibition efficiency vs ethylene glycol
%       concentration at 85 degC (BM and FSW), Colloids paper Table 4.
%   (b) Fatigue inhibitor efficiency eta(T, dK) computed point-wise from
%       the digitized paired FSW.C / FSW.I curves of Fig. 14, with the
%       gravimetric 10% EG points shown for comparison.
%   Writes data/derived/eta_surface.csv for the hybrid model.

clear; close all;
fig_defaults();
S = source_data();

here    = fileparts(mfilename('fullpath'));
datadir = fullfile(here, '..', 'data', 'digitized');
dervdir = fullfile(here, '..', 'data', 'derived');
if ~exist(dervdir, 'dir'), mkdir(dervdir); end

fig = figure('Visible', 'off');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% ---- Panel (a): gravimetric efficiency vs concentration -------------
nexttile; hold on;
c = S.grav.eg_conc_pct(2:end);
plot(c, S.grav.eta_BM_eg(2:end),  '-o', 'MarkerFaceColor', 'auto', ...
    'DisplayName', 'BM (gravimetry, 85 ^{\circ}C)');
plot(c, S.grav.eta_FSW_eg(2:end), '-s', 'MarkerFaceColor', 'auto', ...
    'DisplayName', 'FSW (gravimetry, 85 ^{\circ}C)');
xlabel('Ethylene glycol concentration (%)');
ylabel('Inhibition efficiency (%)');
xlim([5 45]); ylim([35 75]);
legend('Location', 'northwest');
title('(a)', 'FontWeight', 'normal');

% ---- Panel (b): eta(T, dK) from the paired deduplicated curves ------
% Joint bootstrap: both curves of a pair are resampled together and the
% ratio recomputed, giving simultaneous uncertainty bands on eta.
nexttile; hold on;
temps = [25 55 85];
rows = [];
rng(3);
nboot = 300;
cols = get(gca, 'ColorOrder');
eta85hi = [];   % high-dK 85C band, for the gravimetric comparison
for t = 1:numel(temps)
    Tc = readtable(fullfile(datadir, sprintf('fig12_fswc_%dc.csv', temps(t))));
    Ti = readtable(fullfile(datadir, sprintf('fig14_fswi_%dc.csv', temps(t))));
    lo = max(min(Tc.dK), min(Ti.dK));
    hi = min(max(Tc.dK), max(Ti.dK));
    dk = logspace(log10(lo), log10(hi), 40);
    eta = eta_of(Tc, Ti, dk);
    Eb = zeros(nboot, numel(dk));
    for b = 1:nboot
        kc = randi(height(Tc), height(Tc), 1);
        ki = randi(height(Ti), height(Ti), 1);
        Eb(b, :) = eta_of(Tc(kc,:), Ti(ki,:), dk);
    end
    elo = prctile(Eb, 2.5); ehi = prctile(Eb, 97.5);
    fill([dk, fliplr(dk)], [elo, fliplr(ehi)], cols(t,:), ...
        'FaceAlpha', 0.13, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(dk, eta, '-', 'Color', cols(t,:), ...
        'DisplayName', sprintf('%d ^{\\circ}C (fatigue)', temps(t)));
    rows = [rows; [repmat(temps(t), numel(dk), 1), dk(:), eta(:), elo(:), ehi(:)]]; %#ok<AGROW>
    if temps(t) == 85
        iq = dk >= 9;
        eta85hi = [interp1(dk, elo, 10), interp1(dk, ehi, 10)];
    end
end
% gravimetric 10% EG reference point at 85 degC
plot(15.5, S.grav.eta_FSW_eg(2), 'p', 'MarkerSize', 11, ...
    'MarkerFaceColor', [0 0 0], 'Color', [0 0 0], ...
    'DisplayName', 'gravimetry, 10% EG');
yline(0, ':', 'HandleVisibility', 'off');
set(gca, 'XScale', 'log');
xlabel('\DeltaK (MPa\cdotm^{1/2})');
ylabel('Inhibitor efficiency \eta (%)');
ylim([-25 95]);
legend('Location', 'northwest');
title('(b)', 'FontWeight', 'normal');
fprintf('85C bootstrap band at dK = 10: [%.1f, %.1f]%%; gravimetric value %.1f%%\n', ...
    eta85hi(1), eta85hi(2), S.grav.eta_FSW_eg(2));

writetable(array2table(rows, 'VariableNames', {'T_C','dK','eta_pct','eta_lo','eta_hi'}), ...
    fullfile(dervdir, 'eta_surface.csv'));
export_vector(fig, 'fig_inhibitor_efficiency', 180, 76);

fprintf('eta at dK=5:  ');
for t = temps
    r = rows(rows(:,1)==t, :);
    fprintf('%dC: %.1f%%  ', t, interp1(r(:,2), r(:,3), 5));
end
fprintf('\neta at dK=10: ');
for t = temps
    r = rows(rows(:,1)==t, :);
    fprintf('%dC: %.1f%%  ', t, interp1(r(:,2), r(:,3), 10));
end
fprintf('\n');

function eta = eta_of(Tc, Ti, dk)
    % efficiency from a pair of curves, smoothed in log-log space
    [xc, ic] = unique(log10(Tc.dK)); yc = log10(Tc.dadN); yc = yc(ic);
    [xi, ii] = unique(log10(Ti.dK)); yi = log10(Ti.dadN); yi = yi(ii);
    vc = 10.^interp1(xc, smoothloess(xc, yc), log10(dk), 'pchip');
    vi = 10.^interp1(xi, smoothloess(xi, yi), log10(dk), 'pchip');
    eta = 100 * (1 - vi ./ vc);
end

function ys = smoothloess(x, y)
    % simple local linear smoother (window of 7 neighbours)
    n = numel(x); ys = zeros(size(y));
    for i = 1:n
        [~, idx] = sort(abs(x - x(i)));
        k = idx(1:min(7, n));
        p = polyfit(x(k), y(k), 1);
        ys(i) = polyval(p, x(i));
    end
end
