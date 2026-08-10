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

% ---- Panel (b): eta(T, dK) from the digitized paired curves ---------
nexttile; hold on;
temps = [25 55 85];
rows = [];
for t = 1:numel(temps)
    Tc = readtable(fullfile(datadir, sprintf('fig14_fswc_%dc.csv', temps(t))));
    Ti = readtable(fullfile(datadir, sprintf('fig14_fswi_%dc.csv', temps(t))));
    lo = max(min(Tc.dK), min(Ti.dK));
    hi = min(max(Tc.dK), max(Ti.dK));
    dk = logspace(log10(lo), log10(hi), 40);
    % smooth interpolation of each curve in log-log space
    vc = 10.^interp1(log10(Tc.dK), smoothloess(log10(Tc.dK), log10(Tc.dadN)), log10(dk), 'pchip');
    vi = 10.^interp1(log10(Ti.dK), smoothloess(log10(Ti.dK), log10(Ti.dadN)), log10(dk), 'pchip');
    eta = 100 * (1 - vi ./ vc);
    plot(dk, eta, '-', 'DisplayName', sprintf('%d ^{\\circ}C (fatigue)', temps(t)));
    rows = [rows; [repmat(temps(t), numel(dk), 1), dk(:), eta(:)]]; %#ok<AGROW>
end
% gravimetric 10% EG reference points at 85 degC
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

writetable(array2table(rows, 'VariableNames', {'T_C','dK','eta_pct'}), ...
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
