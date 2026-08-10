%S1_PARIS_BASELINE Refit all digitized crack-growth curves + QA gate.
%   Refits every digitized da/dN-dK curve (Figs 12 and 14) in log-log
%   space with bootstrap confidence intervals, applies the QA gate
%   against the Paris parameters printed in the source paper, and writes
%   data/derived/paris_refit.csv used by all later stages.
%
%   QA findings (2026-08-04): the fig12/fig14 duplicated conditions
%   agree to 0.1-0.2% (independent extractions), so the digitization is
%   faithful. BM, FSW.C 55C and 85C recover the printed fits within
%   0.2%. The printed FSW and FSW.C 25C legend fits do not match their
%   own plotted points; the printed FSW.I legend equations appear
%   shifted by one row (FSW.I 25C refits exactly to the constants
%   printed for FSW.C 25C). Refit values are used downstream; printed
%   values are kept for traceability.
%
%   Produces figures/fig_paris_published.pdf (digitized points + fits).

clear; close all;
fig_defaults();
S = source_data();

here    = fileparts(mfilename('fullpath'));
datadir = fullfile(here, '..', 'data', 'digitized');
dervdir = fullfile(here, '..', 'data', 'derived');
if ~exist(dervdir, 'dir'), mkdir(dervdir); end

% condition table: file | label | T (degC, NaN = not pre-corroded) |
% inhibited | published C, m (NaN = none printed / not trusted) | gate
conds = { ...
 'fig12_bm',       'BM',          NaN, 0, 3e-11, 3.7741, true;
 'fig12_fsw',      'FSW',         NaN, 0, 9e-12, 4.246,  true;
 'fig12_fswc_25c', 'FSW.C 25C',   25,  0, 7e-12, 4.5757, true;
 'fig12_fswc_55c', 'FSW.C 55C',   55,  0, 2e-11, 4.0869, true;
 'fig12_fswc_85c', 'FSW.C 85C',   85,  0, 7e-12, 5.0928, true;
 'fig14_fswc_25c', 'FSW.C 25C(b)',25,  0, 7e-12, 4.5757, false;
 'fig14_fswc_55c', 'FSW.C 55C(b)',55,  0, 2e-11, 4.0869, false;
 'fig14_fswc_85c', 'FSW.C 85C(b)',85,  0, 7e-12, 5.0928, false;
 'fig14_fswi_25c', 'FSW.I 25C',   25,  1, 9e-11, 3.7706, false;  % legend suspect
 'fig14_fswi_55c', 'FSW.I 55C',   55,  1, 3e-12, 5.4579, false;  % legend suspect
 'fig14_fswi_85c', 'FSW.I 85C',   85,  1, 2e-11, 4.3437, false};

nboot = 2000;
rng(1);
rows = [];
for i = 1:size(conds, 1)
    T = readtable(fullfile(datadir, [conds{i,1} '.csv']));
    lx = log10(T.dK); ly = log10(T.dadN); n = numel(lx);
    p  = polyfit(lx, ly, 1);
    m_fit = p(1); C_fit = 10^p(2);
    r  = ly - polyval(p, lx);
    R2 = 1 - sum(r.^2) / sum((ly - mean(ly)).^2);
    sigma = std(r);                       % scatter in decades

    mb = zeros(nboot,1); cb = zeros(nboot,1);
    for b = 1:nboot
        k = randi(n, n, 1);
        pb = polyfit(lx(k), ly(k), 1);
        mb(b) = pb(1); cb(b) = pb(2);
    end
    mci = prctile(mb, [2.5 97.5]);
    Cci = 10.^prctile(cb, [2.5 97.5]);
    % bootstrap convergence: CI endpoints from the first half of the
    % draws must agree with the full set
    mci_h = prctile(mb(1:nboot/2), [2.5 97.5]);
    convchk(i) = max(abs(mci_h - mci) ./ abs(mci)); %#ok<SAGROW>

    dev = NaN; verdict = "info";
    if ~isnan(conds{i,6})
        dev = 100 * (m_fit - conds{i,6}) / conds{i,6};
        if conds{i,7}
            verdict = string(ternary(abs(dev) <= 2, "PASS", "FAIL"));
        end
    end
    fprintf('%-13s n=%2d  m=%.4f [%.3f %.3f]  C=%.2e [%.1e %.1e]  R2=%.4f  sd=%.3f dex  dev=%+.2f%% %s\n', ...
        conds{i,2}, n, m_fit, mci, C_fit, Cci, R2, sigma, dev, verdict);

    rows = [rows; {conds{i,1}, conds{i,2}, conds{i,3}, conds{i,4}, n, ...
        C_fit, Cci(1), Cci(2), m_fit, mci(1), mci(2), R2, sigma, ...
        conds{i,5}, conds{i,6}, dev, char(verdict)}]; %#ok<AGROW>
end
Tq = cell2table(rows, 'VariableNames', {'file','label','T_C','inhibited','n', ...
    'C','C_lo','C_hi','m','m_lo','m_hi','R2','sigma_dex', ...
    'C_published','m_published','m_dev_pct','verdict'});
writetable(Tq, fullfile(dervdir, 'paris_refit.csv'));
fprintf('Refit table written to data/derived/paris_refit.csv\n');
fprintf('Bootstrap convergence: max CI-endpoint shift half vs full draws = %.2f%%\n', ...
    100 * max(convchk));

% ---------------- Figure: digitized data + refit lines ---------------
% Panel (a): the five uninhibited conditions (Fig. 12 extractions);
% panel (b): the paired inhibited/uninhibited family (Fig. 14).
fig = figure('Visible', 'off');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
sel = 1:5;
for i = sel
    T = readtable(fullfile(datadir, [conds{i,1} '.csv']));
    hs = scatter(T.dK, T.dadN, 12, 'filled', 'MarkerFaceAlpha', 0.55, ...
        'DisplayName', conds{i,2});
    dk = logspace(log10(min(T.dK)), log10(max(T.dK)), 50);
    plot(dk, Tq.C(i) .* dk.^Tq.m(i), '-', 'Color', hs.CData, ...
        'HandleVisibility', 'off');
end
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('\DeltaK (MPa\cdotm^{1/2})'); ylabel('da/dN (plotted units/cycle)');
legend('Location', 'northwest'); title('(a)', 'FontWeight', 'normal');

nexttile; hold on;
sel = [6 9 7 10 8 11];
for i = sel
    T = readtable(fullfile(datadir, [conds{i,1} '.csv']));
    hs = scatter(T.dK, T.dadN, 12, 'filled', 'MarkerFaceAlpha', 0.55, ...
        'DisplayName', strrep(conds{i,2}, '(b)', ''));
    dk = logspace(log10(min(T.dK)), log10(max(T.dK)), 50);
    plot(dk, Tq.C(i) .* dk.^Tq.m(i), '-', 'Color', hs.CData, ...
        'HandleVisibility', 'off');
end
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('\DeltaK (MPa\cdotm^{1/2})'); ylabel('da/dN (plotted units/cycle)');
legend('Location', 'northwest'); title('(b)', 'FontWeight', 'normal');

export_vector(fig, 'fig_paris_published', 180, 81);

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
