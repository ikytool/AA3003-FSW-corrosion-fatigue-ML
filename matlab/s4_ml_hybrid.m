%S4_ML_HYBRID Machine-learning layer on the digitized dataset.
%   Data: pre-corroded conditions only (FSW.C from Figs 12+14 = two
%   independent digitization passes, FSW.I from Fig 14).
%   Target y = log10(da/dN in plotted units); features [log10 dK, T(K), inh].
%
%   Models: physics baseline (Arrhenius-Paris), multilinear regression,
%   boosted trees, Gaussian process, hybrid (physics + GPR residual).
%
%   Protocols:
%     LOCO    - leave-one-curve-out (9 curves; group CV avoids leakage
%               between the two digitization passes of the same test)
%     LOTO-55 - train at 25 and 85 degC only, predict all 55 degC data
%
%   Outputs: figures/fig_parity.pdf, figures/fig_loto.pdf,
%            data/derived/model_metrics.csv

clear; close all;
fig_defaults();
S = source_data();

here    = fileparts(mfilename('fullpath'));
datadir = fullfile(here, '..', 'data', 'digitized');
dervdir = fullfile(here, '..', 'data', 'derived');
if ~exist(dervdir, 'dir'), mkdir(dervdir); end

files = { ...
 'fig12_fswc_25c', 25, 0; 'fig12_fswc_55c', 55, 0; 'fig12_fswc_85c', 85, 0;
 'fig14_fswc_25c', 25, 0; 'fig14_fswc_55c', 55, 0; 'fig14_fswc_85c', 85, 0;
 'fig14_fswi_25c', 25, 1; 'fig14_fswi_55c', 55, 1; 'fig14_fswi_85c', 85, 1};

X = []; y = []; grp = [];
for k = 1:size(files, 1)
    T = readtable(fullfile(datadir, [files{k,1} '.csv']));
    n = height(T);
    X = [X; log10(T.dK), repmat(S.TK(files{k,2}), n, 1), repmat(files{k,3}, n, 1)]; %#ok<AGROW>
    y = [y; log10(T.dadN)]; %#ok<AGROW>
    % group by CONDITION (T, inh): the fig12/fig14 duplicates of the
    % same test must leave the training set together (no leakage)
    grp = [grp; repmat(files{k,2} * 10 + files{k,3}, n, 1)]; %#ok<AGROW>
end
n = numel(y);
fprintf('Dataset: %d points, %d conditions\n', n, numel(unique(grp)));

models = {'physics', 'mlr', 'gbt', 'gpr', 'hybrid'};
rng(1);

% ---------------- Protocol 1: leave-one-curve-out --------------------
metrics = [];
yhat_cv = nan(n, numel(models));
pi_cv   = nan(n, 2, numel(models));   % 95% intervals where available
for g = unique(grp)'
    te = grp == g; tr = ~te;
    for im = 1:numel(models)
        [yp, pint] = fit_predict(models{im}, X(tr,:), y(tr), X(te,:));
        yhat_cv(te, im) = yp;
        if ~isempty(pint), pi_cv(te, :, im) = pint; end
    end
end
for im = 1:numel(models)
    metrics = [metrics; pack(models{im}, 'LOCO', y, yhat_cv(:,im), pi_cv(:,:,im))]; %#ok<AGROW>
end

% ---------------- Protocol 2: leave-one-temperature-out --------------
T55 = abs(X(:,2) - S.TK(55)) < 1;
yhat55 = nan(nnz(T55), numel(models));
pi55   = nan(nnz(T55), 2, numel(models));
for im = 1:numel(models)
    [yp, pint] = fit_predict(models{im}, X(~T55,:), y(~T55), X(T55,:));
    yhat55(:, im) = yp;
    if ~isempty(pint), pi55(:, :, im) = pint; end
    metrics = [metrics; pack(models{im}, 'LOTO-55C', y(T55), yp, pi55(:,:,im))]; %#ok<AGROW>
end

% ---------------- Conformal calibration of the hybrid intervals ------
% Inflate the hybrid's intervals by the 95th percentile of the
% standardized cross-validation residuals (split-conformal idea): the
% intervals then cover 95% in interpolation by construction, and the
% honest report is the blind LOTO coverage with the SAME factor.
ih = find(strcmp(models, 'hybrid'));
sd_cv0 = (pi_cv(:,2,ih) - pi_cv(:,1,ih)) / (2 * 1.96);
kcal = prctile(abs(y - yhat_cv(:,ih)) ./ sd_cv0, 95) / 1.96;
fprintf('Conformal inflation factor (from LOCO residuals): k = %.2f\n', kcal);

sd_550 = (pi55(:,2,ih) - pi55(:,1,ih)) / (2 * 1.96);
cov55_cal = 100 * mean(abs(y(T55) - yhat55(:,ih)) <= 1.96 * kcal * sd_550);
fprintf('Hybrid LOTO-55 coverage: %.1f%% raw -> %.1f%% conformal\n', ...
    100 * mean(abs(y(T55) - yhat55(:,ih)) <= 1.96 * sd_550), cov55_cal);
row = struct('model', 'hybrid_conformal', 'protocol', 'LOTO-55C', ...
    'RMSE_log', metrics(strcmp({metrics.model}, 'hybrid') & ...
    strcmp({metrics.protocol}, 'LOTO-55C')).RMSE_log, ...
    'R2', metrics(strcmp({metrics.model}, 'hybrid') & ...
    strcmp({metrics.protocol}, 'LOTO-55C')).R2, ...
    'MAPE_pct', NaN, 'PICP95_pct', cov55_cal);
metrics = [metrics; row];

Tm = struct2table(metrics);
writetable(Tm, fullfile(dervdir, 'model_metrics.csv'));
disp(Tm);

% ---------------- Permutation feature importance ---------------------
% Within each LOCO fold: replace one feature of the held-out points by
% random draws from the training marginal, measure the RMSE increase.
% Done for the pure-ML model (gbt) and the hybrid.
featnames = {'log10 dK', 'T', 'inhibitor'};
imp = zeros(2, 3);
rng(2);
for g = unique(grp)'
    te = grp == g; tr = ~te;
    for imn = 1:2
        nm = {'gbt', 'hybrid'};
        yp0 = fit_predict(nm{imn}, X(tr,:), y(tr), X(te,:));
        e0 = mean((y(te) - yp0).^2);
        for j = 1:3
            Xp = X(te,:);
            src = X(tr, j);
            Xp(:, j) = src(randi(numel(src), nnz(te), 1));
            ypj = fit_predict(nm{imn}, X(tr,:), y(tr), Xp);
            imp(imn, j) = imp(imn, j) + mean((y(te) - ypj).^2) - e0;
        end
    end
end
imp = sqrt(max(imp / numel(unique(grp)), 0));   % RMSE-increase scale
Ti = array2table(imp, 'VariableNames', featnames, 'RowNames', {'gbt','hybrid'});
writetable(Ti, fullfile(dervdir, 'feature_importance.csv'), 'WriteRowNames', true);
fprintf('Permutation importance (RMSE increase, decades):\n'); disp(Ti);

% ---------------- Statistical diagnostics figure ---------------------
% (a) QQ plot of standardized LOCO residuals of the hybrid model
% (b) prediction-interval calibration curve, LOCO and LOTO, with the
%     conformally calibrated LOTO curve
sd_cv = sd_cv0 * kcal;
z_cv  = (y - yhat_cv(:,ih)) ./ sd_cv;
sd_55 = sd_550 * kcal;
z_55  = (y(T55) - yhat55(:,ih)) ./ sd_55;

fig = figure('Visible', 'off');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
zs = sort(z_cv);
nq = numel(zs);
th = sqrt(2) * erfinv(2 * ((1:nq)' - 0.5) / nq - 1);
plot(th, zs, 'o', 'MarkerSize', 4);
lim = [min([th; zs]) - 0.3, max([th; zs]) + 0.3];
plot(lim, lim, 'k--'); axis([lim lim]);
xlabel('normal quantiles');
ylabel('standardized residuals (hybrid, LOCO)');
title('(a)', 'FontWeight', 'normal');

nexttile; hold on;
qs = 0.50:0.05:0.95;
cov_cv = zeros(size(qs)); cov_55 = zeros(size(qs));
for iq = 1:numel(qs)
    zq = sqrt(2) * erfinv(qs(iq));
    cov_cv(iq) = mean(abs(z_cv) <= zq);
    cov_55(iq) = mean(abs(z_55) <= zq);
end
plot([0.45 1], [0.45 1], 'k--', 'DisplayName', 'ideal');
plot(qs, cov_cv, '-o', 'MarkerFaceColor', 'auto', 'MarkerSize', 4, ...
    'DisplayName', 'LOCO, calibrated');
plot(qs, cov_55, '-s', 'MarkerFaceColor', 'auto', 'MarkerSize', 4, ...
    'DisplayName', 'LOTO-55 (blind), calibrated');
xlabel('nominal coverage');
ylabel('empirical coverage');
xlim([0.45 1]); ylim([0.3 1.02]);
legend('Location', 'northwest');
title('(b)', 'FontWeight', 'normal');

export_vector(fig, 'fig_diagnostics', 180, 76);

% ---------------- Figure: parity plots (LOCO) ------------------------
fig = figure('Visible', 'off');
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
show = {'physics', 'mlr', 'gbt', 'hybrid'};
ttl  = {'(a) physics (Arrhenius-Paris)', '(b) multilinear regression', ...
        '(c) boosted trees', '(d) hybrid physics + GP'};
for i = 1:4
    im = find(strcmp(models, show{i}));
    nexttile; hold on;
    scatter(y, yhat_cv(:,im), 8, 'filled', 'MarkerFaceAlpha', 0.4);
    lims = [floor(min(y)*2)/2, ceil(max(y)*2)/2];
    plot(lims, lims, 'k--'); axis([lims lims]);
    r = metrics(strcmp({metrics.model}, show{i}) & strcmp({metrics.protocol}, 'LOCO'));
    text(lims(1)+0.1, lims(2)-0.25, sprintf('RMSE = %.3f\nR^2 = %.3f', ...
        r.RMSE_log, r.R2), 'FontSize', 8, 'VerticalAlignment', 'top');
    xlabel('measured log_{10}(da/dN)');
    ylabel('predicted log_{10}(da/dN)');
    title(ttl{i}, 'FontWeight', 'normal');
end
export_vector(fig, 'fig_parity', 180, 168);

% ---------------- Figure: LOTO blind prediction ----------------------
fig = figure('Visible', 'off');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
for panel = 1:2
    inh = panel - 1;
    nexttile; hold on;
    sel = T55 & X(:,3) == inh;
    hd = scatter(10.^X(sel,1), 10.^y(sel), 14, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.45, 'DisplayName', 'measured 55 ^{\circ}C (held out)');
    dkq = logspace(log10(3.4), log10(15.5), 60)';
    Xq = [log10(dkq), repmat(S.TK(55), 60, 1), repmat(inh, 60, 1)];
    sty = {'-', '--', ':', '-.'};
    cn = 0;
    for nm = {'physics', 'gbt', 'hybrid'}
        cn = cn + 1;
        [yq, pint] = fit_predict(nm{1}, X(~T55,:), y(~T55), Xq);
        plot(dkq, 10.^yq, sty{cn}, 'DisplayName', nm{1});
        if strcmp(nm{1}, 'hybrid') && ~isempty(pint)
            fill([dkq; flipud(dkq)], 10.^[pint(:,1); flipud(pint(:,2))], ...
                [0 114 178]/255, 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
                'DisplayName', 'hybrid 95% PI');
        end
    end
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('\DeltaK (MPa\cdotm^{1/2})');
    ylabel('da/dN (plotted units/cycle)');
    if inh == 0, title('(a) FSW.C, blind 55 ^{\circ}C', 'FontWeight', 'normal');
    else,        title('(b) FSW.I, blind 55 ^{\circ}C', 'FontWeight', 'normal'); end
    legend('Location', 'northwest');
end
export_vector(fig, 'fig_loto', 180, 81);

% ======================================================================
% All models are implemented in base MATLAB (no toolboxes):
% OLS, gradient-boosted depth-2 trees, GP with ARD-RBF kernel whose
% hyperparameters maximise the log marginal likelihood (fminsearch).
function [yp, pint] = fit_predict(model, Xtr, ytr, Xte)
    pint = [];
    switch model
        case 'physics'
            yp = physics_baseline(Xtr, ytr, Xte);
        case 'mlr'
            A = [ones(size(Xtr,1),1) Xtr];
            beta = A \ ytr;
            yp = [ones(size(Xte,1),1) Xte] * beta;
        case 'gbt'
            mdl = boost_fit(Xtr, ytr, 400, 0.08);
            yp = boost_predict(mdl, Xte);
        case 'gpr'
            [yp, ps] = gp_fit_predict(Xtr, ytr, Xte);
            pint = [yp - 1.96*ps, yp + 1.96*ps];
        case 'hybrid'
            yb_tr = physics_baseline(Xtr, ytr, Xtr);
            yb_te = physics_baseline(Xtr, ytr, Xte);
            [dr, ps] = gp_fit_predict(Xtr, ytr - yb_tr, Xte);
            yp = yb_te + dr;
            % propagate the physics-baseline parameter uncertainty into
            % the interval: stratified bootstrap of the training data
            % (per condition), refit the baseline, add its variance to
            % the GP predictive variance
            B = 80;
            Yb = zeros(size(Xte, 1), B);
            key = Xtr(:,2) * 10 + Xtr(:,3);
            ukey = unique(key);
            for b = 1:B
                idx = [];
                for g = ukey'
                    rows = find(key == g);
                    idx = [idx; rows(randi(numel(rows), numel(rows), 1))]; %#ok<AGROW>
                end
                Yb(:, b) = physics_baseline(Xtr(idx,:), ytr(idx), Xte);
            end
            tot = sqrt(ps.^2 + var(Yb, 0, 2));
            pint = [yp - 1.96*tot, yp + 1.96*tot];
    end
end

% ---- gradient boosting with depth-2 regression trees ----------------
function mdl = boost_fit(X, y, rounds, lr)
    mdl.f0 = mean(y); mdl.lr = lr; mdl.trees = cell(rounds, 1);
    F = repmat(mdl.f0, size(y));
    for t = 1:rounds
        tr = tree2_fit(X, y - F);
        F = F + lr * tree2_predict(tr, X);
        mdl.trees{t} = tr;
    end
end
function yp = boost_predict(mdl, X)
    yp = repmat(mdl.f0, size(X,1), 1);
    for t = 1:numel(mdl.trees)
        yp = yp + mdl.lr * tree2_predict(mdl.trees{t}, X);
    end
end
function tr = tree2_fit(X, r)
    [tr.j, tr.s] = best_split(X, r);
    L = X(:, tr.j) <= tr.s;
    [tr.jl, tr.sl] = best_split(X(L,:),  r(L));
    [tr.jr, tr.sr] = best_split(X(~L,:), r(~L));
    LL = L  & X(:, tr.jl) <= tr.sl;  LR = L  & X(:, tr.jl) >  tr.sl;
    RL = ~L & X(:, tr.jr) <= tr.sr;  RR = ~L & X(:, tr.jr) >  tr.sr;
    tr.mu = [safemean(r(LL)) safemean(r(LR)) safemean(r(RL)) safemean(r(RR))];
end
function yp = tree2_predict(tr, X)
    L = X(:, tr.j) <= tr.s;
    yp = zeros(size(X,1), 1);
    yp( L & X(:,tr.jl) <= tr.sl) = tr.mu(1);
    yp( L & X(:,tr.jl) >  tr.sl) = tr.mu(2);
    yp(~L & X(:,tr.jr) <= tr.sr) = tr.mu(3);
    yp(~L & X(:,tr.jr) >  tr.sr) = tr.mu(4);
end
function [jb, sb] = best_split(X, r)
    jb = 1; sb = median(X(:,1)); best = inf;
    if numel(r) < 4, return; end
    for j = 1:size(X, 2)
        qs = unique(prctile(X(:,j), 10:5:90));
        for q = qs(:)'
            L = X(:,j) <= q;
            if nnz(L) < 2 || nnz(~L) < 2, continue; end
            sse = sum((r(L) - mean(r(L))).^2) + sum((r(~L) - mean(r(~L))).^2);
            if sse < best, best = sse; jb = j; sb = q; end
        end
    end
end
function m = safemean(v)
    if isempty(v), m = 0; else, m = mean(v); end
end

% ---- Gaussian process, ARD-RBF + noise, base MATLAB -----------------
function [yp, ps] = gp_fit_predict(Xtr, ytr, Xte)
    mux = mean(Xtr); sx = std(Xtr); sx(sx == 0) = 1;
    Xs = (Xtr - mux) ./ sx;  Xq = (Xte - mux) ./ sx;
    muy = mean(ytr); sy = std(ytr); if sy == 0, sy = 1; end
    ys = (ytr - muy) / sy;
    n = numel(ys);

    hyp0 = log([0.2 1 1 1 1]);   % [noise sd, signal sd, lengthscales x3]
    opt  = optimset('Display', 'off', 'MaxFunEvals', 400, 'MaxIter', 400);
    hyp  = fminsearch(@(h) gp_nll(h, Xs, ys), hyp0, opt);

    sn = exp(hyp(1)); sf = exp(hyp(2)); ell = exp(hyp(3:5));
    K  = sf^2 * rbf(Xs, Xs, ell) + (sn^2 + 1e-8) * eye(n);
    L  = chol(K, 'lower');
    alpha = L' \ (L \ ys);
    Ks = sf^2 * rbf(Xq, Xs, ell);
    ypq = Ks * alpha;
    v   = L \ Ks';
    vq  = max(sf^2 + sn^2 - sum(v.^2, 1)', 1e-10);
    yp  = muy + sy * ypq;
    ps  = sy * sqrt(vq);
end
function f = gp_nll(hyp, Xs, ys)
    sn = exp(hyp(1)); sf = exp(hyp(2)); ell = exp(hyp(3:5));
    n = numel(ys);
    K = sf^2 * rbf(Xs, Xs, ell) + (sn^2 + 1e-8) * eye(n);
    [L, flag] = chol(K, 'lower');
    if flag, f = 1e10; return; end
    alpha = L' \ (L \ ys);
    f = 0.5 * (ys' * alpha) + sum(log(diag(L))) + 0.5 * n * log(2*pi);
end
function K = rbf(A, B, ell)
    D = zeros(size(A,1), size(B,1));
    for d = 1:size(A,2)
        D = D + ((A(:,d) - B(:,d)') / ell(d)).^2;
    end
    K = exp(-0.5 * D);
end

function yb = physics_baseline(Xtr, ytr, Xte)
    % Arrhenius-Paris: per (T, inh) Paris fits on the training data,
    % then ln C linear in 1/T and m linear in T within each family.
    yb = zeros(size(Xte, 1), 1);
    for inh = 0:1
        Ts = unique(Xtr(Xtr(:,3) == inh, 2));
        Cf = zeros(size(Ts)); mf = zeros(size(Ts));
        for it = 1:numel(Ts)
            sel = Xtr(:,2) == Ts(it) & Xtr(:,3) == inh;
            p = polyfit(Xtr(sel,1), ytr(sel), 1);
            mf(it) = p(1); Cf(it) = p(2);
        end
        if numel(Ts) >= 2
            pc = polyfit(1 ./ Ts, Cf, 1); pm = polyfit(Ts, mf, 1);
        else
            pc = [0 Cf(1)]; pm = [0 mf(1)];
        end
        q = Xte(:,3) == inh;
        yb(q) = polyval(pc, 1 ./ Xte(q,2)) + polyval(pm, Xte(q,2)) .* Xte(q,1);
    end
end

function row = pack(model, protocol, y, yhat, pint)
    row.model = model; row.protocol = protocol;
    row.RMSE_log = sqrt(mean((y - yhat).^2));
    row.R2 = 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);
    row.MAPE_pct = mean(abs(10.^yhat - 10.^y) ./ 10.^y) * 100;
    if ~isempty(pint) && ~all(isnan(pint(:)))
        row.PICP95_pct = 100 * mean(y >= pint(:,1) & y <= pint(:,2));
    else
        row.PICP95_pct = NaN;
    end
end
