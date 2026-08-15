%S4_ML_HYBRID Machine-learning layer on the consolidated dataset.
%   Data: the six pre-corroded conditions, ONE curve per physical test
%   (FSW.C from Fig. 12, FSW.I from Fig. 14; the Fig. 14 replots of the
%   FSW.C tests are excluded as duplicates of the same specimens).
%   Target y = log10(da/dN, plotted units); features [log10 dK, T(K), inh],
%   where T is the PRE-CORROSION CONDITIONING temperature (fatigue tests
%   ran in air), never a crack-tip temperature.
%
%   Models: physics baseline (condition Paris fits + smooth T trend),
%   multilinear regression, boosted trees (reference only), Gaussian
%   process, hybrid (physics + GP residual).
%
%   Protocols (all elements refit inside each training fold, nothing
%   from a held-out group enters baselines, scaling or calibration):
%     LOCO     - leave-one-condition-out (6 folds)
%     LOTO-T   - leave-one-temperature-out for T = 25, 55, 85 degC
%                (55 is an interior hold-out; 25/85 are endpoint
%                extrapolation folds)
%   Interval inflation is NESTED: for each outer fold the factor comes
%   from an inner leave-one-condition-out pass on the training set only.
%
%   Outputs: figures/fig_parity.pdf, fig_loto.pdf, fig_diagnostics.pdf,
%            data/derived/model_metrics.csv, fold_contrasts.csv,
%            feature_importance.csv

clear; close all;
fig_defaults();
S = source_data();

here    = fileparts(mfilename('fullpath'));
datadir = fullfile(here, '..', 'data', 'digitized');
dervdir = fullfile(here, '..', 'data', 'derived');
if ~exist(dervdir, 'dir'), mkdir(dervdir); end

% one curve per physical specimen (deduplicated)
files = { ...
 'fig12_fswc_25c', 25, 0; 'fig12_fswc_55c', 55, 0; 'fig12_fswc_85c', 85, 0;
 'fig14_fswi_25c', 25, 1; 'fig14_fswi_55c', 55, 1; 'fig14_fswi_85c', 85, 1};

X = []; y = []; grp = [];
for k = 1:size(files, 1)
    T = readtable(fullfile(datadir, [files{k,1} '.csv']));
    n = height(T);
    X = [X; log10(T.dK), repmat(S.TK(files{k,2}), n, 1), repmat(files{k,3}, n, 1)]; %#ok<AGROW>
    y = [y; log10(T.dadN)]; %#ok<AGROW>
    grp = [grp; repmat(files{k,2} * 10 + files{k,3}, n, 1)]; %#ok<AGROW>
end
n = numel(y);
fprintf('Dataset: %d rate points from %d unique curves (6 conditions)\n', ...
    n, size(files, 1));

models = {'physics', 'mlr', 'gbt', 'gpr', 'hybrid'};
rng(1);

% ---------------- Protocol 1: leave-one-condition-out ----------------
metrics = [];
yhat_cv = nan(n, numel(models));
pi_cv   = nan(n, 2, numel(models));
for g = unique(grp)'
    te = grp == g; tr = ~te;
    kc = nested_kcal(X(tr,:), y(tr), grp(tr));        % training-only factor
    for im = 1:numel(models)
        [yp, pint] = fit_predict(models{im}, X(tr,:), y(tr), X(te,:));
        yhat_cv(te, im) = yp;
        if ~isempty(pint)
            mid = mean(pint, 2); half = (pint(:,2) - pint(:,1)) / 2;
            pi_cv(te, :, im) = [mid - kc*half, mid + kc*half];
        end
    end
end
for im = 1:numel(models)
    metrics = [metrics; pack(models{im}, 'LOCO', y, yhat_cv(:,im), pi_cv(:,:,im))]; %#ok<AGROW>
end

% ---------------- Protocol 2: all leave-one-temperature-out folds ----
foldT = [25 55 85];
yhatT = cell(1, 3); piT = cell(1, 3); selT = cell(1, 3);
for f = 1:3
    sel = abs(X(:,2) - S.TK(foldT(f))) < 1;
    selT{f} = sel;
    tr = ~sel;
    kc = nested_kcal(X(tr,:), y(tr), grp(tr));
    yhatT{f} = nan(nnz(sel), numel(models));
    piT{f}   = nan(nnz(sel), 2, numel(models));
    for im = 1:numel(models)
        [yp, pint] = fit_predict(models{im}, X(tr,:), y(tr), X(sel,:));
        yhatT{f}(:, im) = yp;
        if ~isempty(pint)
            mid = mean(pint, 2); half = (pint(:,2) - pint(:,1)) / 2;
            piT{f}(:, :, im) = [mid - kc*half, mid + kc*half];
        end
        metrics = [metrics; pack(models{im}, sprintf('LOTO-%dC', foldT(f)), ...
            y(sel), yp, piT{f}(:,:,im))]; %#ok<AGROW>
    end
end

Tm = struct2table(metrics);
writetable(Tm, fullfile(dervdir, 'model_metrics.csv'));
disp(Tm);

% ---------------- Fixed-dK contrasts per LOTO fold -------------------
% At reference dK values: measured (per-curve Paris fit of the held-out
% curve) versus predicted log rate, error in decades and as a
% multiplicative factor, plus the inhibitor rate ratio.
dKref = [5 8 12 15];
ih = find(strcmp(models, 'hybrid'));
crows = [];
for f = 1:3
    sel = selT{f};
    for inh = 0:1
        cs = sel & X(:,3) == inh;
        pmeas = polyfit(X(cs,1), y(cs), 1);      % held-out curve's own fit
        rng_dk = [min(10.^X(cs,1)), max(10.^X(cs,1))];
        for dk = dKref
            if dk < rng_dk(1) || dk > rng_dk(2), continue; end
            % hybrid prediction at this dK (trained without this T)
            Xq = [log10(dk), S.TK(foldT(f)), inh];
            yq = fit_predict('hybrid', X(~sel,:), y(~sel), Xq);
            ym = polyval(pmeas, log10(dk));
            crows = [crows; {sprintf('LOTO-%dC', foldT(f)), inh, dk, ...
                ym, yq, yq - ym, 10^(yq - ym)}]; %#ok<AGROW>
        end
    end
end
Tc = cell2table(crows, 'VariableNames', ...
    {'fold','inhibited','dK','log_v_measured','log_v_predicted', ...
     'err_decades','factor'});
writetable(Tc, fullfile(dervdir, 'fold_contrasts.csv'));
fprintf('Fixed-dK contrasts written (%d rows)\n', height(Tc));

% ---------------- Permutation feature importance ---------------------
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
imp = sqrt(max(imp / numel(unique(grp)), 0));
Ti = array2table(imp, 'VariableNames', featnames, 'RowNames', {'gbt','hybrid'});
writetable(Ti, fullfile(dervdir, 'feature_importance.csv'), 'WriteRowNames', true);
fprintf('Permutation importance (RMSE increase, decades):\n'); disp(Ti);

% ---------------- Statistical diagnostics figure ---------------------
% (a) QQ plot of standardized LOCO residuals of the hybrid model
% (b) empirical coverage of the nested-calibrated hybrid bands, LOCO
%     and the three temperature hold-out folds
sd_cv = (pi_cv(:,2,ih) - pi_cv(:,1,ih)) / (2 * 1.96);
z_cv  = (y - yhat_cv(:,ih)) ./ sd_cv;

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
plot([0.45 1], [0.45 1], 'k--', 'DisplayName', 'ideal');
covq = zeros(size(qs));
for iq = 1:numel(qs)
    zq = sqrt(2) * erfinv(qs(iq));
    covq(iq) = mean(abs(z_cv) <= zq);
end
plot(qs, covq, '-o', 'MarkerFaceColor', 'auto', 'MarkerSize', 4, ...
    'DisplayName', 'LOCO');
mk = {'-s', '-^', '-d'};
for f = 1:3
    sdf = (piT{f}(:,2,ih) - piT{f}(:,1,ih)) / (2 * 1.96);
    zf  = (y(selT{f}) - yhatT{f}(:,ih)) ./ sdf;
    for iq = 1:numel(qs)
        zq = sqrt(2) * erfinv(qs(iq));
        covq(iq) = mean(abs(zf) <= zq);
    end
    plot(qs, covq, mk{f}, 'MarkerFaceColor', 'auto', 'MarkerSize', 4, ...
        'DisplayName', sprintf('hold-out %d ^{\\circ}C', foldT(f)));
end
xlabel('nominal coverage');
ylabel('empirical coverage');
xlim([0.45 1]); ylim([0 1.02]);
legend('Location', 'southeast');
title('(b)', 'FontWeight', 'normal');

export_vector(fig, 'fig_diagnostics', 180, 76);

% ---------------- Figure: parity plots (LOCO) ------------------------
fig = figure('Visible', 'off');
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
show = {'physics', 'mlr', 'gbt', 'hybrid'};
ttl  = {'(a) physics (condition Paris + T trend)', '(b) multilinear regression', ...
        '(c) boosted trees (reference)', '(d) hybrid physics + GP'};
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

% ---------------- Figure: the three temperature hold-out folds -------
fig = figure('Visible', 'off');
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
pl = 'abcdef'; ip = 0;
for inh = 0:1
    for f = 1:3
        ip = ip + 1;
        nexttile; hold on;
        sel = selT{f} & X(:,3) == inh;
        scatter(10.^X(sel,1), 10.^y(sel), 10, 'k', 'filled', ...
            'MarkerFaceAlpha', 0.4, 'DisplayName', 'measured (held out)');
        dkq = logspace(log10(3.4), log10(15.5), 60)';
        Xq = [log10(dkq), repmat(S.TK(foldT(f)), 60, 1), repmat(inh, 60, 1)];
        tr = ~selT{f};
        cn = 0; sty = {'-', '--'};
        for nm = {'physics', 'hybrid'}
            cn = cn + 1;
            [yq, pint] = fit_predict(nm{1}, X(tr,:), y(tr), Xq);
            plot(dkq, 10.^yq, sty{cn}, 'DisplayName', nm{1});
            if strcmp(nm{1}, 'hybrid') && ~isempty(pint)
                fill([dkq; flipud(dkq)], 10.^[pint(:,1); flipud(pint(:,2))], ...
                    [0 114 178]/255, 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
                    'DisplayName', 'hybrid band');
            end
        end
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlabel('\DeltaK (MPa\cdotm^{1/2})');
        if f == 1, ylabel('da/dN (plotted units/cycle)'); end
        nmst = 'FSW.C'; if inh, nmst = 'FSW.I'; end
        title(sprintf('(%s) %s, hold-out %d ^{\\circ}C', pl(ip), nmst, ...
            foldT(f)), 'FontWeight', 'normal');
        if ip == 1, legend('Location', 'northwest', 'FontSize', 7); end
    end
end
export_vector(fig, 'fig_loto', 180, 120);

% ======================================================================
% Nested calibration: inner leave-one-condition-out on the TRAINING set
% only; returns the 95th-percentile inflation of standardized hybrid
% residuals. No held-out information enters this factor.
function kc = nested_kcal(Xtr, ytr, grptr)
    zs = [];
    for g = unique(grptr)'
        ite = grptr == g; itr = ~ite;
        if nnz(itr) < 20, continue; end
        [yp, pint] = fit_predict('hybrid', Xtr(itr,:), ytr(itr), Xtr(ite,:));
        if isempty(pint), continue; end
        sd = (pint(:,2) - pint(:,1)) / (2 * 1.96);
        zs = [zs; abs(ytr(ite) - yp) ./ sd]; %#ok<AGROW>
    end
    if isempty(zs), kc = 1; else, kc = max(1, prctile(zs, 95) / 1.96); end
end

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
            % add the physics-baseline parameter variance (stratified
            % bootstrap of the training data, refit per draw)
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
    % Condition-wise Paris fits on the TRAINING data only, then
    % ln C linear in 1/T and m linear in T within each inhibitor family.
    % T is the pre-corrosion conditioning temperature.
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
    row.factor = 10 ^ row.RMSE_log;      % multiplicative error factor
    row.R2 = 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);
    row.MAPE_pct = mean(abs(10.^yhat - 10.^y) ./ 10.^y) * 100;
    if ~isempty(pint) && ~all(isnan(pint(:)))
        row.PICP95_pct = 100 * mean(y >= pint(:,1) & y <= pint(:,2));
    else
        row.PICP95_pct = NaN;
    end
end
