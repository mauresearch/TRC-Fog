function [ParetoSet, ParetoFront] = SPEA2(N, T, lb, ub, dim, fobj, numObjectives, opts)
% SPEA2.
if nargin < 8 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'integer'), opts.integer = false; end
pc = 0.9; etac = 20; pm = 1 / dim; etam = 20;
lb = expand(lb, dim); ub = expand(ub, dim);
M = numObjectives; Nbar = N;

P = repmat(lb, N, 1) + rand(N, dim) .* repmat(ub - lb, N, 1);
P = fixx(P, lb, ub, opts);
FP = eval_all(P, fobj, M);
A = zeros(0, dim); FA = zeros(0, M);

for t = 1:T + 1
    U = [P; A]; FU = [FP; FA];
    [U, FU] = unique_rows(U, FU);
    fit = spea2_fitness(FU);
    [A, FA, fitA] = env_selection(U, FU, fit, Nbar);
    if t == T + 1, break; end
    nA = size(A, 1);
    par = zeros(N, 1);
    for k = 1:N
        a = randi(nA); b = randi(nA);
        if fitA(a) < fitA(b), par(k) = a; else, par(k) = b; end
    end
    P = zeros(N, dim);
    for k = 1:2:N
        [c1, c2] = sbx(A(par(k), :), A(par(min(k + 1, N)), :), lb, ub, pc, etac);
        P(k, :) = pmut(c1, lb, ub, pm, etam);
        if k + 1 <= N, P(k + 1, :) = pmut(c2, lb, ub, pm, etam); end
    end
    P = fixx(P, lb, ub, opts);
    FP = eval_all(P, fobj, M);
end
if T < 1
    fit = spea2_fitness(FP); [A, FA] = env_selection(P, FP, fit, Nbar);
end
nd = nondominated(FA);
ParetoSet = A(nd, :); ParetoFront = FA(nd, :);
end

function fit = spea2_fitness(F)
n = size(F, 1);
le = true(n); lt = false(n);
for k = 1:size(F, 2)
    c = F(:, k); le = le & bsxfun(@le, c, c'); lt = lt | bsxfun(@lt, c, c');
end
D = le & lt;
S = sum(D, 2);
R = D' * S;
dist = pdist_sq(F);
dist(1:n+1:end) = inf;
k = max(1, floor(sqrt(n)));
if n == 1
    sigma = 0;
else
    sd = sort(dist, 2);
    sigma = sqrt(sd(:, min(k, n - 1)));
end
fit = R + 1 ./ (sigma + 2);
end

function [A, FA, fitA] = env_selection(U, FU, fit, Nbar)
nd = find(fit < 1);
if numel(nd) <= Nbar
    [~, o] = sort(fit);
    sel = o(1:min(Nbar, numel(o)));
else
    sel = nd;
    Fs = FU(sel, :);
    dist = sqrt(pdist_sq(Fs)); dist(1:numel(sel)+1:end) = inf;
    while numel(sel) > Nbar
        sd = sort(dist, 2);
        [~, order] = sortrows(sd);
        w = order(1);
        sel(w) = []; dist(w, :) = []; dist(:, w) = [];
    end
end
A = U(sel, :); FA = FU(sel, :); fitA = fit(sel);
end

function D = pdist_sq(F)
s = sum(F.^2, 2);
D = max(bsxfun(@plus, s, s') - 2 * (F * F'), 0);
end

function nd = nondominated(F)
n = size(F, 1); nd = true(n, 1);
for i = 1:n
    if any(all(bsxfun(@le, F, F(i, :)), 2) & any(bsxfun(@lt, F, F(i, :)), 2)), nd(i) = false; end
end
end

function [U, FU] = unique_rows(U, FU)
[~, ia] = unique(U, 'rows', 'stable'); U = U(ia, :); FU = FU(ia, :);
end
function v = expand(v, d), v = v(:)'; if numel(v) == 1, v = repmat(v, 1, d); end, end
function X = fixx(X, lb, ub, opts)
X = min(max(X, repmat(lb, size(X, 1), 1)), repmat(ub, size(X, 1), 1));
if opts.integer, X = round(X); end
end
function F = eval_all(X, fobj, M)
F = zeros(size(X, 1), M);
for i = 1:size(X, 1)
    try
        f = real(fobj(X(i, :))); f = f(:)';
        if numel(f) ~= M || any(~isfinite(f)), f = 1e6 * ones(1, M); end
    catch
        f = 1e6 * ones(1, M);
    end
    F(i, :) = f;
end
end
function [c1, c2] = sbx(p1, p2, lb, ub, pc, eta)
c1 = p1; c2 = p2;
if rand() > pc, return; end
for j = 1:numel(p1)
    if rand() > 0.5 || abs(p1(j) - p2(j)) < 1e-14, continue; end
    y1 = min(p1(j), p2(j)); y2 = max(p1(j), p2(j));
    yl = lb(j); yu = ub(j); r = rand();
    beta = 1 + 2 * (y1 - yl) / (y2 - y1); alpha = 2 - beta^-(eta + 1);
    if r <= 1 / alpha, bq = (r * alpha)^(1 / (eta + 1)); else, bq = (1 / (2 - r * alpha))^(1 / (eta + 1)); end
    a1 = 0.5 * ((y1 + y2) - bq * (y2 - y1));
    beta = 1 + 2 * (yu - y2) / (y2 - y1); alpha = 2 - beta^-(eta + 1);
    if r <= 1 / alpha, bq = (r * alpha)^(1 / (eta + 1)); else, bq = (1 / (2 - r * alpha))^(1 / (eta + 1)); end
    a2 = 0.5 * ((y1 + y2) + bq * (y2 - y1));
    a1 = min(max(a1, yl), yu); a2 = min(max(a2, yl), yu);
    if rand() <= 0.5, c1(j) = a2; c2(j) = a1; else, c1(j) = a1; c2(j) = a2; end
end
end
function y = pmut(x, lb, ub, pm, eta)
y = x;
for j = 1:numel(x)
    if rand() >= pm || ub(j) <= lb(j), continue; end
    d1 = (x(j) - lb(j)) / (ub(j) - lb(j)); d2 = (ub(j) - x(j)) / (ub(j) - lb(j));
    r = rand(); mp = 1 / (eta + 1);
    if r <= 0.5
        val = 2 * r + (1 - 2 * r) * (1 - d1)^(eta + 1); dq = val^mp - 1;
    else
        val = 2 * (1 - r) + 2 * (r - 0.5) * (1 - d2)^(eta + 1); dq = 1 - val^mp;
    end
    y(j) = min(max(x(j) + dq * (ub(j) - lb(j)), lb(j)), ub(j));
end
end
