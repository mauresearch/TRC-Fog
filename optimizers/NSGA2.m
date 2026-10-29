function [ParetoSet, ParetoFront] = NSGA2(N, T, lb, ub, dim, fobj, numObjectives, opts)
% NSGA-II.
if nargin < 8 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'integer'), opts.integer = false; end
pc = 0.9; etac = 20; pm = 1 / dim; etam = 20;
lb = expand(lb, dim); ub = expand(ub, dim);
M = numObjectives;

X = repmat(lb, N, 1) + rand(N, dim) .* repmat(ub - lb, N, 1);
X = fixx(X, lb, ub, opts);
F = eval_all(X, fobj, M);
[rk, cd] = rank_cd(F);

for t = 1:T
    par = zeros(N, 1);
    for k = 1:N
        a = randi(N); b = randi(N);
        if rk(a) < rk(b) || (rk(a) == rk(b) && cd(a) > cd(b)), par(k) = a; else, par(k) = b; end
    end
    Y = zeros(N, dim);
    for k = 1:2:N
        p1 = X(par(k), :); p2 = X(par(min(k + 1, N)), :);
        [c1, c2] = sbx(p1, p2, lb, ub, pc, etac);
        Y(k, :) = pmut(c1, lb, ub, pm, etam);
        if k + 1 <= N, Y(k + 1, :) = pmut(c2, lb, ub, pm, etam); end
    end
    Y = fixx(Y, lb, ub, opts);
    FY = eval_all(Y, fobj, M);
    [X, F] = env_select([X; Y], [F; FY], N);
    [rk, cd] = rank_cd(F);
end
nd = rk == 1;
ParetoSet = X(nd, :); ParetoFront = F(nd, :);
[ParetoFront, ia] = unique(ParetoFront, 'rows', 'stable'); ParetoSet = ParetoSet(ia, :);
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
function [rank, fronts] = nd_sort(F)
n = size(F, 1); le = true(n); lt = false(n);
for k = 1:size(F, 2)
    c = F(:, k); le = le & bsxfun(@le, c, c'); lt = lt | bsxfun(@lt, c, c');
end
D = le & lt; cnt = sum(D, 1)'; rank = zeros(n, 1); fronts = {}; r = 0; rem = true(n, 1);
while any(rem)
    r = r + 1; cur = find(rem & cnt == 0); if isempty(cur), cur = find(rem); end
    rank(cur) = r; fronts{r} = cur; rem(cur) = false;
    cnt = cnt - sum(D(cur, :), 1)'; cnt(~rem) = inf;
end
end
function cd = crowding(F)
n = size(F, 1); cd = zeros(n, 1);
if n <= 2, cd(:) = inf; return; end
for k = 1:size(F, 2)
    [v, idx] = sort(F(:, k)); cd(idx([1 end])) = inf; rg = v(end) - v(1);
    if rg <= 0, continue; end
    cd(idx(2:end-1)) = cd(idx(2:end-1)) + (v(3:end) - v(1:end-2)) / rg;
end
end
function [rk, cd] = rank_cd(F)
[rk, fronts] = nd_sort(F); cd = zeros(size(F, 1), 1);
for r = 1:numel(fronts), cd(fronts{r}) = crowding(F(fronts{r}, :)); end
end
function [Xs, Fs] = env_select(X, F, N)
[~, fronts] = nd_sort(F); sel = [];
for r = 1:numel(fronts)
    fr = fronts{r};
    if numel(sel) + numel(fr) <= N, sel = [sel; fr];
    else
        [~, o] = sort(crowding(F(fr, :)), 'descend'); sel = [sel; fr(o(1:N - numel(sel)))]; break;
    end
end
Xs = X(sel, :); Fs = F(sel, :);
end
