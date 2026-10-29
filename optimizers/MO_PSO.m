function [ParetoSet, ParetoFront] = MO_PSO(N, T, lb, ub, dim, fobj, numObjectives, opts)
% Multi-objective particle swarm optimizer.
if nargin < 8 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'integer'), opts.integer = false; end
w = 0.5; wdamp = 0.99; c1 = 1; c2 = 2;
nGrid = 7; gAlpha = 0.1; gBeta = 2; gGamma = 2; mu = 0.1;
lb = expand(lb, dim); ub = expand(ub, dim);
M = numObjectives; nRep = N;

X = repmat(lb, N, 1) + rand(N, dim) .* repmat(ub - lb, N, 1);
X = fixx(X, lb, ub, opts);
V = zeros(N, dim);
F = eval_all(X, fobj, M);
PB = X; FPB = F;
[RX, RF] = rep_update(zeros(0, dim), zeros(0, M), X, F, nRep, nGrid, gAlpha, gGamma);

for t = 1:T
    G = grid_index(RF, nGrid, gAlpha);
    for i = 1:N
        L = RX(select_member(G, gBeta), :);
        V(i, :) = w * V(i, :) + c1 * rand(1, dim) .* (PB(i, :) - X(i, :)) ...
                              + c2 * rand(1, dim) .* (L - X(i, :));
        xn = X(i, :) + V(i, :);
        out = xn < lb | xn > ub;
        V(i, out) = -V(i, out);
        pmu = (1 - (t - 1) / max(T - 1, 1))^(1 / mu);
        if rand() < pmu
            j = randi(dim);
            xn(j) = xn(j) + pmu * (ub(j) - lb(j)) * (2 * rand() - 1);
        end
        X(i, :) = fixx(xn, lb, ub, opts);
    end
    F = eval_all(X, fobj, M);
    for i = 1:N
        if dominates(F(i, :), FPB(i, :)) || (~dominates(FPB(i, :), F(i, :)) && rand() < 0.5)
            PB(i, :) = X(i, :); FPB(i, :) = F(i, :);
        end
    end
    [RX, RF] = rep_update(RX, RF, X, F, nRep, nGrid, gAlpha, gGamma);
    w = w * wdamp;
end
ParetoSet = RX; ParetoFront = RF;
end

function [RX, RF] = rep_update(RX, RF, X, F, nRep, nGrid, a, g)
RX = [RX; X]; RF = [RF; F];
[RF, ia] = unique(RF, 'rows', 'stable'); RX = RX(ia, :);
nd = nondominated(RF); RX = RX(nd, :); RF = RF(nd, :);
while size(RX, 1) > nRep
    G = grid_index(RF, nGrid, a);
    k = select_member(G, -g);
    RX(k, :) = []; RF(k, :) = [];
end
end
function G = grid_index(F, nGrid, a)
[n, M] = size(F); G = ones(n, 1);
mult = 1;
for k = 1:M
    lo = min(F(:, k)); hi = max(F(:, k)); d = hi - lo;
    lo = lo - a * d; hi = hi + a * d;
    if hi <= lo, sub = ones(n, 1);
    else, sub = min(nGrid, floor((F(:, k) - lo) / (hi - lo) * nGrid) + 1); end
    G = G + (sub - 1) * mult; mult = mult * nGrid;
end
end
function k = select_member(G, beta)
[cells, ~, ic] = unique(G);
cnt = accumarray(ic, 1);
p = cnt .^ (-beta); p = p / sum(p);
c = find(rand() <= cumsum(p), 1, 'first'); if isempty(c), c = numel(cells); end
mem = find(ic == c);
k = mem(randi(numel(mem)));
end
function nd = nondominated(F)
n = size(F, 1); nd = true(n, 1);
for i = 1:n
    if any(all(bsxfun(@le, F, F(i, :)), 2) & any(bsxfun(@lt, F, F(i, :)), 2)), nd(i) = false; end
end
end
function d = dominates(a, b), d = all(a <= b) && any(a < b); end
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
