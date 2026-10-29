function [ParetoSet, ParetoFront] = MO_GWO(N, T, lb, ub, dim, fobj, numObjectives, opts)
% Multi-objective grey wolf optimizer.
if nargin < 8 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'integer'), opts.integer = false; end
nGrid = 10; gAlpha = 0.1; gBeta = 4; gGamma = 2;
lb = expand(lb, dim); ub = expand(ub, dim);
M = numObjectives; nArc = N;

X = repmat(lb, N, 1) + rand(N, dim) .* repmat(ub - lb, N, 1);
X = fixx(X, lb, ub, opts);
F = eval_all(X, fobj, M);
[AX, AF] = arc_update(zeros(0, dim), zeros(0, M), X, F, nArc, nGrid, gAlpha, gGamma);

for t = 1:T
    a = 2 - t * (2 / T);
    G = grid_index(AF, nGrid, gAlpha);
    for i = 1:N
        [iA, iB, iD] = pick_three(G, gBeta);
        Xa = AX(iA, :); Xb = AX(iB, :); Xd = AX(iD, :);
        A1 = 2 * a * rand(1, dim) - a; C1 = 2 * rand(1, dim);
        A2 = 2 * a * rand(1, dim) - a; C2 = 2 * rand(1, dim);
        A3 = 2 * a * rand(1, dim) - a; C3 = 2 * rand(1, dim);
        X1 = Xa - A1 .* abs(C1 .* Xa - X(i, :));
        X2 = Xb - A2 .* abs(C2 .* Xb - X(i, :));
        X3 = Xd - A3 .* abs(C3 .* Xd - X(i, :));
        X(i, :) = fixx((X1 + X2 + X3) / 3, lb, ub, opts);
    end
    F = eval_all(X, fobj, M);
    [AX, AF] = arc_update(AX, AF, X, F, nArc, nGrid, gAlpha, gGamma);
end
ParetoSet = AX; ParetoFront = AF;
end

function [iA, iB, iD] = pick_three(G, beta)
iA = select_member(G, beta, []);
ex = find(G == G(iA));
if numel(unique(G)) > 1, iB = select_member(G, beta, G(iA)); else, iB = pick_other(ex, iA); end
if numel(unique(G)) > 2, iD = select_member(G, beta, [G(iA); G(iB)]); else, iD = pick_other(1:numel(G), [iA iB]); end
end
function k = pick_other(pool, notThese)
c = setdiff(pool(:)', notThese(:)');
if isempty(c), k = notThese(1); else, k = c(randi(numel(c))); end
end
function [AX, AF] = arc_update(AX, AF, X, F, nArc, nGrid, a, g)
AX = [AX; X]; AF = [AF; F];
[AF, ia] = unique(AF, 'rows', 'stable'); AX = AX(ia, :);
nd = nondominated(AF); AX = AX(nd, :); AF = AF(nd, :);
while size(AX, 1) > nArc
    G = grid_index(AF, nGrid, a);
    k = select_member(G, -g, []);
    AX(k, :) = []; AF(k, :) = [];
end
end
function G = grid_index(F, nGrid, a)
[n, M] = size(F); G = ones(n, 1); mult = 1;
for k = 1:M
    lo = min(F(:, k)); hi = max(F(:, k)); d = hi - lo;
    lo = lo - a * d; hi = hi + a * d;
    if hi <= lo, sub = ones(n, 1);
    else, sub = min(nGrid, floor((F(:, k) - lo) / (hi - lo) * nGrid) + 1); end
    G = G + (sub - 1) * mult; mult = mult * nGrid;
end
end
function k = select_member(G, beta, excludeCells)
[cells, ~, ic] = unique(G);
cnt = accumarray(ic, 1);
p = cnt .^ (-beta);
p(ismember(cells, excludeCells)) = 0;
if sum(p) <= 0, p = ones(size(p)); end
p = p / sum(p);
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
