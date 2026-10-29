function [ParetoSet, ParetoFront, Archive, ArchiveFitness] = IMOWO(N, T, lb, ub, dim, fobj, numObjectives, opts)
% Improved multi-objective walrus optimizer.

if nargin < 8 || isempty(opts), opts = struct(); end
if isfield(opts, 'useLHSOBL')
    opts = set_default(opts, 'useLHS', opts.useLHSOBL);
    opts = set_default(opts, 'useOBL', opts.useLHSOBL);
end
opts = set_default(opts, 'useLHS', true);
opts = set_default(opts, 'useOBL', true);
opts = set_default(opts, 'archive', 'knn');
opts = set_default(opts, 'leader', 'cd');
opts = set_default(opts, 'useLevy', true);
opts = set_default(opts, 'operators', {});
opts = set_default(opts, 'integer', false);
opts = set_default(opts, 'archiveSize', N);
opts = set_default(opts, 'verbose', false);

lb = lb(:)'; ub = ub(:)';
if numel(lb) == 1, lb = repmat(lb, 1, dim); end
if numel(ub) == 1, ub = repmat(ub, 1, dim); end
M = numObjectives;

if opts.useLHS
    U = zeros(N, dim);
    for j = 1:dim
        U(:, j) = (randperm(N)' - rand(N, 1)) / N;
    end
else
    U = rand(N, dim);
end
X1 = repmat(lb, N, 1) + U .* repmat(ub - lb, N, 1);
if opts.useOBL
    X2 = repmat(lb + ub, N, 1) - X1;
    Xc = fix_candidates([X1; X2], lb, ub, opts);
    Fc = eval_all(Xc, fobj, M);
    [X, F] = env_select(Xc, Fc, N);
else
    X = fix_candidates(X1, lb, ub, opts);
    F = eval_all(X, fobj, M);
end
[Archive, ArchiveFitness, ArchiveCD] = update_archive([], [], X, F, opts.archiveSize, opts.archive);

primesList = first_primes(dim);
pp = 0.4;
Fnum = round(N * pp);  Mnum = Fnum;  Cnum = N - Fnum - Mnum;

Tl = T - double(opts.useOBL);
for t = 1:Tl
    progress = t / Tl;
    Alpha = 1 - t / Tl;
    Beta  = 1 - 1 / (1 + exp((0.5*Tl - t) / Tl * 10));
    A = 2 * Alpha;
    [rk, ~] = nd_sort(F);

    Xadv = zeros(N, dim);  Fadv = zeros(N, M);
    for i = 1:N
        [Best, Second] = pick_leaders(Archive, ArchiveCD, X, rk, opts.leader);
        R = 2 * rand() - 1;
        Danger = A * R;
        Safety = rand();
        xi = X(i, :);

        if abs(Danger) >= 1
            a = randi(N); b = randi(N);
            xn = xi + (Beta * rand()^2) * (X(a, :) - X(b, :));

        elseif Safety >= 0.5
            if i <= Mnum
                h = zeros(1, dim);
                idx = (t - 1) * N + i;
                for d = 1:dim, h(d) = halton_val(idx, primesList(d)); end
                xn = lb + h .* (ub - lb);
            elseif i <= Mnum + Fnum
                male = X(i - Mnum, :);
                xn = xi + Alpha * (male - xi) + (1 - Alpha) * (Best - xi);
            else
                if opts.useLevy
                    stp = levy_step(dim);
                else
                    stp = 0.1 * randn(1, dim);
                end
                O = Best + xi .* stp;
                xn = rand() * (O - xi);
            end

        elseif abs(Danger) >= 0.5
            xn = xi * R - abs(Best - xi) * rand()^2;

        else
            xn = zeros(1, dim);
            for d = 1:dim
                a1 = Beta * rand() - Beta;  b1 = tan(rand() * pi);
                X1 = Best(d)   - a1 * b1 * abs(Best(d)   - xi(d));
                a2 = Beta * rand() - Beta;  b2 = tan(rand() * pi);
                X2 = Second(d) - a2 * b2 * abs(Second(d) - xi(d));
                xn(d) = (X1 + X2) / 2;
            end
        end

        xn = fix_candidates(xn, lb, ub, opts);
        for k = 1:numel(opts.operators)
            op = opts.operators{k};
            if rand() < op.prob
                xn = fix_candidates(op.fn(xn, progress), lb, ub, opts);
            end
        end
        fn = eval_one(xn, fobj, M);

        if dominates(fn, F(i, :))
            Xadv(i, :) = X(i, :);  Fadv(i, :) = F(i, :);
            X(i, :) = xn;          F(i, :) = fn;
        else
            Xadv(i, :) = xn;       Fadv(i, :) = fn;
        end
    end

    Xall = [X; Xadv];  Fall = [F; Fadv];
    [X, F] = env_select(Xall, Fall, N);
    [Archive, ArchiveFitness, ArchiveCD] = update_archive(Archive, ArchiveFitness, Xall, Fall, opts.archiveSize, opts.archive);

    if opts.verbose
        fprintf('IMOWO iter %d/%d: archive = %d\n', t, Tl, size(Archive, 1));
    end
end

ParetoSet = Archive;
ParetoFront = ArchiveFitness;
end

function opts = set_default(opts, name, val)
if ~isfield(opts, name), opts.(name) = val; end
end

function X = fix_candidates(X, lb, ub, opts)
X = min(max(X, repmat(lb, size(X, 1), 1)), repmat(ub, size(X, 1), 1));
if opts.integer, X = round(X); end
end

function F = eval_all(X, fobj, M)
F = zeros(size(X, 1), M);
for i = 1:size(X, 1), F(i, :) = eval_one(X(i, :), fobj, M); end
end

function f = eval_one(x, fobj, M)
try
    f = real(fobj(x));
    f = f(:)';
    if numel(f) ~= M || any(~isfinite(f)), f = 1e6 * ones(1, M); end
catch
    f = 1e6 * ones(1, M);
end
end

function d = dominates(a, b)
d = all(a <= b) && any(a < b);
end

function [rank, fronts] = nd_sort(F)
n = size(F, 1);
le = true(n); lt = false(n);
for k = 1:size(F, 2)
    col = F(:, k);
    le = le & bsxfun(@le, col, col');
    lt = lt | bsxfun(@lt, col, col');
end
D = le & lt;
cnt = sum(D, 1)';
rank = zeros(n, 1); fronts = {};
r = 0; remaining = true(n, 1);
while any(remaining)
    r = r + 1;
    cur = find(remaining & cnt == 0);
    if isempty(cur), cur = find(remaining); end
    rank(cur) = r; fronts{r} = cur;
    remaining(cur) = false;
    cnt = cnt - sum(D(cur, :), 1)';
    cnt(~remaining) = inf;
end
end

function cd = crowding(F)
n = size(F, 1);
cd = zeros(n, 1);
if n <= 2, cd(:) = inf; return; end
for k = 1:size(F, 2)
    [v, idx] = sort(F(:, k));
    cd(idx(1)) = inf; cd(idx(end)) = inf;
    rng_k = v(end) - v(1);
    if rng_k <= 0, continue; end
    cd(idx(2:end-1)) = cd(idx(2:end-1)) + (v(3:end) - v(1:end-2)) / rng_k;
end
end

function [Xs, Fs] = env_select(X, F, N)
[~, fronts] = nd_sort(F);
sel = [];
for r = 1:numel(fronts)
    fr = fronts{r};
    if numel(sel) + numel(fr) <= N
        sel = [sel; fr];
    else
        cdv = crowding(F(fr, :));
        [~, o] = sort(cdv, 'descend');
        sel = [sel; fr(o(1:N - numel(sel)))];
        break;
    end
end
Xs = X(sel, :); Fs = F(sel, :);
end

function [A, AF, ACD] = update_archive(A, AF, X, F, maxSize, mode)
A = [A; X]; AF = [AF; F];
[AF, ia] = unique(AF, 'rows', 'stable'); A = A(ia, :);
[rank, ~] = nd_sort(AF);
A = A(rank == 1, :); AF = AF(rank == 1, :);
if size(A, 1) > maxSize
    if strcmpi(mode, 'cd')
        while size(A, 1) > maxSize
            cdv = crowding(AF);
            [~, w] = min(cdv);
            A(w, :) = []; AF(w, :) = [];
        end
    else
        keep = knn_truncate(AF, maxSize);
        A = A(keep, :); AF = AF(keep, :);
    end
end
ACD = crowding(AF);
end

function keep = knn_truncate(F, maxSize)
lo = min(F, [], 1); rg = max(F, [], 1) - lo; rg(rg <= 0) = 1;
Z = bsxfun(@rdivide, bsxfun(@minus, F, lo), rg);
s = sum(Z.^2, 2);
D = sqrt(max(bsxfun(@plus, s, s') - 2 * (Z * Z'), 0));
n = size(F, 1);
D(1:n+1:end) = inf;
keep = (1:n)';
while numel(keep) > maxSize
    m = size(D, 1);
    [d1, i1] = min(D, [], 2);
    D2 = D; D2(sub2ind([m m], (1:m)', i1)) = inf;
    d2 = min(D2, [], 2);
    cand = find(d1 <= min(d1) + 1e-12);
    [~, q] = min(d2(cand));
    w = cand(q);
    keep(w) = []; D(w, :) = []; D(:, w) = [];
end
end

function [Best, Second] = pick_leaders(A, ACD, X, rk, mode)
if strcmpi(mode, 'first')
    order = [find(rk == 1); find(rk ~= 1)];
    Best = X(order(1), :);
    Second = X(order(min(2, numel(order))), :);
    return;
end
n = size(A, 1);
if n == 1, Best = A(1, :); Second = A(1, :); return; end
w = ACD;
fin = w(isfinite(w));
if isempty(fin) || max(fin) <= 0
    w(:) = 1;
else
    w(~isfinite(w)) = 2 * max(fin);
    w = w + eps;
end
i1 = roulette(w);
w2 = w; w2(i1) = 0;
i2 = roulette(w2);
Best = A(i1, :); Second = A(i2, :);
end

function i = roulette(w)
c = cumsum(w(:)) / sum(w);
i = find(rand() <= c, 1, 'first');
if isempty(i), i = numel(w); end
end

function s = levy_step(d)
beta = 1.5;
sigma = (gamma(1+beta) * sin(pi*beta/2) / (gamma((1+beta)/2) * beta * 2^((beta-1)/2)))^(1/beta);
s = (randn(1, d) * sigma) ./ abs(randn(1, d)).^(1/beta);
end

function h = halton_val(index, base)
h = 0; f = 1 / base; i = index;
while i > 0
    h = h + f * mod(i, base);
    i = floor(i / base);
    f = f / base;
end
end

function p = first_primes(k)
p = primes(max(10, 15 * k));
while numel(p) < k, p = primes(2 * max(p)); end
p = p(1:k);
end
