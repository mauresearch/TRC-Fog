function f = non_domination_sort_mod(x, M, D)
% Non-dominated sorting and crowding distance (vectorised).
N = size(x, 1);
Fo = x(:, D+1:D+M);
le = true(N); eq = true(N);
for k = 1:M
    c = Fo(:, k);
    le = le & bsxfun(@ge, c, c');
    eq = eq & bsxfun(@eq, c, c');
end
B = le & ~eq;
cnt = sum(B, 2);
rank = zeros(N, 1);
cur = find(cnt == 0);
r = 1;
remaining = true(N, 1);
while ~isempty(cur)
    rank(cur) = r;
    remaining(cur) = false;
    cnt = cnt - sum(B(:, cur), 2);
    cnt(~remaining) = -1;
    cur = find(cnt == 0 & remaining);
    r = r + 1;
end
x(:, M + D + 1) = rank;
[~, ord] = sort(rank);
xs = x(ord, 1:M + D + 1);
f = zeros(N, M + D + 2);
pos = 0;
for fr = 1:max(rank)
    n = sum(rank == fr);
    y = xs(pos + 1:pos + n, :);
    dist = zeros(n, 1);
    for k = 1:M
        [v, idx] = sort(y(:, D + k));
        part = zeros(n, 1);
        part(idx(end)) = Inf;
        part(idx(1)) = Inf;
        if n > 2
            if v(end) - v(1) == 0
                part(idx(2:end-1)) = Inf;
            else
                part(idx(2:end-1)) = (v(3:end) - v(1:end-2)) / (v(end) - v(1));
            end
        end
        dist = dist + part;
    end
    f(pos + 1:pos + n, :) = [y, dist];
    pos = pos + n;
end
end
