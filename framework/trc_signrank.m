function p = trc_signrank(a, b)
% Paired Wilcoxon signed-rank test.
d = a(:) - b(:);
d = d(isfinite(d));
d = d(abs(d) > 1e-12);
n = numel(d);
if n < 5, p = NaN; return; end
[s, idx] = sort(abs(d));
r = zeros(n, 1);
k = 1; tieTerm = 0;
while k <= n
    j = k;
    while j < n && abs(s(j+1) - s(k)) < 1e-12, j = j + 1; end
    r(k:j) = (k + j) / 2;
    t = j - k + 1; tieTerm = tieTerm + (t^3 - t);
    k = j + 1;
end
rk = zeros(n, 1); rk(idx) = r;
Wp = sum(rk(d > 0));
mu = n * (n + 1) / 4;
sg = sqrt(n * (n + 1) * (2*n + 1) / 24 - tieTerm / 48);
z = (Wp - mu) / max(sg, eps);
p = erfc(abs(z) / sqrt(2));
end
