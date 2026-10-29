function idx = compromise_point(F, fref)
% Compromise schedule of a Pareto front.
s = fref(:)';
s(~isfinite(s) | s <= 0) = 1;
Z = bsxfun(@rdivide, F, s);
[~, idx] = min(sqrt(sum(Z.^2, 2)));
end
