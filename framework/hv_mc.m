function hv = hv_mc(F, ref, nSamples, seed)
% Monte-Carlo hypervolume.
if nargin < 3, nSamples = 100000; end
if nargin < 4, seed = 12345; end
M = size(F, 2); ref = ref(:)';
if numel(ref) == 1, ref = repmat(ref, 1, M); end
F = F(all(bsxfun(@lt, F, ref), 2), :);
if isempty(F), hv = 0; return; end
s0 = rng(); rng(seed, 'twister');
lo = min(0, min(F, [], 1));
U = bsxfun(@plus, lo, bsxfun(@times, rand(nSamples, M), ref - lo));
rng(s0);
dom = false(nSamples, 1);
for i = 1:size(F, 1)
    dom = dom | all(bsxfun(@ge, U, F(i, :)), 2);
end
hv = mean(dom) * prod(ref - lo);
end
