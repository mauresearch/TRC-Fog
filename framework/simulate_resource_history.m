function historyUtil = simulate_resource_history(numNodes, historyLength, seed)

if nargin < 2 || isempty(historyLength)
    historyLength = 20;
end
if nargin < 3 || isempty(seed)
    seed = 99;
end

s = rng();
rng(seed, 'twister');

t = 1:historyLength;
historyUtil = zeros(numNodes, historyLength);

for j = 1:numNodes
    baseline = 0.2 + 0.5*rand();
    trend = 0.1 * sin(2*pi*t/historyLength + rand()*pi);
    noise = 0.05 * randn(1, historyLength);
    u = baseline + trend + noise;
    historyUtil(j, :) = max(0, min(1, u));
end

rng(s);
end
