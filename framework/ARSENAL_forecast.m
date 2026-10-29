function [forecastUtil, W, rmse, mae] = ARSENAL_forecast(historyUtil, K)
% ARSENAL load forecast of the fog nodes.

if nargin < 2 || isempty(K)
    K = 5;
end

[m, T] = size(historyUtil);
K = min(K, T-1);
if K < 1
    forecastUtil = historyUtil(:, end)';
    W = 1;
    rmse = 0; mae = 0;
    return;
end

numCandidates = 200;
bestW = exp(-(0:K-1)') ;
bestW = bestW / sum(bestW);
bestMAE = compute_mae(historyUtil, bestW, K);

candW = [1; zeros(K-1, 1)];
candMAE = compute_mae(historyUtil, candW, K);
if candMAE < bestMAE
    bestMAE = candMAE;
    bestW = candW;
end

for c = 1:numCandidates
    alpha = linspace(3, 0.5, K)';
    g = sample_gamma(alpha);
    candW = g / sum(g);

    candMAE = compute_mae(historyUtil, candW, K);
    if candMAE < bestMAE
        bestMAE = candMAE;
        bestW = candW;
    end
end

W = bestW;

forecastUtil = zeros(1, m);
for j = 1:m
    recent = historyUtil(j, end:-1:end-K+1);
    forecastUtil(j) = sum(W' .* recent);
end
forecastUtil = max(0, min(1, forecastUtil));

[rmse, mae] = evaluate_accuracy(historyUtil, W, K);

end
