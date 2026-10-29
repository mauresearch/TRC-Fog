function [rmse, mae] = evaluate_accuracy(historyUtil, W, K)
    [m, T] = size(historyUtil);
    sqErrs = [];
    absErrs = [];
    for j = 1:m
        for t = (K+1):T
            pred = sum(W' .* historyUtil(j, t-1:-1:t-K));
            e = historyUtil(j,t) - pred;
            sqErrs(end+1) = e^2;
            absErrs(end+1) = abs(e);
        end
    end
    if isempty(sqErrs)
        rmse = 0; mae = 0;
    else
        rmse = sqrt(mean(sqErrs));
        mae = mean(absErrs);
    end
end
