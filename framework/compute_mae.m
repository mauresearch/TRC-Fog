function mae = compute_mae(historyUtil, W, K)
    [m, T] = size(historyUtil);
    errs = [];
    for j = 1:m
        for t = (K+1):T
            pred = sum(W' .* historyUtil(j, t-1:-1:t-K));
            errs(end+1) = abs(historyUtil(j,t) - pred);
        end
    end
    if isempty(errs)
        mae = 0;
    else
        mae = mean(errs);
    end
end
