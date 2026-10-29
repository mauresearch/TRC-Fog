function g = sample_gamma(alphaVec)
    g = zeros(size(alphaVec));
    for idx = 1:numel(alphaVec)
        a = alphaVec(idx);
        if a < 1
            g(idx) = sample_gamma_ge1(a+1) * rand()^(1/a);
        else
            g(idx) = sample_gamma_ge1(a);
        end
    end
end
