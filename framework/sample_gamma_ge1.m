function x = sample_gamma_ge1(a)
    d = a - 1/3;
    c = 1 / sqrt(9*d);
    while true
        v = -1;
        while v <= 0
            z = randn();
            v = (1 + c*z)^3;
        end
        u = rand();
        if log(u) < 0.5*z^2 + d - d*v + d*log(v)
            x = d*v;
            return;
        end
    end
end
