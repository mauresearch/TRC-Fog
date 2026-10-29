function AF = arrhenius_af(Tcelsius, P)
% Arrhenius acceleration factor relative to the reference temperature.
kB = 8.617e-5;
TK   = Tcelsius + 273.15;
Tref = P.Tref + 273.15;
AF = exp((P.Ea / kB) .* (1 ./ Tref - 1 ./ TK));
end
