function out = puntoFijo(g, f, x0, tolx, tolf, maxit)
% PUNTOFIJO  Iteración de punto fijo x_{n+1} = g(x_n).
%   out = puntoFijo(g, f, x0, tolx, tolf, maxit) vigila el paso |x_{n+1} - x_n| y el
%   residuo |f(x_{n+1})| de la ecuación original f(x) = 0.
%   Estados: "convergio", "divergio" (iterado no finito o |x| > 1e8), "max_iter".
%
%   Correcciones respecto a la versión de laboratorios anteriores:
%   - se vigila el residuo de f y no solo el paso (un paso pequeño no basta);
%   - se detecta divergencia en lugar de iterar hasta Inf/NaN.
out = iniciarSalida(0);
H = nan(maxit, 6); x = x0; nf = 0;
for n = 1:maxit
    gx = g(x);
    fx = f(gx); nf = nf + 1;
    Ea = abs(gx - x); Er = Ea/max(abs(gx), eps);
    H(n,:) = [n x gx abs(fx) Ea Er];
    x = gx;
    if ~isfinite(x) || abs(x) > 1e8
        out.estado = "divergio"; break
    end
    if Ea <= tolx && abs(fx) <= tolf
        out.estado = "convergio"; break
    end
end
out = cerrarSalida(out, H(1:n,:), x, fx, Ea, n, nf, ...
    ["n","x_n","g(x_n)","|f(g(x_n))|","Ea","Er"]);
end
