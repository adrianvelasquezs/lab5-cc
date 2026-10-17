function out = falsaPosicion(f, a, b, tolx, tolf, maxit)
% FALSAPOSICION  Regula falsi con criterio de parada doble (paso y residuo).
%   out = falsaPosicion(f, a, b, tolx, tolf, maxit) busca una raíz de f en [a, b].
%   Estados: "convergio", "intervalo_invalido", "max_iter".
%
%   Correcciones respecto a la versión de laboratorios anteriores:
%   - test de signo con sign() y reutilización de f en los extremos;
%   - el historial guarda a y b para poder detectar un extremo que queda fijo;
%   - se exige paso Y residuo pequeños.
fa = f(a); fb = f(b); nf = 2;
out = iniciarSalida(nf);
nombres = ["n","a","b","x","f(x)","Ea","Er"];
if ~(isfinite(fa) && isfinite(fb)) || sign(fa)*sign(fb) > 0
    out.estado = "intervalo_invalido"; return
end
if fa == 0 || fb == 0                       % un extremo ya es raíz
    xr = a; if fb == 0, xr = b; end
    out = cerrarSalida(out, zeros(0, 7), xr, 0, 0, 0, nf, nombres);
    out.estado = "convergio"; return
end
H = nan(maxit, 7); xold = NaN;
for n = 1:maxit
    x = b - fb*(b - a)/(fb - fa);
    fx = f(x); nf = nf + 1;
    Ea = abs(x - xold); Er = Ea/max(abs(x), eps);
    H(n,:) = [n a b x fx Ea Er];
    if sign(fa)*sign(fx) < 0
        b = x; fb = fx;
    else
        a = x; fa = fx;
    end
    if fx == 0 || (Ea <= tolx && abs(fx) <= tolf)
        out.estado = "convergio"; break
    end
    xold = x;
end
out = cerrarSalida(out, H(1:n,:), x, fx, Ea, n, nf, nombres);
end
