function out = biseccion(f, a, b, tolx, tolf, maxit)
% BISECCION  Método de bisección con criterio de parada doble (paso y residuo).
%   out = biseccion(f, a, b, tolx, tolf, maxit) busca una raíz de f en [a, b].
%   Declara "convergio" si |x_n - x_{n-1}| <= tolx y |f(x_n)| <= tolf.
%   Estados: "convergio", "intervalo_invalido", "max_iter".
%
%   Correcciones respecto a la versión de laboratorios anteriores:
%   - el test de signo usa sign(fa)*sign(fx) en vez de fa*fx (evita underflow/overflow);
%   - el punto medio se calcula como a + (b - a)/2;
%   - se reutilizan f(a), f(b) y se cuentan exactamente las evaluaciones de f;
%   - se exige paso Y residuo pequeños (un paso pequeño no implica convergencia).
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
    x = a + (b - a)/2;
    fx = f(x); nf = nf + 1;
    Ea = abs(x - xold); Er = Ea/max(abs(x), eps);
    H(n,:) = [n a b x fx Ea Er];
    if sign(fa)*sign(fx) < 0
        b = x;
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
