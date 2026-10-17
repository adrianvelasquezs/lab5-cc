function out = newtonRaphson(f, df, x0, tolx, tolf, maxit, tolder)
% NEWTONRAPHSON  Método de Newton-Raphson con control de derivada y de divergencia.
%   out = newtonRaphson(f, df, x0, tolx, tolf, maxit, tolder)
%   Estados: "convergio", "derivada_pequena" (|f'(x_n)| < tolder), "divergio", "max_iter".
%   out.nevalf y out.nevaldf cuentan por separado las evaluaciones de f y f'.
%
%   Correcciones respecto a la versión de laboratorios anteriores:
%   - se verifica |f'(x_n)| ANTES de dividir (evita pasos de magnitud ~1/eps o Inf);
%   - se reutiliza f(x_{n+1}) en la siguiente iteración (una evaluación de f por paso);
%   - se exige paso Y residuo pequeños y se detecta divergencia.
x = x0; fx = f(x); nf = 1; nd = 0; Ea = NaN;
out = iniciarSalida(nf);
H = nan(maxit, 8); n = 0;
while n < maxit
    dfx = df(x); nd = nd + 1;
    if abs(dfx) < tolder
        out.estado = "derivada_pequena"; break
    end
    n = n + 1;
    xn = x - fx/dfx;
    if ~isfinite(xn) || abs(xn) > 1e8
        H(n,:) = [n x fx dfx xn NaN NaN NaN];
        x = xn; fx = NaN; out.estado = "divergio"; break
    end
    fxn = f(xn); nf = nf + 1;
    Ea = abs(xn - x); Er = Ea/max(abs(xn), eps);
    H(n,:) = [n x fx dfx xn abs(fxn) Ea Er];
    x = xn; fx = fxn;
    if Ea <= tolx && abs(fx) <= tolf
        out.estado = "convergio"; break
    end
end
out = cerrarSalida(out, H(1:n,:), x, fx, Ea, n, nf, ...
    ["n","x_n","f(x_n)","f'(x_n)","x_{n+1}","r_{n+1}","Ea","Er"]);
out.nevaldf = nd;
end
