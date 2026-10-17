function out = secante(f, x0, x1, tolx, tolf, maxit)
% SECANTE  Método de la secante con control de denominador nulo y de divergencia.
%   out = secante(f, x0, x1, tolx, tolf, maxit)
%   Estados: "convergio", "denominador_nulo" (f(x_n) = f(x_{n-1})), "divergio", "max_iter".
%
%   Correcciones respecto a la versión de laboratorios anteriores:
%   - se verifica el denominador f(x_n) - f(x_{n-1}) antes de dividir;
%   - una sola evaluación de f por iteración (se reutiliza la anterior);
%   - se exige paso Y residuo pequeños y se detecta divergencia.
f0 = f(x0); f1 = f(x1); nf = 2;
out = iniciarSalida(nf);
H = nan(maxit, 7); x = x1; fx = f1; Ea = NaN; n = 0;
while n < maxit
    if f1 == f0
        out.estado = "denominador_nulo"; break
    end
    n = n + 1;
    x = x1 - f1*(x1 - x0)/(f1 - f0);
    if ~isfinite(x) || abs(x) > 1e8
        H(n,:) = [n x0 x1 x NaN NaN NaN];
        fx = NaN; out.estado = "divergio"; break
    end
    fx = f(x); nf = nf + 1;
    Ea = abs(x - x1); Er = Ea/max(abs(x), eps);
    H(n,:) = [n x0 x1 x abs(fx) Ea Er];
    x0 = x1; f0 = f1; x1 = x; f1 = fx;
    if Ea <= tolx && abs(fx) <= tolf
        out.estado = "convergio"; break
    end
end
out = cerrarSalida(out, H(1:n,:), x, fx, Ea, n, nf, ...
    ["n","x_{n-1}","x_n","x_{n+1}","r_{n+1}","Ea","Er"]);
end
