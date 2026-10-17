function [p, razon] = ordenConvergencia(xs, r)
% ORDENCONVERGENCIA  Orden p y razón media de convergencia de una sucesión xs hacia r.
%   p: pendiente (mínimos cuadrados) de ln e_{n+1} frente a ln e_n en los últimos pares.
%   razon: media geométrica de e_{n+1}/e_n (constante asintótica si el método es lineal).
e = abs(xs(:) - r);
e = e(1:find(e > 1e-13, 1, 'last'));  % descarta la cola dominada por el redondeo
if numel(e) < 3, p = NaN; razon = NaN; return, end
u = log(e(1:end-1));  v = log(e(2:end));
k = max(1, numel(u) - 5):numel(u);    % últimos pares (régimen asintótico)
coef = polyfit(u(k), v(k), 1);
p = coef(1);
razon = exp(mean(v - u));
end
