function T = filaResumen(caso, o, rRef, raices)
% FILARESUMEN  Fila de comparación de una ejecución de un método para f(x) = 0.
%   T = filaResumen(caso, o, rRef, raices) usa la salida o de un método, la raíz objetivo
%   rRef y las raíces de referencia ordenadas (negativa, central, positiva).
%   Evals: evaluaciones de f (o "f + f'" en Newton). ErrorRef = |raiz - rRef|.
[~, k] = min(abs(raices - o.raiz));
etiquetas = ["negativa", "central", "positiva"];
alcanzada = etiquetas(k);
if ~isfinite(o.raiz) || abs(o.raiz) > 1e6 || o.estado ~= "convergio"
    alcanzada = "ninguna";
end
evals = string(o.nevalf);
if o.nevaldf > 0
    evals = sprintf("%d+%d", o.nevalf, o.nevaldf);
end
T = table(string(caso), o.raiz, o.residuo, o.errPaso, o.iter, evals, o.estado, ...
    abs(o.raiz - rRef), alcanzada, 'VariableNames', {'Caso','Raiz','Residuo','ErrorPaso', ...
    'Iter','Evals','Estado','ErrorRef','RaizAlcanzada'});
end
