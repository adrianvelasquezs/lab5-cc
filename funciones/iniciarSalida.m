function out = iniciarSalida(nf)
% INICIARSALIDA  Estructura de salida común a todos los métodos para f(x) = 0.
%   Campos: raiz, residuo, errPaso, iter, nevalf, nevaldf, estado, hist (tabla).
out = struct('raiz', NaN, 'residuo', NaN, 'errPaso', NaN, 'iter', 0, ...
    'nevalf', nf, 'nevaldf', 0, 'estado', "max_iter", 'hist', table());
end
