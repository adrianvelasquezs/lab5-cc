function out = cerrarSalida(out, H, x, fx, Ea, n, nf, nombres)
% CERRARSALIDA  Completa la estructura de salida y arma la tabla de historial.
out.raiz = x; out.residuo = abs(fx); out.errPaso = Ea;
out.iter = n; out.nevalf = nf;
out.hist = array2table(H, 'VariableNames', nombres(1:size(H,2)));
end
