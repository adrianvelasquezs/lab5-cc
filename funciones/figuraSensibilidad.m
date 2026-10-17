function fig = figuraSensibilidad(f, df, x0N, oN, oSA, oSB, oG1, oG2, raicesRef, rRef)
% FIGURASENSIBILIDAD  Figura 2 del punto 2: sensibilidad y fallas por datos iniciales.
%   (a) tangentes/secante casi horizontales cerca de f' = 0; (b) trayectorias x_n;
%   (c) residuos de las semillas sensibles; (d) punto fijo g1 (contractiva) vs g2 (repulsiva).
%   oN{k} son las salidas de Newton para las semillas x0N(k); los residuos nulos se
%   dibujan en 1e-16 para poder usar escala logarítmica.
piso = 1e-16;
fig = figure('Name', 'Figura 2', 'Units', 'pixels');
fig.Position(3:4) = [900 680];          % más espacio para los cuatro paneles
theme(fig, 'light');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact');

nexttile
xa = linspace(-2.6, 2.6, 400);
plot(xa, f(xa), 'LineWidth', 1.5); hold on; grid on; yline(0, '-');
for k = 2:3
    plot(xa, f(x0N(k)) + df(x0N(k))*(xa - x0N(k)), '--', 'LineWidth', 1.2);
end
plot(xa, f(1.0) + (f(1.3) - f(1.0))/0.3*(xa - 1.0), ':', 'LineWidth', 1.5);
xlabel('x'); ylabel('f_N(x)'); ylim([-8 8]);
legend('f_N', '', sprintf('tangente x_0 = %.2f', x0N(2)), sprintf('tangente x_0 = %.2f', x0N(3)), ...
    'secante (1.0, 1.3)', 'Location', 'south', 'FontSize', 7);
title('(a) Tangente/secante casi horizontal'); hold off

nexttile
for k = 1:3
    plot(0:oN{k}.iter, [x0N(k); oN{k}.hist.("x_{n+1}")], '.-', 'LineWidth', 1.2, 'MarkerSize', 12);
    hold on
end
plot(0:oSB.iter+1, [1.0; 1.3; oSB.hist.("x_{n+1}")], 'v-', 'LineWidth', 1.2);
grid on; yline(raicesRef, ':');
xlabel('iteración n'); ylabel('x_n');
legend([compose("Newton %.2f", x0N(1:3)), "Secante (1.0, 1.3)"], 'Location', 'southeast', 'FontSize', 7);
title('(b) Trayectorias de los iterados'); hold off

nexttile
for k = 1:3
    semilogy(oN{k}.hist.n, max(oN{k}.hist.("r_{n+1}"), piso), '.-', 'LineWidth', 1.2, 'MarkerSize', 12);
    hold on
end
semilogy(oSA.hist.n, max(oSA.hist.("r_{n+1}"), piso), '^-', 'LineWidth', 1.2);
semilogy(oSB.hist.n, max(oSB.hist.("r_{n+1}"), piso), 'v-', 'LineWidth', 1.2);
grid on; xlabel('iteración n'); ylabel('r_n');
legend('Newton 2.5', 'Newton 1.16', 'Newton 1.15 (otra raíz)', 'Secante A', 'Secante B', ...
    'Location', 'southwest', 'FontSize', 7);
title('(c) Residuo de los casos sensibles'); hold off

nexttile
semilogy(0:oG1.iter, abs([oG1.hist.x_n(1); oG1.hist.("g(x_n)")] - rRef), '^-', 'LineWidth', 1.2);
hold on; grid on
semilogy(0:oG2.iter, abs([oG2.hist.x_n(1); oG2.hist.("g(x_n)")] - rRef), 'x-', 'LineWidth', 1.2);
xlabel('iteración n'); ylabel('|x_n - x^*|');
legend('g_1 (contractiva)', 'g_2 (repulsiva)', 'Location', 'west', 'FontSize', 7);
title('(d) Punto fijo desde x_0 = 2'); hold off
end
