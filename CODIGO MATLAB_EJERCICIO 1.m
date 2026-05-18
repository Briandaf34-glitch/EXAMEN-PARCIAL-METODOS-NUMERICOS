% INFORME - EXAMEN PARCIAL 
% =========================================================================
clc; clear; close all;
format longG; % Configuración para visualizar con máxima precisión numérica 

% =========================================================================
% DATOS INICIALES DEL EXPERIMENTO (Nodos de control a 37 °C)
% =========================================================================
f = [100; 120; 145; 170; 200; 235; 270; 310; 355; 405; 460; 520; 585; ...
     655; 730; 810; 895; 985; 1080; 1180; 1290; 1410; 1540; 1680; 1830; ...
     1990; 2160; 2340; 2530; 2730];

Z = [152.3; 149.1; 146.8; 144.9; 142.0; 139.5; 137.9; 136.1; 134.8; 133.6; ...
     132.7; 131.9; 131.4; 131.1; 130.9; 131.0; 131.3; 131.9; 132.7; 133.8; ...
     135.2; 136.9; 138.9; 141.1; 143.5; 146.1; 149.0; 152.2; 155.6; 159.2]; 

n = length(f);

% =========================================================================
% PARTE A: ANÁLISIS EXPLORATORIO DE DATOS
% =========================================================================
figure('Name', 'Parte A - Análisis Exploratorio');
plot(f, Z, 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
grid on;
title('Análisis Exploratorio: Magnitud de Impedancia |Z| vs Frecuencia f'); 
xlabel('Frecuencia f (Hz)');
ylabel('Impedancia |Z| (\Omega)');

% Identificación visual del mínimo indexando el valor más bajo medido
[valMinVisual, idxMinVisual] = min(Z);
fprintf('--- PARTE A: MÍNIMO VISUAL ---\n');
fprintf('Frecuencia aproximada del mínimo: %d Hz\n', f(idxMinVisual)); 
fprintf('Valor de impedancia mínimo observado: %.4f Ohm\n\n', valMinVisual);

% =========================================================================
% PARTE B1: INTERPOLACIÓN POLINÓMICA Y VALIDACIÓN LOO
% =========================================================================
% Ajuste de un polinomio estable de grado 5 para ilustrar la mitigación de Runge
grado_pol = 5; 
p5 = polyfit(f, Z, grado_pol); 

% Evaluar la interpolación en el punto de interés f = 1000 Hz
Z_1000_poly = polyval(p5, 1000); 

% --- Procedimiento Leave-One-Out (LOO) ---
rng(42); % Fijar semilla aleatoria para garantizar reproducibilidad
puntos_azar = randperm(n, 5); % Selección de 5 índices al azar 
errores_loo = zeros(5, 1);

for i = 1:5
    idx_omito = puntos_azar(i);
    
    % Construir conjuntos temporales omitiendo el punto i
    f_loo = f; f_loo(idx_omito) = [];
    Z_loo = Z; Z_loo(idx_omito) = [];
    
    % Ajustar polinomio con los 29 puntos restantes
    p_temp = polyfit(f_loo, Z_loo, grado_pol);
    
    % Predecir el valor en el nodo omitido y calcular error relativo
    z_predicho = polyval(p_temp, f(idx_omito));
    errores_loo(i) = abs((Z(idx_omito) - z_predicho) / Z(idx_omito)); 
end
error_loo_promedio = mean(errores_loo) * 100;

fprintf('--- PARTE B1: INTERPOLACIÓN POLINÓMICA ---\n');
fprintf('Valor interpolado polinómico en f = 1000 Hz: %.4f Ohm\n', Z_1000_poly); 
fprintf('Error relativo estimado por LOO (Promedio): %.4f%%\n\n', error_loo_promedio); 

% =========================================================================
% PARTE B2: SPLINES CÚBICOS NATURALES
% =========================================================================
% Construcción del spline cúbico natural ('variational' define s''(x) = 0 en extremos)
spline_fit = csape(f, Z, 'variational'); 

% Evaluación en una malla fina para análisis comparativo visual
f_fina = linspace(min(f), max(f), 1000); 
Z_spline_fina = fnval(spline_fit, f_fina);
Z_poly_fina = polyval(p5, f_fina);

% Graficar comparación para la discusión técnica del informe
figure('Name', 'Parte B2 - Comparativa Polinomio vs Spline');
plot(f, Z, 'ro', 'DisplayName', 'Datos Originales'); hold on;
plot(f_fina, Z_poly_fina, 'b--', 'LineWidth', 1.2, 'DisplayName', 'Polinomio Grado 5');
plot(f_fina, Z_spline_fina, 'g-', 'LineWidth', 1.8, 'DisplayName', 'Spline Cúbico Natural'); 
grid on; legend('Location', 'best');
title('Comparativa de Modelos: Estabilidad frente a Oscilaciones');
xlabel('Frecuencia f (Hz)'); ylabel('Impedancia |Z| (\Omega)');

% Evaluar Spline en f = 1000 Hz
Z_1000_spline = fnval(spline_fit, 1000); 

fprintf('--- PARTE B2: SPLINES CÚBICOS ---\n');
fprintf('Valor interpolado por Spline en f = 1000 Hz: %.4f Ohm\n\n', Z_1000_spline); 
% =========================================================================
% PARTE C: DERIVACIÓN NUMÉRICA Y ESTABILIDAD DEL MÍNIMO
% =========================================================================
% Derivación analítica del objeto spline (Evita errores de diferencias finitas)
spline_deriv1 = fnder(spline_fit, 1); 
spline_deriv2 = fnder(spline_fit, 2); 

% Evaluar la primera derivada en la malla fina para graficar
Z_deriv1_fina = fnval(spline_deriv1, f_fina); 
figure('Name', 'Parte C - Primera Derivada');
plot(f_fina, Z_deriv1_fina, 'k-', 'LineWidth', 1.5); hold on;
plot(f, fnval(spline_deriv1, f), 'bo', 'DisplayName', 'Nodos'); 
grid on;
title('Primera Derivada Analítica d|Z|/df vs Frecuencia'); 
xlabel('Frecuencia f (Hz)'); ylabel('d|Z|/df (\Omega/Hz)');

% Encontrar la raíz exacta de la primera derivada (donde d|Z|/df = 0)
% Usamos como aproximación inicial el mínimo visual de 730 Hz
f_min_exacto = fzero(@(x) fnval(spline_deriv1, x), 730); 
Z_min_exacto = fnval(spline_fit, f_min_exacto);

% Evaluar el signo de la segunda derivada en ese punto crítico
d2Z_critico = fnval(spline_deriv2, f_min_exacto); 

fprintf('--- PARTE C: DERIVACIÓN Y PUNTO CRÍTICO ---\n');
fprintf('Frecuencia exacta del mínimo local: %.4f Hz\n', f_min_exacto); 
fprintf('Valor de impedancia mínimo calculado: %.4f Ohm\n', Z_min_exacto);
fprintf('Segunda derivada en el mínimo (d^2|Z|/df^2): %.4f Ohm/Hz^2\n', d2Z_critico);

% =========================================================================
% PARTE D: BÚSQUEDA DE RAÍCES 
% =========================================================================
Z_th = 150; 
% Definición de la función objetivo: |Z|(f) - 150 = 0
obj_raiz = @(x) fnval(spline_fit, x) - Z_th; 
deriv_raiz = @(x) fnval(spline_deriv1, x); 

tolerancia = 1e-5; % Margen para asegurar la convergencia con alta precisión 

% --- RAÍZ 1: Zona de baja frecuencia (Intervalo [100, 200] Hz) ---
[r1_bis, it1_bis] = mi_biseccion(obj_raiz, 100, 200, tolerancia);
[r1_nr, it1_nr]   = newton_raphson(obj_raiz, deriv_raiz, 100, tolerancia);

% --- RAÍZ 2: Zona de alta frecuencia (Intervalo [2000, 2400] Hz) ---
[r2_bis, it2_bis] = mi_biseccion(obj_raiz, 2000, 2400, tolerancia);
[r2_nr, it2_nr]   = newton_raphson(obj_raiz, deriv_raiz, 2300, tolerancia);

dfdZ_sensibilidad = 1 / fnval(spline_deriv1, r2_nr);

fprintf('--- PARTE D: BÚSQUEDA DE RAÍCES (Umbral = 150 Ohm) ---\n');
fprintf('Raíz 1 (Baja Frecuencia):\n');
fprintf('  > Bisección:      %.4f Hz  Iteraciones: %d\n', r1_bis, it1_bis);
fprintf('  > Newton-Raphson: %.4f Hz  Iteraciones: %d\n', r1_nr, it1_nr);
fprintf('Raíz 2 (Alta Frecuencia):\n');
fprintf('  > Bisección:      %.4f Hz  Iteraciones: %d\n', r2_bis, it2_bis);
fprintf('  > Newton-Raphson: %.4f Hz  Iteraciones: %d\n', r2_nr, it2_nr);
fprintf('Análisis de Sensibilidad en Raíz 2 (cercana a 2200 Hz):\n');
fprintf('  > Derivada inversa df/d|Z|: %.4f Hz/Ohm\n\n', dfdZ_sensibilidad); 

fprintf('--- PARTE D: BÚSQUEDA DE RAÍCES (Umbral = 150 Ohm) ---\n'); 
fprintf('Raíz 1 (Baja Frecuencia):\n');
fprintf('  > Bisección:      %.4f Hz  Iteraciones: %d\n', r1_bis, it1_bis); 
fprintf('  > Newton-Raphson: %.4f Hz  Iteraciones: %d\n', r1_nr, it1_nr); 
fprintf('Raíz 2 (Alta Frecuencia):\n');
fprintf('  > Bisección:      %.4f Hz  Iteraciones: %d\n', r2_bis, it2_bis); 
fprintf('  > Newton-Raphson: %.4f Hz  Iteraciones: %d\n', r2_nr, it2_nr); 
fprintf('Análisis de Sensibilidad en Raíz 2 (cercana a 2200 Hz):\n'); 
fprintf('  > Derivada inversa df/d|Z|: %.4f Hz/Ohm\n\n', dfdZ_sensibilidad); 

% =========================================================================
% LAS FUNCIONES AUXILIARES 
% =========================================================================
function [c, iter] = mi_biseccion(func, a, b, tol)
iter = 0;
while (b - a)/2 > tol
    iter = iter + 1;
    c = (a + b) / 2;
    if func(c) == 0
        return;
    elseif func(a)*func(c) < 0
        b = c;
    else
        a = c;
    end
end
c = (a + b) / 2;
end

function [x1, iter] = newton_raphson(func, deriv, x0, tol)
iter = 0;
max_iter = 100;
x1 = x0;
error = 1;
while error > tol && iter < max_iter
    iter = iter + 1;
    y = func(x1);
    dy = deriv(x1);
    if dy == 0
        error('Derivada cero encontrada.');
    end
    x_new = x1 - y/dy;
    error = abs(x_new - x1);
    x1 = x_new;
end
end
