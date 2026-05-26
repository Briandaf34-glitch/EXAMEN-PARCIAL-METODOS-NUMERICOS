% EXAMEN PARCIAL DE MÉTODOS NUMÉRICOS
% EJERCICIO 2
clear; clc; format shortG;

% --- DATOS DEL ENSAYO ---
f = [10.0:2.5:107.5]; % Vector de frecuencias de 10 a 107.5
V = [0.842, 0.911, 0.986, 1.062, 1.143, 1.227, 1.314, 1.401, 1.482, 1.551, ...
    1.216, 1.048, 0.866, 0.689, 0.521, 0.364, 0.223, 0.103, 0.012, -0.041, ...
    -0.057,-0.034, 0.018, 0.096, 0.197, 0.318, 0.452, 0.579, 0.700, 0.809, ...
    0.611, 0.688, 0.756, 0.811, 0.856, 0.894, 0.926, 0.954, 0.980, 1.004, 1.02, 1.03]; 
% Ajuste de dimensiones a los 40 puntos de las tablas proporcionadas
V = V(1:40); 
Z = [182.4, 178.9, 175.1, 171.0, 166.8, 162.7, 158.9, 155.4, 152.0, 149.0, ...
    146.1, 145.2, 145.8, 147.3, 149.9, 153.5, 158.0, 163.2, 168.9, 174.8, ...
    180.5, 186.2, 191.5, 196.2, 200.1, 203.1, 205.2, 206.3, 206.1, 204.7, ...
    198.0, 194.4, 190.9, 187.8, 185.1, 183.0, 181.6, 180.8, 180.6, 180.9];

% --- PARTE 1: INTERPOLACIÓN ---
fprintf('--- PARTE 1: INTERPOLACIÓN ---\n');
f_obj = [41.0, 73.0];
disp('Resultados usando Spline Cúbico:');
V_spline = spline(f, V, f_obj);
Z_spline = spline(f, Z, f_obj);
fprintf('V(41.0) = %.4f V | Z(41.0) = %.4f Ohm\n', V_spline(1), Z_spline(1));
fprintf('V(73.0) = %.4f V | Z(73.0) = %.4f Ohm\n\n', V_spline(2), Z_spline(2));

% --- PARTE 2: DERIVACIÓN NUMÉRICA ---
fprintf('--- PARTE 2: DERIVACIÓN ---\n');
h = 2.5;

% Diferencia progresiva O(h^2) para f=10
df10_prog = (-3*V(1) + 4*V(2) - V(3)) / (2*h);
fprintf('dV/df en 10.0 kHz (Progresiva O2): %.4f\n', df10_prog);

% Puntos a evaluar y sus índices en el arreglo:
puntos = [40.0, 70.0, 100.0];
indices = [13, 25, 37]; % Índices correspondientes a esas frecuencias en el arreglo

for i = 1:length(indices)
    idx = indices(i);
    % O(h^2) Centrada
    df_o2 = (V(idx+1) - V(idx-1)) / (2*h);
    % O(h^4) Centrada
    df_o4 = (-V(idx+2) + 8*V(idx+1) - 8*V(idx-1) + V(idx-2)) / (12*h);

    fprintf('dV/df en %.1f kHz -> Ord 2: %.4f | Ord 4: %.4f\n', puntos(i), df_o2, df_o4);
end

% Comparación con la derivada del spline
% Creamos la estructura del spline y usamos la función de derivación
spline_V = csapi(f, V); % Requiere Curve Fitting Toolbox
spline_deriv = fnder(spline_V, 1);
df_spline = fnval(spline_deriv, [10.0, 40.0, 70.0, 100.0]);
fprintf('\nDerivadas usando Spline: 10kHz=%.4f, 40kHz=%.4f, 70kHz=%.4f, 100kHz=%.4f\n\n', df_spline);


% --- PARTE 3: RAÍCES Y BISECCIÓN ---
fprintf('--- PARTE 3: BÚSQUEDA DE RAÍCES ---\n');
% Encontrar intervalos de cambio de signo
cambios = find(V(1:end-1) .* V(2:end) < 0);

for i = 1:length(cambios)
    f_izq = f(cambios(i));
    f_der = f(cambios(i)+1);
    fprintf('Raíz %d detectada en el intervalo [%.1f, %.1f]\n', i, f_izq, f_der);

    % Bisección usando el spline como función continua
    a = f_izq; b = f_der;
    tol = 1e-4; max_iter = 50; iter = 0;
    while (b - a)/2 > tol && iter < max_iter
        c = (a + b) / 2;
        % Evaluamos el spline en 'c' y en 'a'
        if spline(f, V, c) == 0
            break;
        elseif spline(f, V, a) * spline(f, V, c) < 0
            b = c;
        else
            a = c;
        end
        iter = iter + 1;
    end
    fprintf(' -> Raíz estimada por bisección (Spline): %.4f kHz\n', (a+b)/2);
end
