%% TRABAJO MECÁNICA Y ONDAS I - APARTADO B.
% Autoras: Sara y Claudia
% Grado en Física - Universidad de Córdoba


clear; clc; close all;

%% 1. PARÁMETROS DEL SISTEMA Y CONDICIONES INICIALES
g = 9.81;       % Aceleración de la gravedad (m/s^2)
m = 100;        % Masa del vehículo (kg)
a = 0.25;       % Parámetro de curvatura z = a*x^2 (1/m) [Ajustado al dibujo]

x0 = 5;         % Posición horizontal inicial (m)
v0 = 0;         % Velocidad horizontal inicial (m/s)
y0 = [x0; v0];  % Vector de estado inicial

tspan = [0 15]; % Intervalo de tiempo para la simulación (s)

%% 2. INTEGRACIÓN NUMÉRICA CON ODE45
[t, Y] = ode45(@(t, y) ec_movimiento_b(t, y, g, a), tspan, y0);

x = Y(:, 1);    % Posición x(t)
vx = Y(:, 2);   % Velocidad v_x(t)
z = a * x.^2;   % Posición z(t) por la ligadura z = a*x^2

%% 3. CÁLCULOS FÍSICOS AVANZADOS (REACCIÓN NORMAL, FUERZAS G Y ENERGÍA)
% Multiplicador de Lagrange lambda(t) y componentes cartesianas de N
lambda = (m * (g + 2 * a * vx.^2)) ./ (1 + 4 * a^2 * x.^2);
Nx = -2 * a * x .* lambda;                   % Componente horizontal N_x
Nz = lambda;                                 % Componente vertical N_z
N_modulo = lambda .* sqrt(1 + 4 * a^2 * x.^2); % Módulo de la Normal |N|

% Fuerza G sentida por los pasajeros (N / (m*g))
fuerza_g = N_modulo / (m * g);

% Balance de Energía Mecánica
E_cinetica = 0.5 * m * (1 + 4 * a^2 * x.^2) .* vx.^2;
E_potencial = m * g * a * x.^2;
E_total = E_cinetica + E_potencial;

% Error relativo de conservación de energía: |E(t) - E(0)| / E(0)
error_energia = abs(E_total - E_total(1)) / E_total(1);


%% 4. FIGURA 1: ESTUDIO ESTÁTICO COMPLETO (6 SUBPLOTS DE RIGOR FÍSICO)
figure('Name', 'Estudio Dinámico y Análisis de Precisión - Caída en el Abismo Cóncavo', ...
       'NumberTitle', 'off', 'Color', 'w');

% Subplot 1: Posición y Velocidad vs Tiempo
subplot(3, 2, 1);
plot(t, x, 'b-', 'LineWidth', 1.5); hold on;
plot(t, vx, 'r--', 'LineWidth', 1.5);
xlabel('Tiempo (s)'); ylabel('Magnitud');
title('a) Posición x(t) y Velocidad v_x(t)');
legend('x (m)', 'v_x (m/s)'); grid on;

% Subplot 2: Espacio / Retrato de Fases (v_x vs x)
subplot(3, 2, 2);
plot(x, vx, 'm-', 'LineWidth', 1.5);
xlabel('Posición x (m)'); ylabel('Velocidad v_x (m/s)');
title('b) Retrato de Fase (v_x vs x)'); grid on;

% Subplot 3: Módulo de la Reacción Normal N(t)
subplot(3, 2, 3);
plot(t, N_modulo, 'k-', 'LineWidth', 1.5);
xlabel('Tiempo (s)'); ylabel('Fuerza Normal (N)');
title('c) Módulo de la Reacción Normal N(t)'); grid on;

% Subplot 4: Aceleración Soportada en "Fuerzas G"
subplot(3, 2, 4);
plot(t, fuerza_g, 'c-', 'LineWidth', 1.5); hold on;
yline(1, 'k--', '1g (Normal)');
xlabel('Tiempo (s)'); ylabel('Sostenido (g)');
title('d) Carga de Aceleración ("Fuerza G")'); grid on;

% Subplot 5: Conservación de la Energía Mecánica Total
subplot(3, 2, 5);
plot(t, E_total, 'g-', 'LineWidth', 1.5);
xlabel('Tiempo (s)'); ylabel('Energía (J)');
title('e) Energía Mecánica Total E(t)');
ylim([min(E_total)*0.95, max(E_total)*1.05]); grid on;

% Subplot 6: Error Relativo de la Energía (Rigor Numérico)
subplot(3, 2, 6);
plot(t, error_energia, 'r-', 'LineWidth', 1.5);
xlabel('Tiempo (s)'); ylabel('Error Relativo');
title('f) Error Numérico de Energía \DeltaE(t)/E_0'); grid on;


%% 5. FIGURA 2: MODELADO TRIDIMENSIONAL DE LA CANALETA (CANALETA 3D)
figure('Name', 'Vista Tridimensional de la Atracción (Canaleta 3D)', ...
       'NumberTitle', 'off', 'Color', 'w');

% Malla para la canaleta 3D
x_mesh = linspace(-1.2*abs(x0), 1.2*abs(x0), 100);
y_mesh = linspace(-2, 2, 30);
[X_mesh, Y_mesh] = meshgrid(x_mesh, y_mesh);
Z_mesh = a * X_mesh.^2;

% Dibujar superficie de la canaleta usando un colormap válido ('cool' o 'winter')
surf(X_mesh, Y_mesh, Z_mesh, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
colormap cool; hold on; 

% Dibujar la línea central de la vía y el punto inicial
plot3(x, zeros(size(x)), z, 'r-', 'LineWidth', 2.5); 
plot3(x(1), 0, z(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

grid on; axis equal; view(-37.5, 30);
xlabel('x (m)'); ylabel('y (m) [Ancho de la vía]'); zlabel('z (m)');
title('Estructura 3D de la Atracción "Caída en el Abismo Cóncavo"');

%% 6. FIGURA 3: ANIMACIÓN 2D Y GRABACIÓN A VÍDEO (.MP4)
figure('Name', 'Animación con Vector Normal Dinámico - Caída en la Parábola', ...
       'NumberTitle', 'off', 'Color', 'w');

x_rail = linspace(-1.2*abs(x0), 1.2*abs(x0), 300);
z_rail = a * x_rail.^2;

% Dibujar la vía parabólica
plot(x_rail, z_rail, 'k--', 'LineWidth', 2); hold on;
grid on; axis equal;
xlabel('x (m)'); ylabel('z (m)');
title('Animación en Tiempo Real con Vector de Fuerza Normal \vec{N}');

% Objetos gráficos dinámicos
h_vehiculo = plot(x(1), z(1), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
h_texto = text(x(1), z(1) + 1.2, sprintf('t = %.2f s | G = %.2f g', t(1), fuerza_g(1)), ...
               'FontSize', 10, 'FontWeight', 'bold');

% Vector de Fuerza Normal N_vector
escala_N = 0.001; 
h_quiver = quiver(x(1), z(1), Nx(1)*escala_N, Nz(1)*escala_N, 0, ...
                  'g', 'LineWidth', 2.5, 'MaxHeadSize', 0.8);

legend([h_vehiculo, h_quiver], {'Vehículo', 'Fuerza Normal \vec{N}'}, 'Location', 'northeast');

%% --- CONFIGURACIÓN DEL GRABADOR DE VÍDEO (VideoWriter) ---
nombre_video = 'animacion_parabola.mp4';
v = VideoWriter(nombre_video, 'MPEG-4');
v.FrameRate = 30; % Fotogramas por segundo (30 fps)
open(v);          % Abrir el archivo de vídeo para empezar a escribir

% Bucle de reproducción y captura de fotogramas
for k = 1:length(t)
    % Actualizar posición del vehículo, texto y vector N
    set(h_vehiculo, 'XData', x(k), 'YData', z(k));
    set(h_texto, 'Position', [x(k), z(k) + 1.2, 0], ...
                 'String', sprintf('t = %.2f s | G = %.2f g', t(k), fuerza_g(k)));
    set(h_quiver, 'XData', x(k), 'YData', z(k), ...
                  'UData', Nx(k)*escala_N, 'VData', Nz(k)*escala_N);
    
    drawnow;
    
    % CAPTURAR EL FOTOGRAMA ACTUAL Y ESCRIBIRLO EN EL VÍDEO
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v); % Cerrar el archivo de vídeo al terminar
disp(['¡Vídeo guardado con éxito como: ', nombre_video, '!']);


%% 7. FUNCIÓN INTERNA DE ECUACIONES DIFERENCIALES
function dydt = ec_movimiento_b(~, y, g, a)
    % Vector de estado: y(1) = x, y(2) = v_x
    x = y(1);
    vx = y(2);
    
    % Aceleración d(v_x)/dt
    ax = (-4 * a^2 * x * vx^2 - 2 * g * a * x) / (1 + 4 * a^2 * x^2);
    
    dydt = [vx; ax];
end