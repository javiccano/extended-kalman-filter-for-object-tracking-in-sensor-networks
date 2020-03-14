%%%%%%%%%%%% Authors: Javier Carnerero Cano   %%%%%%%%%%%%
%%%%%%%%%%%%          Vicente Gallardo Cabrera %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc



x_0_real = [-0.5; 0]; % posicion inicial
u = [0.02; 0.01]; % vector de desplazamiento
 
% posicion sensores
s_1 = [0; 0];
s_2 = [0; 1];
s_3 = [1; 0];
s_4 = [1; 1];
s = [s_1 s_2 s_3 s_4];
sigma_g = 0.02; % varianza del ruido, se puede variar
instantes = 100; % numero de mediciones

x_t_real = zeros(2, instantes);
A = eye(size(x_t_real, 1)); % matriz de transicion de la posicion

y_t = zeros(size(s, 2), instantes);

for t = 1 : instantes
    if t == 1
        x_t_real(:, t) = x_0_real; 
    else
        x_t_real(:, t) = A * x_t_real(:, t - 1) + u;    
    end
    g = sigma_g * randn(size(s, 2), 1);    
    for i = 1 : size(y_t, 1)
        y_t(i, t) = exp(-0.5 * norm(x_t_real(:, t) - s(:, i))^2) + g(i);
    end
end

y_t_1 = y_t(1, :);
y_t_2 = y_t(2, :);
y_t_3 = y_t(3, :);
y_t_4 = y_t(4, :);

figure, plot(s_1(1), s_1(2), 'o', 'color', 'r', 'LineWidth', 2);
hold on, plot(s_2(1), s_2(2), 'o', 'color', 'g', 'LineWidth', 2);
hold on, plot(s_3(1), s_3(2), 'o', 'color', 'y', 'LineWidth', 2);
hold on, plot(s_4(1), s_4(2), 'o', 'color', 'c', 'LineWidth', 2);
hold on, plot(x_t_real(1, :), x_t_real(2, :), '.', 'color', 'b', 'LineWidth', 2);
xlabel('{\itx}_{1}'), ylabel('{\itx}_{2}');
legend('sensor 1', 'sensor 2', 'sensor 3', 'sensor 4', 'posición real del objeto', 'Location', 'best');
title('Plano del movimiento del objeto');
set(gcf, 'color', 'w'), grid;

figure, plot(y_t_1, 'color', 'r', 'LineWidth', 2);
hold on, plot(y_t_2, 'color', 'g', 'LineWidth', 2);
hold on, plot(y_t_3, 'color', 'y', 'LineWidth', 2);
hold on, plot(y_t_4, 'color', 'c', 'LineWidth', 2);
xlabel('instante de muestreo'), ylabel('{\ity}');
title('Valor de las observaciones medidas');
legend('sensor 1', 'sensor 2', 'sensor 3', 'sensor 4', 'Location', 'best');
set(gcf, 'color', 'w'), grid;








% Filtro de Kalman %


mi_0 = [-0.5; 0];
sigma_0 = sqrt(0.04);

x_0 = mi_0 + sigma_0 * randn(2, 1); % posicion inicial estimada
% x_0 = [-0.446070271656680; 0.098857411075882]; % posicion inicial estimada empleada en la memoria

x_t = zeros(2, instantes);
x_t(:, 1) = x_0; 
P_t = sigma_0^2 * eye(size(x_t, 1)); % covarianza de la distribucion de la posicion inicial
Q = 0; % matriz de covarianzas de u (determinista)
R = sigma_g^2 * eye(size(s, 2)); % matriz de covarianzas de g
f = zeros(size(s, 2), 1);

for t = 2 : instantes 
    % prediccion
    x_t_t_1 = A * x_t(:, t - 1) + u; 
    P_t_t_1 = Q + A * P_t * A';
    
    % actualizacion
    J = jacob(x_t_t_1, s);
    K_t = P_t_t_1 * J'/(J * P_t_t_1 * J' + R);    
    for i = 1 : length(f)
        f(i) = exp(-0.5 * norm(x_t_t_1 - s(:, i))^2);
    end

    x_t(:, t) = x_t_t_1 + K_t * (y_t(:, t) - f);
    P_t = P_t_t_1 - K_t * (J * P_t_t_1 * J') * K_t';
end

error_1_2_3_4 = 2 * immse(x_t, x_t_real);








% Todos los sensores

figure, plot(s_1(1), s_1(2), 'o', 'color', 'r', 'LineWidth', 2);
hold on, plot(s_2(1), s_2(2), 'o', 'color', 'g', 'LineWidth', 2);
hold on, plot(s_3(1), s_3(2), 'o', 'color', 'y', 'LineWidth', 2);
hold on, plot(s_4(1), s_4(2), 'o', 'color', 'c', 'LineWidth', 2);
hold on, plot(x_t_real(1, :), x_t_real(2, :), '.', 'color', 'b', 'LineWidth', 2);
hold on, plot(x_t(1, :), x_t(2, :), '.','color' , 'm', 'LineWidth', 2);
xlabel('{\itx}_{1}'), ylabel('{\itx}_{2}');
legend('sensor 1', 'sensor 2', 'sensor 3', 'sensor 4', 'posición real del objeto', 'posición media estimada del objeto', 'Location', 'best');
title('Plano del movimiento del objeto');
set(gcf, 'color', 'w'), grid;




 
% solo el sensor 1

s = s_1;
y_t = y_t_1;
R = sigma_g^2 * eye(size(s, 2));
f = zeros(size(s, 2), 1);

for t = 2 : instantes 
    x_t_t_1 = A * x_t(:, t - 1) + u; 
    P_t_t_1 = Q + A * P_t * A';
    
    J = jacob(x_t_t_1, s);
    K_t = P_t_t_1 * J'/(J * P_t_t_1 * J' + R);    
    for i = 1 : length(f)
        f(i) = exp(-0.5 * norm(x_t_t_1 - s(:, i))^2);
    end

    x_t(:, t) = x_t_t_1 + K_t * (y_t(:, t) - f);
    P_t = P_t_t_1 - K_t * (J * P_t_t_1 * J') * K_t';
end

error_1 = 2 * immse(x_t, x_t_real);

figure, plot(s_1(1), s_1(2), 'o', 'color', 'r', 'LineWidth', 2);
hold on, plot(x_t_real(1, :), x_t_real(2, :), '.', 'color', 'b', 'LineWidth', 2);
hold on, plot(x_t(1, :), x_t(2, :), '.','color' , 'm', 'LineWidth', 2);
xlabel('{\itx}_{1}'), ylabel('{\itx}_{2}');
legend('sensor 1', 'posición real del objeto', 'posición media estimada del objeto', 'Location', 'best');
title('Plano del movimiento del objeto');
set(gcf, 'color', 'w'), grid;




 
% sensor 1 y sensor 2

s = [s_1 s_2];
y_t = [y_t_1; y_t_2];
R = sigma_g^2 * eye(size(s, 2));
f = zeros(size(s, 2), 1);

for t = 2 : instantes 
    x_t_t_1 = A * x_t(:, t - 1) + u; 
    P_t_t_1 = Q + A * P_t * A';
    
    J = jacob(x_t_t_1, s);
    K_t = P_t_t_1 * J'/(J * P_t_t_1 * J' + R);    
    for i = 1 : length(f)
        f(i) = exp(-0.5 * norm(x_t_t_1 - s(:, i))^2);
    end

    x_t(:, t) = x_t_t_1 + K_t * (y_t(:, t) - f);
    P_t = P_t_t_1 - K_t * (J * P_t_t_1 * J') * K_t';
end

error_1_2 = 2 * immse(x_t, x_t_real);

figure, plot(s_1(1), s_1(2), 'o', 'color', 'r', 'LineWidth', 2);
hold on, plot(s_2(1), s_2(2), 'o', 'color', 'g', 'LineWidth', 2);
hold on, plot(x_t_real(1, :), x_t_real(2, :), '.', 'color', 'b', 'LineWidth', 2);
hold on, plot(x_t(1, :), x_t(2, :), '.','color' , 'm', 'LineWidth', 2);
xlabel('{\itx}_{1}'), ylabel('{\itx}_{2}');
legend('sensor 1', 'sensor 2', 'posición real del objeto', 'posición media estimada del objeto', 'Location', 'best');
title('Plano del movimiento del objeto');
set(gcf, 'color', 'w'), grid;


%%%%% mas ejemplos %%%%%%

 
% sensor 1, sensor 2 y sensor 3

s = [s_1 s_2 s_3];
y_t = [y_t_1; y_t_2; y_t_3];
R = sigma_g^2 * eye(size(s, 2));
f = zeros(size(s, 2), 1);

for t = 2 : instantes 
    x_t_t_1 = A * x_t(:, t - 1) + u; 
    P_t_t_1 = Q + A * P_t * A';
    
    J = jacob(x_t_t_1, s);
    K_t = P_t_t_1 * J'/(J * P_t_t_1 * J' + R);    
    for i = 1 : length(f)
        f(i) = exp(-0.5 * norm(x_t_t_1 - s(:, i))^2);
    end

    x_t(:, t) = x_t_t_1 + K_t * (y_t(:, t) - f);
    P_t = P_t_t_1 - K_t * (J * P_t_t_1 * J') * K_t';
end

error_1_2_3 = 2 * immse(x_t, x_t_real);

figure, plot(s_1(1), s_1(2), 'o', 'color', 'r', 'LineWidth', 2);
hold on, plot(s_2(1), s_2(2), 'o', 'color', 'g', 'LineWidth', 2);
hold on, plot(s_3(1), s_3(2), 'o', 'color', 'y', 'LineWidth', 2);
hold on, plot(x_t_real(1, :), x_t_real(2, :), '.', 'color', 'b', 'LineWidth', 2);
hold on, plot(x_t(1, :), x_t(2, :), '.','color' , 'm', 'LineWidth', 2);
xlabel('{\itx}_{1}'), ylabel('{\itx}_{2}');
legend('sensor 1', 'sensor 2', 'sensor 3', 'posición real del objeto', 'posición media estimada del objeto', 'Location', 'best');
title('Plano del movimiento del objeto');
set(gcf, 'color', 'w'), grid;




 
% sensor 1, sensor 3 y sensor 4

s = [s_1 s_3 s_4];
y_t = [y_t_1; y_t_3; y_t_4];
R = sigma_g^2 * eye(size(s, 2));
f = zeros(size(s, 2), 1);

for t = 2 : instantes 
    x_t_t_1 = A * x_t(:, t - 1) + u; 
    P_t_t_1 = Q + A * P_t * A';
    
    J = jacob(x_t_t_1, s);
    K_t = P_t_t_1 * J'/(J * P_t_t_1 * J' + R);    
    for i = 1 : length(f)
        f(i) = exp(-0.5 * norm(x_t_t_1 - s(:, i))^2);
    end

    x_t(:, t) = x_t_t_1 + K_t * (y_t(:, t) - f);
    P_t = P_t_t_1 - K_t * (J * P_t_t_1 * J') * K_t';
end

error_1_3_4 = 2 * immse(x_t, x_t_real);

figure, plot(s_1(1), s_1(2), 'o', 'color', 'r', 'LineWidth', 2);
hold on, plot(s_3(1), s_3(2), 'o', 'color', 'y', 'LineWidth', 2);
hold on, plot(s_4(1), s_4(2), 'o', 'color', 'c', 'LineWidth', 2);
hold on, plot(x_t_real(1, :), x_t_real(2, :), '.', 'color', 'b', 'LineWidth', 2);
hold on, plot(x_t(1, :), x_t(2, :), '.','color' , 'm', 'LineWidth', 2);
xlabel('{\itx}_{1}'), ylabel('{\itx}_{2}');
legend('sensor 1', 'sensor 3', 'sensor 4', 'posición real del objeto', 'posición media estimada del objeto', 'Location', 'best');
title('Plano del movimiento del objeto');
set(gcf, 'color', 'w'), grid;




 
% solo el sensor 4

s = s_4;
y_t = y_t_4;
R = sigma_g^2 * eye(size(s, 2));
f = zeros(size(s, 2), 1);

for t = 2 : instantes 
    x_t_t_1 = A * x_t(:, t - 1) +  u; 
    P_t_t_1 = Q + A * P_t * A';
    
    J = jacob(x_t_t_1, s);
    K_t = P_t_t_1 * J'/(J * P_t_t_1 * J' + R);    
    for i = 1 : length(f)
        f(i) = exp(-0.5 * norm(x_t_t_1 - s(:, i))^2);
    end

    x_t(:, t) = x_t_t_1 + K_t * (y_t(:, t) - f);
    P_t = P_t_t_1 - K_t * (J * P_t_t_1 * J') * K_t';
end

error_4 = 2 * immse(x_t, x_t_real);


figure, plot(s_4(1), s_4(2), 'o', 'color', 'c', 'LineWidth', 2);
hold on, plot(x_t_real(1, :), x_t_real(2, :), '.', 'color', 'b', 'LineWidth', 2);
hold on, plot(x_t(1, :), x_t(2, :), '.','color' , 'm', 'LineWidth', 2);
xlabel('{\itx}_{1}'), ylabel('{\itx}_{2}');
legend('sensor 4', 'posición real del objeto', 'posición media estimada del objeto', 'Location', 'best');
title('Plano del movimiento del objeto');
set(gcf, 'color', 'w'), grid;


