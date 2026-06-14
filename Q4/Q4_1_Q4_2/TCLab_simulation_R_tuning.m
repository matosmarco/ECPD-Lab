% Simulation of the TCLab linear model previously identified
% Afonso Botelho and J. Miranda Lemos, IST, May 2023
% Adapted for Multiple Horizon Testing
%__________________________________________________________________________

% Initialization
clear all
close all
clc

% 
R_values = [0.001, 0.1, 1,  5, 10, 50, 100]; % Valores de R em escala logarítmica
H_fixed = 50;                   % Horizonte fixado da alínea anterior
N_sim = 800;               % Número de iterações (800 passos * 5s = 4000 segundos)

% Paleta de cores para distinguir as várias linhas (A LINHA QUE FALTAVA!)
colors = lines(length(R_values)); 

% 1. Criar as Figuras
fig_abs = figure('Units','normalized','Position',[0.2 0.5 0.3 0.4]);
fig_inc = figure('Units','normalized','Position',[0.5 0.5 0.3 0.4]);

% Variáveis para guardar os gráficos e facilitar a legenda depois
p_y_abs = gobjects(1, length(R_values));
p_u_abs = gobjects(1, length(R_values));
p_y_inc = gobjects(1, length(R_values));
p_u_inc = gobjects(1, length(R_values));

% Load model
load('singleheater_model.mat','A','B','C','Ke','e_var','y_ss','u_ss','Ts');
n = size(A,1);
e_std = 0; % MUDAR PARA O VALOR ANTERIOR DEPOIS (Sem ruído para este teste)

% Build the functions
x_ss = [eye(n)-A; C]\[B*u_ss; y_ss];
c1 = ((eye(n)-A)*x_ss - B*u_ss);
c2 = (y_ss - C*x_ss);
h1 = @(x,u) A*x + B*u + Ke*e_std*randn + c1; % apply control
T1C = @(x) C*x + e_std*randn + c2; % read temperature

% Initial conditions (start at ambient temperature)
Dx0Dy0 = [eye(n)-A, zeros(n,1); C, -1]\[-B*u_ss; 0];
Dx0 = Dx0Dy0(1:n);

fprintf('Running simulation...\n')
for i = 1:length(R_values)
    H = H_fixed;
    R = R_values(i);
    
    % INICIALIZAR VARIÁVEIS DENTRO DO CICLO 
    % (Crucial para a simulação reiniciar limpa a cada valor de H)
    t = nan(1,N_sim);
    x = nan(n,N_sim);
    y = nan(1,N_sim);
    Dy = nan(1,N_sim);
    Du = nan(1,N_sim);
    Dx = nan(n,N_sim);
    u = nan(1,N_sim);
    
    % Estado Inicial
    x(:,1) = Dx0 + x_ss;
    
    for k = 1:N_sim-1
        % Computes analog time
        t(k) = (k-1)*Ts;
    
        % Reads the sensor temperature
        y(:,k) = T1C(x(:,k));
    
        % Compute incremental variables
        Dy(:,k) = y(:,k) - y_ss;
        Dx(:,k) = x(:,k) - x_ss;
        
        % Chamar função do MPC NÃO RESTRINGIDO (Questão 4.2)
        Du(:,k) = mpc_solve(Dx(:,k), H, R, A, B, C); 
        
        u(:,k) = u_ss + Du(:,k);
        % Applies the control variable to the plant
        x(:,k+1) = h1(x(:,k),u(:,k));
    end

    % Preencher os últimos pontos para fechar as linhas no gráfico
    t(N_sim) = (N_sim-1)*Ts;
    y(:,N_sim) = T1C(x(:,N_sim));
    Dy(:,N_sim) = y(:,N_sim) - y_ss;
    Du(:,N_sim) = Du(:,N_sim-1);
    u(:,N_sim) = u(:,N_sim-1);
    
    % String para a legenda
    leg_str = sprintf('R = %.3f', R);
    
    % --- DESENHAR NA FIGURA 1: VARIÁVEIS ABSOLUTAS ---
    figure(fig_abs);
    
    subplot(2,1,1); hold on; grid on;
    p_y_abs(i) = plot(t, y, '.-', 'MarkerSize', 5, 'Color', colors(i,:), 'DisplayName', leg_str);
    
    subplot(2,1,2); hold on; grid on;
    p_u_abs(i) = stairs(t, u, 'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', leg_str);
    
    % --- DESENHAR NA FIGURA 2: VARIÁVEIS INCREMENTAIS ---
    figure(fig_inc);
    
    subplot(2,1,1); hold on; grid on;
    p_y_inc(i) = plot(t, Dy, '.-', 'MarkerSize', 5, 'Color', colors(i,:), 'DisplayName', leg_str);
    
    subplot(2,1,2); hold on; grid on;
    p_u_inc(i) = stairs(t, Du, 'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', leg_str);
end
fprintf(' Done.\n');

% Figura 1: Absolutas
figure(fig_abs);
subplot(2,1,1);
title('Absolute input/output');
xlabel('Time [s]');
ylabel('y [°C]');
yl_y = yline(y_ss, 'k--', 'DisplayName', '$\bar{y}$');
legend([p_y_abs, yl_y], 'Interpreter', 'latex', 'Location', 'best');

subplot(2,1,2);
xlabel('Time [s]');
ylabel('u [%]');
yline(0, 'r--', 'HandleVisibility', 'off');   
yline(100, 'r--', 'HandleVisibility', 'off'); 
yl_u = yline(u_ss, 'k--', 'DisplayName', '$\bar{u}$');
legend([p_u_abs, yl_u], 'Interpreter', 'latex', 'Location', 'best');

% Figura 2: Incrementais
figure(fig_inc);
subplot(2,1,1);
title('Incremental input/output');
xlabel('Time [s]');
ylabel('\Delta{y} [°C]');
legend(p_y_inc, 'Location', 'best');

subplot(2,1,2);
xlabel('Time [s]');
ylabel('\Delta{u} [%]');
yl_du_min = yline(-u_ss, 'r--', 'DisplayName', 'Limit (-$\bar{u}$)');
yl_du_max = yline(100-u_ss, 'r--', 'DisplayName', 'Limit (100-$\bar{u}$)');
legend([p_u_inc, yl_du_min, yl_du_max], 'Interpreter', 'latex', 'Location', 'best');