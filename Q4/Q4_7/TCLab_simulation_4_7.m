% Simulation of the TCLab linear model previously identified
%
% Loads the model identified in the TCLab_identification script, creates
% the h1 and T1C functions that mimick the TCLab interface, and performs a
% simulation starting at ambient temperature.
% You will be developing and testing your MPC controller and Kalman filter
% in this simulation environment. 
%
% Afonso Botelho and J. Miranda Lemos, IST, May 2023
%__________________________________________________________________________

% Initialization
clear all
close all
clc

H = 50;
R = 0.1;

% Load model
load('singleheater_model.mat','A','B','C','Ke','e_var','y_ss','u_ss','Ts');
n = size(A,1);
e_std = sqrt(e_var); % input disturbance standard deviation
%e_std = 0; % MUDAR PARA O VALOR ANTERIOR; SEGUIMOS A SUGESTÃO DE P4.0


% Build the functions for applying the control and reading the temperature,
% mimicking the TCLab interface
x_ss = [eye(n)-A; C]\[B*u_ss; y_ss];
c1 = ((eye(n)-A)*x_ss - B*u_ss);
c1 = c1 * 1.1; % ver com aumento de 10% e sem
c2 = (y_ss - C*x_ss);
h1 = @(x,u) A*x + B*u + Ke*e_std*randn + c1; % apply control
T1C = @(x) C*x + e_std*randn + c2; % read temperature

% Simulation parameters
T = 4000; % Experiment duration [s]
N = T/Ts; % Number of samples to collect

% Open-loop control profile
u = zeros(1,N);


% Perfil de referência absoluta
r = zeros(1,N);
r(1:200)   = 50;
r(201:400) = 40;
r(401:600) = 60;  % O MPC não deve conseguir atingir devido à restrição de segurança!
r(601:end) = 45;

% Initial conditions (start at ambient temperature, i.e. equilibrium for u = 0)
Dx0Dy0 = [eye(n)-A, zeros(n,1); C, -1]\[-B*u_ss; 0];
Dx0 = Dx0Dy0(1:n);

% Initialize signals
t = nan(1,N);
x = nan(n,N);
y = nan(1,N);
Dy = nan(1,N);
Du = nan(1,N);
Dx = nan(n,N+1);
x(:,1) = Dx0 + x_ss;

% Variables added in this question
Dr = r - y_ss;
%fprintf("Value of Dr is: %.2f", Dr)

A_tff = [(eye(n)-A),-B; C, 0]; %matrix obtained in the system on the tracking feedforward 


% Exitflag variables
MPC_success = nan(1,N);


% Kalman filter model
Ad = [A, B; zeros(1,n), 1];
Bd = [B;0];
Cd = [C, 0];
Qe = Ke * e_std * (Ke * e_std)';
Re = e_var;
de = 0.5*e_std;
QEd = [Qe zeros(n,1); zeros(1,n) de];
L = dlqe(Ad, eye(n+1), Cd, QEd, Re);

% ---- NOVA INICIALIZAÇÃO DO FILTRO ----
% A variável passa a chamar-se xdhat_pred (é a predição para k=1)
% E introduzimos aqui o erro exigido na alínea 6
xdhat_pred = zeros(n+1,1);
state_error = C\5;
xdhat_pred(1:n) = Dx0 + state_error;

% variables to save estimates
Dxhat_hist = nan(n,N);
dhat_hist = nan(1,N);
yhat_hist = nan(1,N);

fprintf('Running simulation...\n');
for k = 1:N
    % Computes analog time
    t(k) = (k-1)*Ts;

    % Reads the sensor temperature
    y(:,k) = T1C(x(:,k));

    % Compute incremental variables
    Dy(:,k) = y(:,k) - y_ss;
    Dx(:,k) = x(:,k) - x_ss;
    %Du(:,k) = u(:,k) - u_ss;
    % Kalman Filter
    % Usa a predição atual e corrige-a com a medição y(:,k)
    xdhat = xdhat_pred + L*(Dy(:,k) - Cd*xdhat_pred);

    % Extrair estimativas corrigidas
    Dxhat = xdhat(1:n);
    dhat = xdhat(end);
    Dyhat = C*Dxhat;
    yhat = y_ss + Dyhat;

    %  Feedforward computation
    % Usa a estimativa de perturbação (dhat) acabada de corrigir
    rhs_ss = [B * dhat; Dr(k)];
    ss_vals = A_tff \ rhs_ss;
    Dx_bar = ss_vals(1:n);
    Du_bar = ss_vals(n+1);

    % MPC Control
    % Mudança de variáveis usando o estado estimado (Dxhat)
    dx_hat = Dxhat - Dx_bar; 
    
    % Erro real (apenas para verificação, o MPC não vê isto)
    dx(:,k) = Dx(:,k) - Dx_bar; 

    % Cálculo do Controlo (Ajusta os argumentos do mpc_solve se necessário)
    [du, MPC_success(k)] = mpc_solve(dx_hat, H, R, A, B, C, u_ss + Du_bar, y_ss, r(k));

    Du(:,k) = du + Du_bar;
    u(:,k) = u_ss + Du(:,k);

    % Guardar variáveis no histórico para os plots
    Dxhat_hist(:,k) = Dxhat;
    dhat_hist(:,k) = dhat;
    yhat_hist(:,k) = yhat;

    % Aplica o controlo à planta
    x(:,k+1) = h1(x(:,k),u(:,k));
    
    % Kalman filter prediction (Para k+1)
    xdhat_pred = Ad*xdhat + Bd*Du(:,k);
end

fprintf(' Done.\n');

% Check number of instances that the problem had feasible solution of not

counter = 0;
for ii = 1: length(MPC_success)
    if MPC_success(ii) == 1
        counter= counter +1; 
    end
end

fprintf("MPC was sucessful in %d out of %d samples", counter,length(MPC_success))

%% Plots
% Plot absolute variables
figure('Units','normalized','Position',[0.2 0.5 0.3 0.4])
subplot(2,1,1), hold on, grid on   
title('Absolute input/output')
plot(t,y,'.','MarkerSize',5)
stairs(t,r,'r--','LineWidth',1.5)
yl=yline(y_ss,'k--');
xlabel('Time [s]')
ylabel('y [°C]')
legend('Output y', 'Reference r', '$\bar{y}$', 'Interpreter','latex','Location','best')

subplot(2,1,2), hold on, grid on   
stairs(t,u,'LineWidth',2)
yl=yline(u_ss,'k--');
yline(0,'r--')
yline(100,'r--')
xlabel('Time [s]')
ylabel('u [%]')
legend(yl,'$\bar{u}$','Interpreter','latex','Location','best');

% Plot incremental variables
figure('Units','normalized','Position',[0.5 0.5 0.3 0.4])
subplot(2,1,1), hold on, grid on   
title('Incremental input/output')
plot(t,Dy,'.','MarkerSize',5)
xlabel('Time [s]')
ylabel('\Delta{y} [°C]')
subplot(2,1,2), hold on, grid on   
stairs(t,Du,'LineWidth',2)
yline(-u_ss,'r--')
yline(100-u_ss,'r--')
xlabel('Time [s]')
ylabel('\Delta{u} [%]')

%---QE = e_var * eye(n);-----------------------------------------------------------------------
% End of File

%% Kalman filter plots

figure

subplot(2,1,1)
hold on
grid on

plot(t,y,'b')
plot(t,yhat_hist,'r--')

xlabel('Time [s]')
ylabel('Temperature [°C]')

legend('Measured y','Estimated yhat')

subplot(2,1,2)
hold on
grid on

plot(t,dhat_hist,'LineWidth',2)

xlabel('Time [s]')
ylabel('Estimated disturbance dhat')