% Simulation of the TCLab linear model previously identified
%
% Loads the model identified in the TCLab_identification script, creates
% the h1 and T1C functions that mimic the TCLab interface, and performs a
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

% Build the functions for applying the control and reading the temperature,
% mimicking the TCLab interface
x_ss = [eye(n)-A; C]\[B*u_ss; y_ss];
c1 = ((eye(n)-A)*x_ss - B*u_ss);
c1 = c1 * 1.1; % Test with and without the 10% increase
c2 = (y_ss - C*x_ss);
h1 = @(x,u) A*x + B*u + Ke*e_std*randn + c1; % apply control
T1C = @(x) C*x + e_std*randn + c2; % read temperature

% Simulation parameters
T = 4000; % Experiment duration [s]
N = T/Ts; % Number of samples to collect

% Open-loop control profile
u = zeros(1,N);
u(:,1:200)   = u_ss;
u(:,201:400) = u_ss+5;
u(:,401:600) = u_ss-5;
u(:,601:800) = u_ss;

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
Dr = (55-y_ss)+10;
%fprintf("Value of Dr is: %.2f", Dr)
A_tff = [(eye(n)-A),-B; C, 0]; % matrix obtained in the system on the tracking feedforward 
rhs_ss = [zeros(n,1); Dr];
ss_vals = A_tff \ rhs_ss;
Dx_bar = ss_vals(1:n);
Du_bar = ss_vals(n+1);

% Exitflag variables
MPC_success = nan(1,N);

% Kalman filter model
Ad = [A, B; zeros(1,n), 1];
Bd = [B;0];
Cd = [C, 0];
Qe = Ke * e_std * (Ke*e_std)';
Re = e_var;
de = 0.1*e_std;
QEd = [Qe zeros(n,1); zeros(1,n) de];
L = dlqe(Ad, eye(n+1), Cd, QEd, Re);

% Initial observer estimate
xdhat = zeros(n+1,1);

% initial estimation error (~5 °C)
state_error = C\5;
xdhat(1:n) = Dx0 + state_error;

% variables to save estimates
Dxhat_hist = nan(n,N);
dhat_hist = nan(1,N);
yhat_hist = nan(1,N);

% Simulate incremental model
fprintf('Running simulation...')
for k = 1:N
    % Computes analog time
    t(k) = (k-1)*Ts;

    % Reads the sensor temperature
    y(:,k) = T1C(x(:,k));

    % Compute incremental variables
    Dy(:,k) = y(:,k) - y_ss;
    Dx(:,k) = x(:,k) - x_ss;
    
    % Kalman filter
    if k == 1
        Du_prev = 0;
    else
        Du_prev = Du(:,k-1);
    end

    % Prediction
    xdhat_pred = Ad*xdhat + Bd*Du_prev;

    % Correction
    xdhat = xdhat_pred + L*(Dy(:,k) - Cd*xdhat_pred);

    % Extract estimates
    Dxhat = xdhat(1:n);
    dhat = xdhat(end);

    % Estimated output
    Dyhat = C*Dxhat;
    yhat = y_ss + Dyhat;

    % Save variables
    Dxhat_hist(:,k) = Dxhat;
    dhat_hist(:,k) = dhat;
    yhat_hist(:,k) = yhat;

    dx(:,k) = Dx(:,k) - Dx_bar;

    % [du, MPC_success(k)] = mpc_solve(dx(:,k), H, R, A, B, C, u_ss + Du_bar, y_ss, Dr);
    % Du(:,k) = du + Du_bar;

    % u(:,k) = u_ss + Du(:,k);
    Du(:,k) = 0;
    u(:,k) = u_ss;
    
    % Applies the control variable to the plant
    x(:,k+1) = h1(x(:,k),u(:,k));
end
fprintf(' Done.\n');

%% Plots
% Plot absolute variables
figure('Units','normalized','Position',[0.2 0.5 0.3 0.4])
subplot(2,1,1), hold on, grid on   
title('Absolute input/output')
plot(t,y,'.','MarkerSize',5)
yl=yline(y_ss,'k--');
xlabel('Time [s]')
ylabel('y [°C]')
legend(yl,'$\bar{y}$','Interpreter','latex','Location','best')
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

%% Kalman filter baseline plots
figure('Name', 'Baseline Kalman Filter')
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

%% Kalman Filter Tuning: Effect of delta_E (P4.6)
fprintf('Testing different values for delta_E...\n');

% 1. Define the delta_E values to test
% Testing different orders of magnitude to clearly observe the difference
%delta_e_vals = [0.001, 0.05, 0.1, 0.5, 1, 5, 10]*e_std; 
delta_e_vals = [0.1, 0.5, 1]*e_std;
colors = ['g', 'm', 'c']; % Colors to distinguish the lines
%colors = ['g', 'm', 'c', 'k', 'r']; % Agora tem 5 cores para 5 valores
%colors = ['g', 'm', 'c', 'k', 'r', 'y', 'b']; % Verde, Magenta, Ciano, Preto, Vermelho, Amarelo, Azul
% 2. Prepare the figures
figure('Name', 'Comparison of \delta_E Effect on Kalman Filter', 'Units','normalized','Position',[0.3 0.3 0.4 0.6])

% Subplot 1: Temperature
ax1 = subplot(2,1,1);
hold on; grid on;
plot(t, y, 'b', 'LineWidth', 1.5, 'DisplayName', 'Measured y (Real + Noise)');
xlabel('Time [s]');
ylabel('Temperature [$^\circ$C]', 'Interpreter', 'latex');
title('Temperature Estimate ($\hat{y}$) for different $\delta_E$ values', 'Interpreter', 'latex');

% Subplot 2: Estimated Disturbance
ax2 = subplot(2,1,2);
hold on; grid on;
xlabel('Time [s]');
ylabel('Estimated $\hat{d}$', 'Interpreter', 'latex');
title('Disturbance Estimate ($\hat{d}$) for different $\delta_E$ values', 'Interpreter', 'latex');


% 3. Loop to test each delta_E
for idx = 1:length(delta_e_vals)
    de_test = delta_e_vals(idx);
    
    % a) Recalculate the QEd matrix and Gain L for this delta_E
    QEd_test = [Qe zeros(n,1); zeros(1,n) de_test];
    L_test = dlqe(Ad, eye(n+1), Cd, QEd_test, Re);
    
    % b) Initialize the estimate (with the ~5 °C error requested in the guide)
    xdhat_test = zeros(n+1,1);
    state_error_test = C\5;
    xdhat_test(1:n) = Dx0 + state_error_test;
    
    % Matrices to save the results of this iteration
    yhat_test_hist = nan(1,N);
    dhat_test_hist = nan(1,N);
    
    % c) Run ONLY the filter equations over the recorded data (Dy, Du)
    for k = 1:N
        if k == 1
            Du_prev_test = 0;
        else
            Du_prev_test = Du(:,k-1); % Same as the main simulation
        end

        % Step 1: Prediction
        xdhat_pred_test = Ad * xdhat_test + Bd * Du_prev_test;

        % Step 2: Correction (using the 'Dy' measurement recorded from the plant)
        xdhat_test = xdhat_pred_test + L_test * (Dy(:,k) - Cd * xdhat_pred_test);

        % Extract estimates
        Dxhat_k = xdhat_test(1:n);
        dhat_k = xdhat_test(end);
        
        yhat_test_hist(k) = y_ss + C * Dxhat_k;
        dhat_test_hist(k) = dhat_k;
    end
    
    % d) Add to the plots
    subplot(2,1,1);
    plot(t, yhat_test_hist, colors(idx), 'LineStyle', '--', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('\\delta_E = %.3f', de_test));
    
    subplot(2,1,2);
    plot(t, dhat_test_hist, colors(idx), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('\\delta_E = %.3f', de_test));
end

% Configure Legends
subplot(2,1,1); legend('Location', 'best');
subplot(2,1,2); legend('Location', 'best');

fprintf('Tuning complete. Please check the generated plots.\n');