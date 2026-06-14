% Open-loop experiment for data collection in TCLab
%
% Initializes TCLab, applies a sequence of open-loop controls and records
% the corresponding temperature.
%
% If you see the warning 'Computation time exceeded sampling time by x
% seconds at sample k', it is because the computation time in a given
% loop took more than the sampling period Ts. Try to disable the rt_plot
% flag to fix it or increase Ts.
%
% Functions called: tclab.
%
% J. Miranda Lemos and Afonso Botelho, IST, May 2023
%__________________________________________________________________________

% Initialization
clear all
close all
clc
tclab;

% Experiment parameters
T = 8000; % experiment duration [s]
Ts = 5; % sampling period [s]3000
N = T/Ts; % number of samples to collect

% Open-loop profile
% Test 1; Q3.1: T=3000 (lab)
% u = zeros(2,N);
% u(1,2:120) = 30;
% u(1,121:240) = 33;
% u(1,241:360) = 30;
% u(1,361:480) = 27;
% u(1,481:600) = 35;
%u(2,75:104) = 80;

% Test 2; Q3.1: T = 4200 (home)
% u = zeros(2,N);
% u(1,2:120) = 30;
% u(1,121:240) = 33;
% u(1,241:360) = 30;
% u(1,361:480) = 27;
% u(1,481:600) = 35;
% u(1,601:720) = 25;
% u(1,721:840) = 31;

% Test 3; Q3.1: T = 8000 (home)
% u = zeros(2,N);
% u(1,2:200) = 30;
% u(1,201:400) = 33;
% u(1,401:600) = 30;
% u(1,601:800) = 27;
% u(1,801:1000) = 35;
% u(1,1001:1200) = 25;
% u(1,1201:1400) = 31;
% u(1,1401:1600) = 28;

% Test 3; Q3.1: T = 8000 (home)
u = zeros(2,N);
u(1,2:200) = 30;
u(1,201:400) = 30;
u(1,401:600) = 33;
u(1,601:800) = 30;
u(1,801:1000) = 27;
u(1,1001:1200) = 35;
u(1,1201:1400) = 23;
u(1,1401:1600) = 28;

%Test 1 Q3.2: < 1000 seconds (lab)
% u = zeros(2,N);
% u(1,2:30) = 50;
% u(1,31:50)= 25;
% u(1,51:60)= 0;
% u(1,61:90)= 35;
% u(1,91:100) = 70;
% u(1,101:120)= 45;
% u(1,121:150)= 20;
% u(1,151:170)= 60;
% u(1,171:190)= 0;
% u(1,191:200)= 80;

%Test 2 Q3.2: T = 950 seconds (home)
% u = zeros(2,N);
% u(1,2:20) = 45;
% u(1,21:30)= 10;
% u(1,31:45)= 30;
% u(1,46:55)= 70;
% u(1,56:65) = 0;
% u(1,66:80)= 60;
% u(1,81:90)= 25;
% u(1,91:105)= 50;
% u(1,106:115)= 15;
% u(1,116:130)= 35;
% u(1,131:140)= 0;
% u(1,141:155)= 20;
% u(1,156:165)= 0;
% u(1,166:180)= 80;








% Real-time plot flag. If true, plots the input and measured temperature in
% real time. If false, only plots at the end of the experiment and instead
% prints the results in the command window.
rt_plot = true;

% Initialize figure and signals
if rt_plot
    figure
    drawnow;
end
t = nan(1,N);
y = nan(2,N);

% String with date for saving results
timestr = char(datetime('now','Format','yyMMdd_HHmmSS'));

% Signals the start of the experiment by lighting the LED
led(1)
disp('Temperature test started.')

for k=1:N
    tic;

    % Computes analog time
    t(k) = (k-1)*Ts;

    % Reads the sensor temperatures
    y(1,k) = T1C();
    y(2,k) = T2C();

    % Applies the control variables to the plant
    h1(u(1,k));
    h2(u(2,k));

    if rt_plot
        % Plots results
        clf
        subplot(2,1,1), hold on, grid on   
        plot(t(1:k),y(1,1:k),'.','MarkerSize',10)
        plot(t(1:k),y(2,1:k),'.','MarkerSize',10)
        legend('Temperature 1','Temperature 2','Location','northwest')
        xlabel('Time [s]')
        ylabel('Temperature [°C]')
        subplot(2,1,2), hold on, grid on   
        stairs(t(1:k),u(1,1:k),'LineWidth',2)
        stairs(t(1:k),u(2,1:k),'LineWidth',2)
        legend('Heater 1','Heater 2','Location','northwest')
        xlabel('Time [s]')
        ylabel('Heater [%]')
        ylim([0 100]);
        drawnow;
    else
        fprintf('t = %d, y1 = %.1f C, y2 = %.1f C, u1 = %.1f, u2 = %.1f\n',t(k),y(1,k),y(2,k),u(1,k),u(2,k))
    end

    % Check if computation time did not exceed sampling time
    if toc > Ts
        warning('Computation time exceeded sampling time by %.2f s at sample %d.',toc-Ts,k)
    end
    % Waits for the begining of the new sampling interval
    pause(max(0,Ts-toc));
end

% Turns off both heaters at the end of the experiment
h1(0);
h2(0);

% Signals the end of the experiment by shutting off the LED
led(0)

disp('Temperature test complete.')

if ~rt_plot
    figure
    subplot(2,1,1), hold on, grid on   
    plot(t,y(1,:),'.','MarkerSize',10)
    plot(t,y(2,:),'.','MarkerSize',10)
    legend('Temperature 1','Temperature 2','Location','best')
    xlabel('Time [s]')
    ylabel('Temperature [°C]')
    subplot(2,1,2), hold on, grid on   
    stairs(t,u(1,:),'LineWidth',2)
    stairs(t,u(2,:),'LineWidth',2)
    legend('Heater 1','Heater 2','Location','best')
    xlabel('Time [s]')
    ylabel('Heater control [%]')
    ylim([0 100]);
end

%--------------------------------------------------------------------------

% Save figure and experiment data to file
exportgraphics(gcf,['openloop_plot_',timestr,'.png'],'Resolution',300)
save(['openloop_data_',timestr,'.mat'],'y','u','t');

%--------------------------------------------------------------------------
% End of File


