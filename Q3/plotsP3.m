% 1. Load the .mat file into a structure (replace 'data.mat' with your actual file name)
data = load('openloop_data_2.mat');

% 2. Extract the variables from the loaded structure
t = data.t;
y = data.y;
u = data.u;

% 3. Plot the data
figure

subplot(2,1,1), hold on, grid on   
plot(t,y(1,:),'.','MarkerSize',10)
plot(t,y(2,:),'.','MarkerSize',10)
xlabel('Time [s]', 'FontSize', 18)          
ylabel('Temperature [°C]', 'FontSize', 18)  
lgd1 = legend('Temperature 1','Temperature 2','Location','best');
set(lgd1, 'FontSize', 18);                  
set(gca, 'FontSize', 18);                   

subplot(2,1,2), hold on, grid on   
stairs(t,u(1,:),'LineWidth',2)
stairs(t,u(2,:),'LineWidth',2)
xlabel('Time [s]', 'FontSize', 18)         
ylabel('Heater control [%]', 'FontSize', 18) 
ylim([0 100]);
lgd2 = legend('Heater 1','Heater 2','Location','best');
set(lgd2, 'FontSize', 18);                
set(gca, 'FontSize', 18);                   