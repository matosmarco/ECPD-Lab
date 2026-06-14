% Question 2 – Basics on receding horizon control
% Authors:
% - Rafael Galo (103755)
% - Marco Matos (105932)
% - Tomás Modesto (105944)
% - André Salvaterra (106175)

clear;
clc;
% Considering the canonical form
% x(t+1) = A* x(t) + B u(t)
A = 1.2;
B = 1;
C = 1;
Q = 1; % Define 1 as the project statement, but Q = C'*C

R = [0.01, 0.1, 1, 10, 100];

% Q2.1: Compute the optimal LQ state feedback gain.
for ii=1:length(R)
    [KLQ,S,lambda] = dlqr(A,B,Q,R(ii));
    fprintf("LQ state feedback gain (K_LQ) (of R = %.2f):\n", R(ii))
    disp(KLQ)
end

%Q2.2: Compute the optimal receding-horizon gain for different values of the horizon H.
H = 1:30;

% Define an array of distinct colors for each R value
colors = lines(length(R));

% Figure 1: Gains
figure(1); 
clf; % Clears the current figure
hold on 
grid on
grid minor
title('Asymptotic Convergence of Receding-Horizon Gain K_{RH} towards K_{LQ}', 'FontSize', 13)
xlabel('Prediction Horizon, H [discrete steps]', 'FontSize', 12)
ylabel('State Feedback Gain, K', 'FontSize', 12)
ylim([0,1.25])

% Figure 2: Stability
figure(2); 
clf; % Clears the current figure 
hold on
grid on
grid minor
title('Closed-Loop Pole Modulus |\lambda| as a Function of Horizon H', 'FontSize', 18)
xlabel('Prediction Horizon, H [discrete steps]', 'FontSize', 12)
ylabel('Closed-Loop Pole Magnitude, |\lambda|', 'FontSize', 12)
ylim([0,1.25])

% Q2.2 and Q2.3: Compute the optimal receding-horizon gain for different H
for jj=1:length(R)
    K_RH_H = zeros(1, length(H));
    for ii=1:length(H)
        % K_RH computation : saved in the K_RH_H vector to plot
        W = build_w(A,B,C,H(ii));
        M = (W')*W + R(jj)*eye(length(W));
        pi = build_pi(A,C,H(ii));
        e1 = [1, zeros(1, H(ii)-1)];
        K_RH_H(ii) = (e1 / M) * W' * pi;
    end

    [KLQ,~,~] = dlqr(A,B,Q,R(jj));

    % Plot Figure 1: Line + Marker for K_RH
    figure(1)
    plot(H, K_RH_H, '-o', 'Color', colors(jj,:), 'MarkerSize', 4, ...
        'LineWidth', 2, 'MarkerFaceColor', colors(jj,:), ...
        'DisplayName', sprintf('$K_{RH}$ ($R = %g$)', R(jj)));    % Line for KLQ
    yline(KLQ, '--', 'Color', colors(jj,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('$K_{LQ}$ ($R = %g$)', R(jj)));
    % Plot Figure 2: Stability
    figure(2)
    eigenvalues = abs(A - B*K_RH_H);
    plot(H, eigenvalues, '-s', 'Color', colors(jj,:), 'MarkerSize', 4, ...
        'LineWidth', 2.2, 'MarkerFaceColor', colors(jj,:), ...
        'DisplayName', sprintf('$|\\lambda|$ ($R = %g$)', R(jj)));
end

% Legend entry Evaluation of Optimal Control Gain K_RH
figure(1)
lgd1 = legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 18, 'NumColumns', 2); 
set(gca, 'Box', 'on', 'TickDir', 'in', 'LineWidth', 1, 'FontSize', 18)
title('Asymptotic Convergence of Receding-Horizon Gain K_{RH} towards K_{LQ}', 'FontSize', 20)
xlabel('Prediction Horizon, H [discrete steps]', 'FontSize', 20)
ylabel('State Feedback Gain, K', 'FontSize', 20) 

% Legend entry Stability of A - B K_{RH} vs. Horizon H
figure(2)
yline(1, 'r--', 'LineWidth', 2, 'DisplayName', 'Circle Boundary ($|\lambda| = 1$)')
lgd2 = legend('Interpreter','latex','Location', 'best', 'FontSize', 18, 'NumColumns', 3); 
set(gca, 'FontSize', 18)
title('Closed-Loop Pole Modulus |\lambda| as a Function of Horizon H', 'FontSize', 20) 
xlabel('Prediction Horizon, H [discrete steps]', 'FontSize', 20) 
ylabel('Closed-Loop Pole Magnitude, |\lambda|', 'FontSize', 20)