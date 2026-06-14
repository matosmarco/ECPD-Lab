
% Question 4.5
% function [u0, N] = mpc_solve(x0,H,R,A,B,C, u_ss, y_ss, Dr)
% % K_RH computation : saved in the K_RH_H vector to plot
% N = 0; % counter of MPC exitflag
% n = size(A, 1); % State size
% W = build_w(A,B,C,H);
% M = (W')*W + R*eye(length(W));
% pi = build_pi(A,C,H);
% F = 2*M;
% f = 2*(x0')*(pi')*W;
% f = f';
% lb = ones(H, 1) * (0-u_ss); % u >= 0%
% ub = ones(H, 1)* (100-u_ss); % u <= 100%
% %U=quadprog(F,f, [], [], [], [], lb, ub);
% y_max = ones(H,1)*55 - ones(H,1)*y_ss- ones(H,1)*Dr; % Limite incremental
% %y_max = ones(H,1)*55;
% %G = eye(n);
% A_ineq = W; % G = 1 (G*W)
% b_ineq = y_max- pi*x0;
% [U,~ , exitflag,~]= quadprog(F,f,A_ineq, b_ineq, [], [], lb, ub);
% 
% if exitflag == 1
%     fprintf("MPC optimization was successful.");
%     u0 = U(1);
%     N = N+1; 
% else
%     fprintf("MPC optimization failed.");
%     u0 = 0;
% end
% 
% end

% % With soft constraints
function [u0, N] = mpc_solve(x0,H,R,A,B,C, u_ss, y_ss, Dr)
% K_RH computation : saved in the K_RH_H vector to plot
N = 0; % counter of MPC exitflag
alpha_weight = 1e3; % TODO: verificar valor
n = size(A, 1); % State size
W = build_w(A,B,C,H);
M = (W')*W + R*eye(length(W));
pi = build_pi(A,C,H);
%F = 2*M;
F = 2 *blkdiag(M,alpha_weight*eye(H));
f = 2*(x0')*(pi')*W;
f = f';
f = [f; zeros(H,1)];
lb_eta = zeros(H,1);
ub_eta = inf(H,1);
lb = ones(H, 1) * (0-u_ss); % u >= 0%
ub = ones(H, 1)* (100-u_ss); % u <= 100%
lb = [lb; lb_eta];
ub = [ub; ub_eta];
%U=quadprog(F,f, [], [], [], [], lb, ub);

y_max = ones(H,1)*55 - ones(H,1)*y_ss- ones(H,1)*Dr; % Limite incremental
%y_max = ones(H,1)*55;
%G = eye(n);
A_ineq = [W, -eye(H)]; % G = 1 (G*W)
b_ineq = y_max- pi*x0;
[U,~ , exitflag,~]= quadprog(F,f,A_ineq, b_ineq, [], [], lb, ub);

if exitflag == 1
    fprintf("MPC optimization was successful.");
    u0 = U(1);
    N = N+1; 
else
    fprintf("MPC optimization failed.");
    u0 = 0;
end

end
