% Question 5
function [u0, N] = mpc_solve(x0, H, R, A, B, C, u_bar, y_ss, ref)
% x0: Estado estimado incremental dX
% u_bar: Controlo feedforward absoluto (u_ss + Du_bar)
% ref: Referência absoluta r(k)

N = 0; % counter of MPC exitflag
alpha_weight = 1e3; % Peso da soft constraint (podes ajustar se necessário)
n = size(A, 1); % State size

W = build_w(A,B,C,H);
M = (W')*W + R*eye(length(W));
pi = build_pi(A,C,H);

% Matrizes de custo (Quadprog)
F = 2 * blkdiag(M, alpha_weight*eye(H));
f = 2 * (x0') * (pi') * W;
f = f';
f = [f; zeros(H,1)];

% Limites de Controlo Incremental: 0 <= du + u_bar <= 100
lb_eta = zeros(H,1);
ub_eta = inf(H,1);
lb = ones(H, 1) * (0 - u_bar);   % u >= 0%
ub = ones(H, 1) * (100 - u_bar); % u <= 100%
lb = [lb; lb_eta];
ub = [ub; ub_eta];

% Limite de Output (Soft Constraint): y_pred <= 55 
% Como y_pred = delta_y + ref, temos delta_y <= 55 - ref
y_max = ones(H,1)*55 - ones(H,1)*ref; 

A_ineq = [W, -eye(H)]; 
b_ineq = y_max - pi*x0;

% Resolver o problema de otimização
[U, ~, exitflag, ~] = quadprog(F, f, A_ineq, b_ineq, [], [], lb, ub);

if exitflag == 1
    u0 = U(1);
    N = 1; 
else
    fprintf("MPC optimization failed.\n");
    u0 = 0;
end

end