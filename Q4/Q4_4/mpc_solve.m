% Question 4.4
function u0 = mpc_solve(x0,H,R,A,B,C, u_ss)
    % K_RH computation : saved in the K_RH_H vector to plot
    W = build_w(A,B,C,H);
    M = (W')*W + R*eye(length(W));
    pi = build_pi(A,C,H);
    F = 2*M;
    f = 2*(x0')*(pi')*W;
    f = f';
    lb = ones(H, 1) * (0-u_ss); % u >= 0%
    ub = ones(H, 1)* (100-u_ss); % u <= 100%
    U=quadprog(F,f, [], [], [], [], lb, ub);
    u0 = U(1);
end