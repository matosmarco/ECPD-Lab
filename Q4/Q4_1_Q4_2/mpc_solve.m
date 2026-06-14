% mpc_solve for questions 4.1 and 4.2
function u0 = mpc_solve(x0,H,R,A,B,C)
    % K_RH computation : saved in the K_RH_H vector to plot
    W = build_w(A,B,C,H);
    M = (W')*W + R*eye(length(W));
    pi = build_pi(A,C,H);
    F = 2*M;
    f = 2*(x0')*(pi')*W;
    U=quadprog(F,f);
    u0 = U(1);
end