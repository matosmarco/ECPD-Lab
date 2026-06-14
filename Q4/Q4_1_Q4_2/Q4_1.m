A = 1.2;
B = 1;
C = 1;
Q= 1;
H = 5;
x0 = 1;
R = 10;

W = build_w(A,B,C,H);
M = (W')*W + R*eye(length(W));
pi = build_pi(A,C,H);
e1 = [1, zeros(1, H-1)];

K_RH_H = (e1 / M) * W' * pi;

u0_q2 = -K_RH_H*x0;


u0 = mpc_solve(x0,H,R,A,B,C);

fprintf('%.2f vs. %.2f', u0_q2, u0)