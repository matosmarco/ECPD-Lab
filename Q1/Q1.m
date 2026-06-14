% Question 1 – Basics on constrained optimization
% Objective: Get acquainted with the basics of numerical 
% function minimization using MATLAB Optimization Toolbox
% Code based on the given file 'ProbBasic.m'
% Authors:
% - Rafael Galo (103755)
% - Marco Matos (105932)
% - Tomás Modesto (105944)
% - André Salvaterra (106175)

clear;
clc;

% Plots function level curves

% Range of independent variables to consider 
x1min=-2;
x1max=2;
x2min=-4;
x2max=4;

% Number of intervals in the mesh grid
N1=100;
N2=100;

xv1 = linspace(x1min,x1max,N1);
xv2 = linspace(x2min,x2max,N2);
[xx1,xx2] = meshgrid(xv1,xv2);

% Computes the function at the different points of the mesh grid
for ii=1:N1
    for jj=1:N2
        x=[xx1(ii,jj); xx2(ii,jj)];
        ff(ii,jj)=Rosenbrock_function(x);
    end
end


% Plots the level curves using the Matlab function contour
Nlevel=10;  % Number of level curves in the contour plot
LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize';
figure(1), contour(xv1,xv2,ff,Nlevel,LW,1.2), colorbar
axis([x1min x1max x2min x2max]), axis square
hold on


% Initial estimate of the minimum
x0=[-1; 1];

% Define the options to be used with the fminunc solver:
% The quasi-newton algorithm is used because it does not require the
% gradient; the default algorithm requires the gradient
options = optimoptions('fminunc','Algorithm','quasi-newton');

% Uses the solver fminunc to compute the minimum of the function defined in
% the Matlab function defined in the file Rosenbrock_function.m
xopt=fminunc(@Rosenbrock_function,x0,options)

% Computes the constrained minimum associated to the constraint 
% x(1) <= 1

Ac=[1 0];
Bc=0.5;
xoptconstr=fmincon(@Rosenbrock_function,x0,Ac,Bc)

%--------------------------------------------------------------------------
% Plots the initial point as a red circle
gg=plot(x0(1),x0(2),'or');
set(gg,'Linewidth',1.5);

% Plots the final estimate of the unconstrained minimum as a red cross
gg=plot(xopt(1),xopt(2),'xr');
set(gg,'Linewidth',1.5);

% Plots the final estimate of the constrained minimum as a red star
gg=plot(xoptconstr(1),xoptconstr(2),'*r');
set(gg,'Linewidth',1.5);

%plots the constraint boundary

%z1c=[0:0.1:1];
%z2c=[1:-0.1:0];

%gg=plot(z1c,z2c,'k');
gg = xline(0.5, 'k')
set(gg,'Linewidth',1.5);

% Identifies axis
gg=xlabel('x_1');
set(gg,'FontSize',14);

gg=ylabel('x_2');
set(gg,'FontSize',14);

% Feasible region
x_patch = [x1min, 0.5, 0.5, x1min];      % Corners X coordinates
y_patch = [x2min, x2min, x2max, x2max];  % Corners Y coordinates
pg = patch(x_patch, y_patch, [0.5 0.8 0.8]); % Gray
set(pg, 'FaceAlpha', 0.1); % Transparency (30%)
set(pg, 'EdgeColor', 'none');           
uistack(pg, 'bottom'); % Put behind level curves
title('Level Curves: Unconstrained and Constrained Minimization', 'FontSize',16)
legend('Feasible Region','Level Curves', 'Initial Estimate', 'Unconstrained minimum', 'Constrained minimum' ...
    , 'FontSize', 12)
hold off

%--------------------------------------------------------------------------
% Plots the 3d view of the function

figure(2)
surf(xx1,xx2,ff);

% Identifies axis
gg=xlabel('x_1');
set(gg,'FontSize',14);

gg=ylabel('x_2');
set(gg,'FontSize',14);

gg=zlabel('f(x)');
set(gg,'FontSize',14);
gg= title('3D Surface Plot of the Rosenbrock Objective Function ');
set(gg, 'FontSize', 14)
%--------------------------------------------------------------------------
% End of File