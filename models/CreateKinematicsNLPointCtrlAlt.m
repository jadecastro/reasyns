function [Xdot] = CreateKinematicsNLPointCtrlAlt(t,X)

global x_look e xyzG

x = X(1); y = X(2); s = X(3); c = X(4);

% Goal point
xd = xyzG(1);
yd = xyzG(2);

Vx = xd - x;
Vy = yd - y;

uV = 2; 
uW = 1/e*(-Vx*s + Vy*c);  % Commanded omega due to feedback linearizaion policy

% Governing equations
xdot = uV*c;
ydot = uV*s;
sdot = uW*c;
cdot = -uW*s;

Xdot = [xdot; ydot; sdot; cdot];
