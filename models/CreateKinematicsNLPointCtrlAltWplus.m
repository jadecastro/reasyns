function [Xdot] = CreateKinematicsNLPointCtrlAltWplus(t,X)

global x_look e xyzG

x = X(1); y = X(2); s = X(3); c = X(4);

% Goal point
xd = xyzG(1);
yd = xyzG(2);

Vx = xd - x;
Vy = yd - y;

uV = 2; 

% Limit on velocity
uW = 3;  % Limits angular velocity on (-1,1)

% Governing equations
xdot = uV*c;
ydot = uV*s;
sdot = uW*c;
cdot = -uW*s;

Xdot = [xdot; ydot; sdot; cdot];
