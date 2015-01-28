function [Xdot] = CreateKinematicsPointCtrlPoly(t,X)

global x_look e

n = 1;  % order of the Taylor expansion

x = X(1); y = X(2); th = X(3);

% Follow a line
xd = 0;%x + x_look;
yd = -3;

Vx = xd - x;
Vy = yd - y;

% Governing equations
xdot = Vx*cosapprox(th,n)^2 + Vy*cosapprox(th,n)*sinapprox(th,n);
ydot = Vx*cosapprox(th,n)*sinapprox(th,n) + Vy*sinapprox(th,n)^2;
thdot = 1/e*(-Vx*sinapprox(th,n) + Vy*cosapprox(th,n));

Xdot = [xdot; ydot; thdot];
