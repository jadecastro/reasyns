function [Xdot] = CreateKinematicsPathCtrlPoly(t,X)

global x_look e

n = 1;  % order of the Taylor expansion

x = X(1); y = X(2); th = X(3);

% Follow a line

% Governing equations
xdot = x_look*cosapprox(th,n)^2 - y*cosapprox(th,n)*sinapprox(th,n);
ydot = x_look*cosapprox(th,n)*sinapprox(th,n) - y*sinapprox(th,n)^2;
thdot = -1/e*(x_look*sinapprox(th,n) + y*cosapprox(th,n));

Xdot = [xdot; ydot; thdot];
