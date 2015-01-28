function [Xdot] = CreateKinematicsNLPoly(t,X,U)

n = 3;  % order of the Taylor expansion

x = X(1); y = X(2); th = X(3);
uW = U(1);

uV = 2;

% Governing equations
xdot = uV*cosapprox(th,n);
ydot = uV*sinapprox(th,n);
thdot = uW;

Xdot = [xdot; ydot; thdot];
