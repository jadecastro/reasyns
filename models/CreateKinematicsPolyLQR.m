function [Xdot] = CreateKinematicsPolyLQR(t,X,U)

global K

n = 3;  % order of the Taylor expansion

x = X(1); y = X(2); th = X(3);
uV = U(1); uW = U(2);

u = [uV;uW] - K*[x;y;th];
uV = u(1);
uW = u(2);

% Governing equations
xdot = uV*cosapprox(th,n);
ydot = uV*sinapprox(th,n);
thdot = uW;

Xdot = [xdot; ydot; thdot];
