function [Xdot] = CreateKinematicsPolyAlt(t,X,U)

x = X(1); y = X(2); s = X(3); c = X(4);
uV = U(1);  uW = U(2);

% Governing equations
xdot = uV*c;
ydot = uV*s;
sdot = uW*c;
cdot = -uW*s;

Xdot = [xdot; ydot; sdot; cdot];
