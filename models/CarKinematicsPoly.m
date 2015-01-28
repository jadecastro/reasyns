function [Xdot] = CarKinematicsPoly(t,X,U)

global l

n = 3;  % order of the Taylor expansion

x = X(1); y = X(2); th = X(3); phi = X(4);
uV = U(1);  uW = U(2);

% Governing equations
xdot = uV*cosapprox(th,n);
ydot = uV*sinapprox(th,n);
thdot = uV*tanapprox(phi,n)/l;
phidot = uW;

Xdot = [xdot; ydot; thdot; phidot];
