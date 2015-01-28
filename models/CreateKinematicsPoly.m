function [Xdot] = CreateKinematicsPoly(t,X,U)

n = 3;  % order of the Taylor expansion

couplingCoeff = 0e-1;  % small coupling coefficient to make more amenable to sos computations

x = X(1); y = X(2); th = X(3);
uV = U(1);  uW = U(2);

% Governing equations
xdot = uV*cosapprox(th,n);
ydot = uV*sinapprox(th,n);
thdot = uW;

% Introduce a tiny amount of coupling
xdot = xdot + couplingCoeff*x - couplingCoeff*y;
ydot = ydot + couplingCoeff*y - couplingCoeff*x;
thdot = thdot - couplingCoeff*th - couplingCoeff*x - couplingCoeff*y;

Xdot = [xdot; ydot; thdot];
