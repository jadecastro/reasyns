function [Xdot] = LinearModelWithInput(t,X,U)

% chain of integrators

x = X(1); y = X(2);
ux = U(1); uy = U(2);

% Governing equations
xdot = ux;
ydot = x + uy;

Xdot = [xdot; ydot];
