function [Xdot] = VDPModelWithInput(t,X,U)

global mu

% chain of integrators

x = X(1); y = X(2);
ux = U(1); uy = U(2);

% Governing equations
xdot = y + ux;
ydot = mu*(1 - x^2)*y - x + uy;

Xdot = [xdot; ydot];
