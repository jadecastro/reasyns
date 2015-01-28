function [Xdot] = VDPModel(t,X)

global xyPath mu
% chain of integrators

x = X(1); y = X(2);

vx = 0;
vy = xyPath(2,1) - x;

% Governing equations
xdot = y + vx;
ydot = mu*(1 - x^2)*y - x + vy;

Xdot = [xdot; ydot];
