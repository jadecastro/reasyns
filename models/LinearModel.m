function [Xdot] = LinearModel(t,X)

global xyPath
% chain of integrators

x = X(1); y = X(2);

vx = xyPath(2,1) - y;
%vy = xyPath(2,2) - y;

% Governing equations
xdot = vx;
ydot = x;

Xdot = [xdot; ydot];
