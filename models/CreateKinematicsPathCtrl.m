function [Xdot] = CreateKinematicsPathCtrl(t,X)

global x_look e

x = X(1); y = X(2); th = X(3);

% Follow a trajectory along a line
xd = x + x_look;
yd = 0;

Vx = xd - x;
Vy = yd - y;
[uV, uW] = feedbackLin(Vx, Vy, th, e);

% Governing equations
xdot = uV*cos(th);
ydot = uV*sin(th);
thdot = uW;

Xdot = [xdot; ydot; thdot];
