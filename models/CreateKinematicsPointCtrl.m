function [Xdot] = CreateKinematicsPointCtrl(t,X)

global x_look e goalPoint

x = X(1); y = X(2); th = X(3);

xd = goalPoint(1);
yd = goalPoint(2);

Vx = xd - x;
Vy = yd - y;
[uV, uW] = feedbackLin(Vx, Vy, th, e);

% Governing equations
xdot = uV*cos(th);
ydot = uV*sin(th);
thdot = uW;

Xdot = [xdot; ydot; thdot];
