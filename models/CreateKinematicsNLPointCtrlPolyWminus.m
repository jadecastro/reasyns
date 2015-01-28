function [Xdot] = CreateKinematicsNLPointCtrlPolyWminus(t,X)

global x_look e

n = 3;  % order of the Taylor expansion

x = X(1); y = X(2); th = X(3);

% Goal point
xd = -4;%x + x_look;
yd = 0;

Vx = xd - x;
Vy = yd - y;

uV = 2; % Vx*cosapprox(th,n) + Vy*sinapprox(th,n);
% uW = 1/e*(-Vx*sinapprox(th,n) + Vy*cosapprox(th,n));

% Limit on velocity
% uV = max(uV,0);  % Limits velocity on (1,0)
uW = -1;  % Limits angular velocity on (-1,1)

% Governing equations
xdot = uV*cosapprox(th,n);
ydot = uV*sinapprox(th,n);
thdot = uW;

Xdot = [xdot; ydot; thdot];
