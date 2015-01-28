function [Xdot] = CreateKinematicsNLPathCtrl(t,X)

global x_look e

x = X(1); y = X(2); th = X(3);

% % Wrap angle
% if th > pi
%     th = th - 2*pi;
% elseif th < -pi
%     th = th + 2*pi;
% end

% Follow a trajectory along a line
xd = x + x_look;
yd = 0;

Vx = xd - x;
Vy = yd - y;
[uV, uW] = feedbackLin(Vx, Vy, th, e);

% Limit on velocity
uV = 2;%max(uV,0);  % Limits velocity on (1,0)
uW = max(min(uW,1),-1);  % Limits angular velocity on (-1,1)

% Governing equations
xdot = uV*cos(th);
ydot = uV*sin(th);
thdot = uW;

Xdot = [xdot; ydot; thdot];
