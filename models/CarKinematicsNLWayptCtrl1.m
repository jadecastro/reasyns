function [Xdot] = CarKinematicsNLWayptCtrl1(t,X)

global e gotopt xyPath closeEnough l

x = X(1); y = X(2); th = X(3); phi = X(4);

% xyPath
% gotopt

[U,gotopt] = CarKinematicsNLWayptCtrl1_ctrl(X,gotopt,l,e,xyPath,closeEnough);

uV = U(1); uW = U(2);

% Governing equations
xdot = uV*cos(th);
ydot = uV*sin(th);
thdot = uV*tan(phi)/l;
phidot = uW;

Xdot = [xdot; ydot; thdot; phidot];
