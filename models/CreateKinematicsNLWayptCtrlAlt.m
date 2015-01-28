function [Xdot] = CreateKinematicsNLWayptCtrlAlt(t,X)

global e gotopt xyPath closeEnough

x = X(1); y = X(2); s = X(3); c = X(4);

[U,gotopt] = CreateKinematicsNLWayptCtrl1_ctrl(X,gotopt,e,xyPath,closeEnough);

uV = U(1); uW = U(2);

% Governing equations
xdot = uV*c;
ydot = uV*s;
sdot = uW*c;
cdot = -uW*s;

Xdot = [xdot; ydot; sdot; cdot];
